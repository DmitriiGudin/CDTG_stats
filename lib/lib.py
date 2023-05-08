from __future__ import division
# This file contains the functions used throughout the program.


import numpy as np
import pandas as pd
import random
import pickle
import multiprocessing
import h5py
import itertools
from joblib import Parallel, delayed
from astropy.stats import biweight_scale, biweight_location
from scipy.stats import binom
import sys 
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import settings
import structure
import output
from utils import prob, remove_dict_list_NaNs




# Generates a pickle file containing the dataset (Pandas data frame) and a set of CDTGs (Pandas data frame).
def preprocess_data(dataset_name, data_file=structure.dataset_csv_file, CDTGs_file=structure.CDTGs_csv_file):
    
    data_file, CDTGs_file = 'data/'+dataset_name+'/'+data_file, 'data/'+dataset_name+'/'+CDTGs_file

    # Determine which columns in the dataset to consider.
    flag_columns = [structure.abundance_limit_flags[k][0] for k in structure.abundance_limit_flags.keys() if not (structure.abundance_limit_flags[k][0] == None)]
    dataset_columns = [structure.name_column] + list(structure.abundances.keys()) + flag_columns
    
    # Get the dataset as Pandas data frame.
    dataset = pd.read_csv(data_file, encoding='utf8', delimiter=',', index_col=structure.name_column, usecols=dataset_columns)
    dataset.index = dataset.index.map(str)
    
    # Cast data types properly. Replace missing values with NaNs.
    for c in dataset_columns:
        if c in structure.abundances.keys():
            dataset[c] = pd.to_numeric(dataset[c], errors='coerce')
            
    # Replace upper/lower limits with NaNs.
    for k in structure.abundance_limit_flags.keys():
        if structure.abundance_limit_flags[k][0] in dataset_columns and not (structure.abundance_limit_flags[k][1]==None):
            if not (structure.abundance_limit_flags[k][1]==None):
                dataset.loc[dataset[structure.abundance_limit_flags[k][0]] != structure.abundance_limit_flags[k][1],k] = np.nan
                
    # Remove now unnecessary upper/lower limit flag columns.
    dataset.drop(flag_columns, axis=1, inplace=True)
                
    
    # Obtain the list of CDTGs (list of lists of star names).
    CDTGs_df = pd.read_csv(CDTGs_file, encoding='utf8', delimiter=',', header=None, usecols=[0,0], names=['col'])
    vals = CDTGs_df['col'].values
    CDTGs = []
    star_names = []
    is_Reading = False
    for v in vals:
        if v!=structure.text_CDTG_end:
            if is_Reading:
                star_names.append(v)
            if v==structure.text_CDTG_begin:
                star_names = []
                is_Reading = True
        else:
            is_Reading = False
            CDTGs.append(star_names)
            star_names=[]
    
    # Create a Pandas dataframe for each CDTG with scale and location measurements and numbers of measurements for each abundance.
    column_list = [[k+'_location', k+'_scale', k+'_N_values'] for k in structure.abundances.keys()]
    column_list = [i for l in column_list for i in l]
    CDTGs_df = pd.DataFrame(columns=column_list)
    for C in CDTGs:
        locations, scales, numbers = {}, {}, {}
        for k in structure.abundances.keys():
            abundances = dataset[dataset.index.isin(C)][k].values
            abundances = abundances[~np.isnan(abundances)]
            numbers[k+'_N_values'] = len(abundances)
            if numbers[k+'_N_values']>=settings.min_cluster_size:
                if numbers[k+'_N_values']>=settings.biweight_estimator_min_cluster_size:
                    scales[k+'_scale'] = biweight_scale(abundances)
                    locations[k+'_location'] = biweight_location(abundances)
                else:
                    scales[k+'_scale'] = np.std(abundances)
                    locations[k+'_location'] = np.mean(abundances)
            else:
                scales[k+'_scale'] = np.nan 
                locations[k+'_location'] = np.nan
        Dict = {}
        Dict.update(locations)
        Dict.update(scales)
        Dict.update(numbers)
        CDTGs_df = pd.concat([CDTGs_df, pd.DataFrame(Dict, index=[0])], ignore_index=True)
        
    # Save the results as a pickle file.
    file = open('data/'+dataset_name+'/'+structure.pickle_file, 'wb')
    pickle.dump({'dataset':dataset, 'CDTGs':CDTGs_df}, file)
    file.close()
    
    
    

# Generates N_MC Monte Carlo samples of subsets of fixed length of the specified abundance array and calculates the scale measure of each. Returns the array of scale values. Works in parallel.
def gen_scale_distrib(abundance_array, abundance_subset_len, N_MC=1000):
    
    # First remove all NaNs. Then convert the array of abundances into a list (otherwise the random.sample() function will not work).
    abundance_array = abundance_array[~np.isnan(abundance_array)]
    abundance_array = list(abundance_array)
    
    # Set the scale function depending on the size of the subset. 
    def f_scale(arr):
        return np.std(arr)
    if abundance_subset_len>=settings.biweight_estimator_min_cluster_size:
        def f_scale(arr):
            return biweight_scale(arr)
    
    # Fill the return array (of scale values) in parallel.
    return_array = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(f_scale)(random.sample(abundance_array, abundance_subset_len)) for i in range(N_MC))

    return return_array




# Generates the dispersion distributions as required by the dataset. Saves the result as an *.hdf5 file. Pickle file containing the pre-processed data is required (CDTG CDF information is added there).
def gen_dispersion_distribs(dataset_name):
    
    # First, let us build a dictionary of CDTG sizes for each abundance we need to simulate.
    pickle_file = 'data/'+dataset_name+'/'+structure.pickle_file
    CDTG_sizes_dict = {}
    file = open(pickle_file, 'rb')
    data = pickle.load(file)
    file.close()
    dataset, CDTGs_df = data['dataset'], data['CDTGs']
    for k in structure.abundances.keys():
        N_values_list = list(set(CDTGs_df[k+'_N_values'].values))
        CDTG_sizes_dict[k] = [l for l in N_values_list if l>=settings.min_cluster_size]
        
    # Now we create the *.HDF5 file and fill it with values; the dataset names are Abundance_size, for example 'FeH_5'.
    f = h5py.File('data/'+dataset_name+'/'+structure.CDF_file,'w')
    N_generators = len([N for k in CDTG_sizes_dict.keys() for N in CDTG_sizes_dict[k]])
    i=0
    for k in CDTG_sizes_dict.keys():
        for N in CDTG_sizes_dict[k]:
            i+=1
            print ("Generating CDF "+str(i)+" out of "+str(N_generators)+" ...")
            f.create_dataset(structure.abundances[k][0]+'_'+str(N), (settings.N_MC_samples,), dtype='f')
            f[structure.abundances[k][0]+'_'+str(N)][:] = gen_scale_distrib(dataset[k].values, N, settings.N_MC_samples)
    
    # Finally, add the CDTG CDF values to the pickle file.
    CDF_dict = {}
    N_CDFs = len([scale for k in structure.abundances.keys() for scale in CDTGs_df[k+'_scale'].values if (~np.isnan(scale))])
    i=0
    for k in CDTG_sizes_dict.keys():
        CDF_list = []
        for x, N in zip(CDTGs_df[k+'_scale'].values, CDTGs_df[k+'_N_values'].values):
            if N<settings.min_cluster_size:
                CDF_list.append(np.nan)
            else:
                i+=1
                print ("Calculating CDTG CDF "+str(i)+" out of "+str(N_CDFs)+" ...")
                cdf_array = f[structure.abundances[k][0]+'_'+str(N)][:]
                cdf = len(cdf_array[cdf_array<=x])/settings.N_MC_samples
                CDF_list.append(cdf)
        CDF_dict[k] = CDF_list
    f.close()
        
    file = open(pickle_file, 'wb')
    pickle.dump({'dataset':dataset, 'CDTGs':CDTGs_df, 'CDFs':CDF_dict}, file)
    file.close()
    
    
    
    
# Calculates the IEAD probability for the list of CDFs L and the alpha value a.
def get_IEAD (L, a):
    L = np.array(L)
    N_below = len(L[L<=a])
    return 1-binom.cdf(N_below-1, len(L), a)




# Calculates the GEAD probability for the list of CDFs L and the list of alpha values a.
def get_GEAD (L, alpha_list):
    
    # Get the numbers of CDTGs below different alpha levels.
    L = np.array(L)
    IEADs = [get_IEAD(L, a) for a in alpha_list]
    N = [len(L[L<=a]) for a in alpha_list]
    
    # Get the list of combinations of N-s to go through.
    iter_lists = [np.arange(n, len(L)+1, 1) for n in N]
    iter_list = np.array(list(itertools.product(*iter_lists)))
    iter_list = np.array([l for l in iter_list if np.all(l[:-1] <= l[1:])])
    
    # Transform alpha values and numbers of CDTGs into marginal values. Transform the found above combinations accordingly.
    alpha_list_t = alpha_list[:]
    for i in range(len(alpha_list_t)-1,0,-1):
        alpha_list_t[i] -= alpha_list_t[i-1]
        N[i] -= N[i-1]
        for j in range(len(iter_list)):
            iter_list[j][i] -= iter_list[j][i-1]
            
    # Calculate the GEAD probability.
    GEAD_prob = 0
    for il in iter_list:
        probs = []
        N_left, alpha_left = len(L), 1
        for a, n in zip(alpha_list_t, il):
            probs.append(binom.pmf(n, N_left, a/alpha_left))
            N_left -= n
            alpha_left -= a
        GEAD_prob += np.prod(probs)
        
    return GEAD_prob




# Calculates the FEAD probability for the list of lists of CDFs L and the alpha value a.
def get_FEAD (L, a):
    IEAD_list = [get_IEAD(l, a) for l in L]
    return get_GEAD (sorted(IEAD_list),sorted(IEAD_list))




# Calculates the OEAD probability for the list of lists of CDFs L and the list of alpha values a.
def get_OEAD (L, alpha_list):
    return np.prod([get_GEAD(l,alpha_list) for l in L])




# Converts a fraction into percentage with specified precision.
def prob(prob, prec):
    prec_str = '.'+str(prec)+'f'
    return str(format (prob*100, prec_str))+'%'
    
    
    
    
# Performs the binomial statistical analysis. Saves the result into the pickle file. Also produces a LaTeX file for the binomial table.
def run_binom_analysis(dataset_name):
    
    # Obtain the list of CDTG CDF values from the pickle file. Remove all NaNs.
    pickle_file = 'data/'+dataset_name+'/'+structure.pickle_file
    file = open(pickle_file, 'rb')
    data = pickle.load(file)
    CDFs = remove_dict_list_NaNs(data['CDFs'])
        
    # Load other pickle data to save back later.
    dataset = data['dataset']
    CDTGs_df = data['CDTGs']
    CDF_dict = data['CDFs']    
    file.close()
        
    # Built dictionaries for numbers below, IEAD, GEAD and FEAD probabilities. Calculate the OEAD probability.
    print ('Calculating the numbers below the alpha levels...')
    num_dict = {k : {a : len(np.array(CDFs[k])[np.array(CDFs[k])<=a]) for a in settings.binom_thresholds+[1]} for k in structure.binom_abundances}
    print ('Calculating the IEAD probabilities...')
    IEAD_dict = {k : {a : get_IEAD(CDFs[k], a) for a in settings.binom_thresholds} for k in structure.binom_abundances}   
    print ('Calculating the GEAD probabilities...')
    GEAD_dict = {k : get_GEAD(CDFs[k], settings.binom_thresholds) for k in structure.binom_abundances}
    print ('Calculating the FEAD probabilities...')
    FEAD_dict = {a : get_FEAD([CDFs[k] for k in structure.binom_abundances], a) for a in settings.binom_thresholds}
    print ('Calculating the OEAD probability...')
    OEAD = get_OEAD([CDFs[k] for k in CDFs.keys()], settings.binom_thresholds)
    
    # Reorder the entries in the dictionaries as per settings.py.
    if settings.reverse_binom_threshold_order:
        num_dict = {k : {a : num_dict[k][a] for a in list(reversed(list(num_dict[k].keys())[:-1])) + [list(num_dict[k].keys())[-1]] } for k in num_dict.keys()}
        IEAD_dict = {k : {a : IEAD_dict[k][a] for a in reversed(list(IEAD_dict[k].keys()))} for k in IEAD_dict.keys()}
        FEAD_dict = {a : FEAD_dict[a] for a in reversed(list(FEAD_dict.keys()))}
    prob_dict = {'num_dict': num_dict, 'IEAD_dict' : IEAD_dict, 'GEAD_dict' : GEAD_dict, 'FEAD_dict' : FEAD_dict, 'OEAD' : OEAD}
    
    # Write the resulting probabilities into the pickle file.
    file = open(pickle_file, 'wb')
    pickle.dump({'dataset':dataset, 'CDTGs':CDTGs_df, 'CDFs':CDF_dict, 'binom_probs':prob_dict}, file)
    file.close()    
    
    # Produce the LaTeX table.
    output.get_binom_tex_table(num_dict, IEAD_dict, GEAD_dict, FEAD_dict, OEAD, dataset_name)