from __future__ import division
# This file contains procedures for human-readable output generation.


import sys 
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import settings
import structure
from utils import prob




# Given a prob_dict dictionary (see lib.py), outputs the binomial probability results and saves them as a *.txt-file.
def gen_binom_txt(prob_dict):
    return None
    
    
    
    
# Given the probability dictionaries + the OEAD value (see lib.py), outputs the binomial probability results and saves them as a *.tex-table.
def get_binom_tex_table(num_dict, IEAD_dict, GEAD_dict, FEAD_dict, OEAD, dataset_name):
    
    # Adds a padding (\phantom{1}) to a number if it is below 10, otherwise returns the original number. The output is a string.
    def pad(N):
        is_percentage = True if str(N)[-1]=='%' else False
        if is_percentage:
            N = N[:-2]
            if float(N)<10:
                return '\\phantom{1}'+N+'\\%'
            else:
                return N+'\\%'
        else:
            if N<10:
                return '\\phantom{1}'+str(N)
            else:
                return str(N)
        
    # Order the binomial thresholds properly for display, as per settings.py.
    binom_thresholds = list(reversed(settings.binom_thresholds)) if settings.reverse_binom_threshold_order else settings.binom_thresholds[:]
        
    # Generate the table LaTeX code.
    f = open('data/'+dataset_name+'/'+structure.binom_table_file,'w')
    f.write('\\begin{deluxetable*}{cccccc}[!]\n')
    f.write('\\tabletypesize{\\footnotesize}\n')
    f.write('\\tablecaption{CDTG Elemental-Abundance Statistics\\label{tab:binomial_probability_'+dataset_name+'}}\n')
    alpha_list = ''.join(['$'+str(a)+'$, ' for a in binom_thresholds])[1:-2]
    f.write('\\tablehead{\\colhead{Abundance} & \\colhead{$\\#$ CDTGs} & \\colhead{$  N < ' + alpha_list + '} & \\colhead{IEAD Probabilities} & \\colhead{GEAD Probabilities} & \\colhead{OEAD Probability}}\n')
    f.write('\\startdata\n')
    for i, k in enumerate(structure.binom_abundances):
        string = '\\text{' + k + '} & '
        num_list = pad(list(num_dict[k].values())[-1]) + ' & ' + ''.join(['$'+pad(N)+'$, ' for N in list(num_dict[k].values())[:-1]])[:-2]+'& '
        percentage_list = ''.join(['$'+pad(prob(p,settings.binom_prec[0],LaTeX=True))+'$, ' for p in IEAD_dict[k].values()])
        percentage_list += ' & $' + pad(prob(GEAD_dict[k],settings.binom_prec[0],LaTeX=True)) + '$ & '
        if i==0:
            percentage_list += '\\hfil\\multirow{' + str(len(structure.binom_abundances)) + '}{*}{$' + prob(OEAD,settings.binom_prec[2],LaTeX=True) + '$}\\hfill'
        percentage_list += '\\\\'      
        string = string + num_list + percentage_list  +'\n'
        f.write(string)
    f.write('\\hline\n')
    percentage_list = ''.join(['$'+prob(FEAD_dict[k],settings.binom_prec[1],LaTeX=True)+'$, ' for k in binom_thresholds])[:-2]
    f.write('\\multicolumn{2}{c}{FEAD Probabilities} & & '+percentage_list+' & &\n')
    f.write('\\enddata\n')
    if len(structure.binom_abundances) == 1:
        alpha_list = 'level $v$ = ' + str(binom_thresholds[0])
    elif len(structure.binom_abundances) == 2:
        alpha_list = 'levels $v$ = ' + str(binom_thresholds[0]) + ' and ' + str(binom_thresholds[1]) + ', respectively'
    else:
        alpha_list = ''.join([str(a)+', ' for a in binom_thresholds[:-1]]) + 'and ' + str(binom_thresholds[-1]) + ', respectively'
    f.write('\\tablecomments{The Individual Elemental-Abundance Dispersion (IEAD) probabilities represent the binomial probabilities for each element for the ' + alpha_list + '. The Full Elemental-Abundance Dispersion (FEAD) probabilities represent the probabilities (across {\\it all} elements) for the ' + alpha_list + '. The Global Elemental-Abundance Dispersion (GEAD) probabilities represent the probabilities for the triplet of CDF levels for each element. The Overall Elemental-Abundance Dispersion (OEAD) probability represents the probability (across {\\it all} elements) resulting from random draws from the full CDF.  See text for details.}\n')
    f.write('\\end{deluxetable*}')
    f.close()