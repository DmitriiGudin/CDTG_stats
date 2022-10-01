from __future__ import division
# This file contains useful utilities that simplify the code elsewhere.


import numpy as np




# Converts a fraction into percentage with specified precision. Set LaTeX=True to use in LaTeX.
def prob(prob, prec, LaTeX=False):
    prec_str = '.'+str(prec)+'f'
    if LaTeX:
        return str(format (prob*100, prec_str))+'\\%'    
    return str(format (prob*100, prec_str))+'%'




# Takes in a dictionary with values being lists. Removes all the NaN values from those lists and returns the updated dictionary.
def remove_dict_list_NaNs(dict_list):
    dict_list = dict(dict_list)
    for k in dict_list.keys():
        cur_list = np.array(dict_list[k])
        cur_list = cur_list[~np.isnan(cur_list)]
        dict_list[k] = list(cur_list)
    return dict_list