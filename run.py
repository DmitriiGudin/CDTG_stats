from __future__ import division
# This file runs CDF generators.


import sys
import os 
sys.path.append(os.path.join(os.path.dirname(__file__), "lib"))
import lib.lib as lib
import structure

    
    
    
if __name__ == "__main__":
    
    if sys.argv[1] == '--all':
        dataset_list = [d for d in os.listdir('data') if d[0]!='.']
    else:
        dataset_list = sys.argv[1:]
    
    for d in dataset_list:
        print ("Working on dataset " + d + "...")
        lib.preprocess_data(d)
        lib.gen_dispersion_distribs(d)
        lib.run_binom_analysis(d)