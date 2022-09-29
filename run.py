from __future__ import division
# This file runs CDF generators.


import sys
import os 
sys.path.append(os.path.join(os.path.dirname(__file__), "lib"))
import lib.lib as lib
import structure

    
    
    
if __name__ == "__main__":
    
    lib.preprocess_data(sys.argv[1])
    lib.gen_dispersion_distribs(sys.argv[1])
    lib.run_binom_analysis(sys.argv[1])