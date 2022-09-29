import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import settings

dataset_csv_file = 'data.csv' # Dataset file (list of stars).
CDTGs_csv_file = 'CDTG.csv' # File containing the CDTG information.
pickle_file = 'CDTG.pickle' # File containing the Dataset and CDTGs class objects.
CDF_file = 'CDFs.hdf5' # File containing the CDFs of the abundance scales.

name_column = 'NAME' # Header of the column of star names in the datafile.

# List of tuples of abundances. Each tuple has abundance name ('[Fe/H]'), abundance codename ('FeH'), and the column name for said abundance in the datafile ('[Fe/H]').
abundances = {\
'[Fe/H]' : ('FeH','[Fe/H]'),\
'[C/Fe]' : ('CFe','[C/Fe]_Corrected'),\
'[Mg/Fe]' : ('MgFe', '[Mg/Fe]'),\
'[Sr/Fe]' : ('SrFe', '[Sr/Fe]'),\
'[Y/Fe]' : ('YFe', '[Y/Fe]'),\
'[Ba/Fe]' : ('BaFe', '[Ba/Fe]'),\
'[Eu/Fe]' : ('EuFe', '[Eu/Fe]')\
}

# Dictionary of abundances ('[Fe/H]') : tuples of abundance flag (upper limit / lower limit/ etc.) column names ('L[Fe/H]') and flag values for precise measurements ('D'). For abundances with no flag columns, use None for the latter two.
abundance_limit_flags = {\
'[Fe/H]' : (None,None),\
'[C/Fe]' : ('L[C/Fe]', 'D'),\
'[Mg/Fe]' : ('L[Mg/Fe]', 'D'),\
'[Sr/Fe]' : ('L[Sr/Fe]', 'D'),\
'[Y/Fe]' : ('L[Y/Fe]', 'D'),\
'[Ba/Fe]' : ('L[Ba/Fe]', 'D'),\
'[Eu/Fe]' : ('L[Eu/Fe]', 'D')\
}
    
# Only these abundances will be used in the binomial analysis.    
binom_abundances = ['[Fe/H]', '[C/Fe]', '[Mg/Fe]', '[Sr/Fe]', '[Y/Fe]', '[Ba/Fe]', '[Eu/Fe]']