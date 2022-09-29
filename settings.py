# This file sets the running parameters of the program. See data/structure.py for advanced settings.

# Minimum number of abundance measurements for which location and scale are estimated.
min_cluster_size = 3

# Minimum number of stars in a cluster to use biweight location/scale estimates. For smaller clusters, the mean and standard deviation are used instead.
biweight_estimator_min_cluster_size = 4

# Number of Monte-Carlo samples to simulate.
N_MC_samples = 1000

# List of binomial thresholds to use when calculating binomial probabilities, in ascending order.
binom_thresholds = [0.25, 0.33, 0.5]

# Number of significant digits in percentages for reported binomial probabilities. 3 numbers: IEAD/GEAD, FEAD, and OEAD probabilities.
binom_prec = [1, 2, 5]