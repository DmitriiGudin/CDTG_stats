## What is this?
This is the latest version of the implementation of the statistical method used in the following astrophysics paper: Dmitrii Gudin et al 2021 ApJ 908 79 (https://iopscience.iop.org/article/10.3847/1538-4357/abd7ed/meta).

The scientific question of interest is: giving a dataset and a set of subgroups in it, how much more similarity there is between the units within the groups, compared to the units between different groups? For instance, we may have a population of students grouped by the high school they are attending, and we may be interested in comparing the spread of grades within each school compared to the spread across the general population.

In the paper linked above, a dataset of chemically peculiar stars in the Milky Way was subjected to this question. The hypothesis was that the stars with enhanced amounts of heavy chemical elements were born in a set of environments with similar astrophysical events responsible for that enrichment. The stars were put into groups based on the similarity between their orbits around the Milky Way, and the expectation was that the stars within groups, having preferentially come from the same environments, would exhibit reduced spread of most chemical abundances.

So how do we measure this similarity? The approach used in the paper compared the distribution of measures of scale of various chemical abundances within the groups to the null distribution - the latter was estimated by performing a large number of random permutations of the full dataset, randomly assigning groups to the stars within it. Cumulative distribution functions (CDF) of each chemical abundance's scale within a group were evaluated, allowing calculation of p-values of the scale of each measured abundance in each group of stars.

Having those $p$-values, we can ask: how many of those value lay below, say, the value of 0.5? The probability of each individual one laying below 0.5 has the Bernoulli distribution with parameter 0.5, and the probability of a given number of (independent) $p$-values laying below 0.5 has the binomial distribution with parameter 0.5. Expanding this analysis, we can combine the values from multiple groups and multiple abundances to predict the overall statistical significance of the reduction of the scale of chemical abundances compared to what would be expected from random chance.

This procedure takes in a full dataset and a list of groups formed from that dataset (the groups must not have intersections; not all dataset elements have to be featured in any of the groups). It outputs a \*.txt and a TeX table files listing the following probabilities:

* *Individual Elemental Abundance Dispersion (IEAD) probability*: the measure of significance of reduction in scale for the given abundance and a given $p$-value.
* *Global Elemental Abundance Dispersion (GEAD) probability*: the measure of overall significance of reduction in scale for the given abundance.
* *Full Elemental Abundance Dispersion (FEAD) probability*: the measure of overall significance of reduction in scale for a given $p$-value.
* *Overall Elemental Abundance Dispersion (OEAD) probability*: the overall significance of results.

The lower these values are, the more statistically significance the reduction in scale is, supporting the hypothesis on the similarity of origins of stars within groups.


## Installation.
No installation required: simply put all files into a folder, and voila!


## Procedure.
1. Create a folder for your data in the data/ directory. For example, in Linux:
```bash
mkdir data/example
```
2. Put files *data.csv* and *CDTG.csv* in the created folder. The first file must contain a column of star names, a column of values for each abundance of interest, and a column of flags for the values to be excluded from analysis. The second file must contain a list of groups, each group listing star names on consequitive lines, with a \textit{text_CDTG_begin} cell at the beginning of each group, and with the \textit{text_CDTG_end} cell at the end of it (set these values in the *data/structure.py* file).
3. Run the procedure with *python run.py*, with the folder name as an attribute:
```python
python run.py example
```
Files *binom.txt* and *binom_table.tex* will be created in the folder.


## Program parameters.
The parameters of the program are set in two files. *settings.py* contains the procedural parameters (this is probably what you will want to edit), while *lib/structure.py* relates to structure of the data.
