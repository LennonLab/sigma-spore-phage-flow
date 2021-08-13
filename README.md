# spore-phage-sigma-flow
Analysis of flow cytometry data relating to aporulation with phage sigma factors

The Pipeline to analyze flow cytometry data of SYBR green stained cultures of Bacillus subtilis is modified from the work of:
Karava M, Bracharz F, Kabisch J (2019) Quantification and isolation of Bacillus subtilis spores using cell sorting and automated gating. PLoS ONE 14(7): e0219892. https://doi.org/10.1371/journal.pone.0219892
https://gitlab.com/sporesort/


Quantification relies on sample volume measurement by NovCyte machine.


To start analysis:
1. place the flow data, exprtedfrom NovoExpress as fcs files, in the data/FCM folder.
Experimental data (strain, treatment, sample dilution) are encoded in the names of folders and files.

2. Setting up the environment
* use R version 3.6.3
https://support.rstudio.com/hc/en-us/articles/200486138-Changing-R-versions-for-RStudio-desktop

* After switching the R version, open spore-phage-sigma-flow.Rproj in Rsudio

* This pipeline has many dependencies.
I have used the "renv" package to track all the dependencies.
learn more at: https://rstudio.github.io/renv/articles/renv.html

First time you will need to install the packages by running the following commands (within project)
library(renv)
renv::restore

* I also use the "here" package to call files using relative paths.
learn more at: https://here.r-lib.org/


3. Running the analysis

* from within Rstudio execute the src/quant-spore_veg.R code. Packages and functions used for this analysis are at src/FCM_functions.R.

*Results*
for each group of replicates the pipeline should make 4 plots in figs/gate_plots folder:
> singlet: FSC height vs. area. used to exclude doublets
> noise: distribution if FSC are values, with line marking the noise filtering threshold.
> scatterNoise: events scatter plot on SYBR fluorescence vs SSC area with events filtered out by noise marked in gray, and events passed to clustering in blue.
> cluster: Events scatter plot on SYBR fluorescence vs FSC area with assignment to spore and veg clusters given by color. Blue dots represent the model based cluster center for veg and spore.

The actual counts. along with pipeline values collected along the way, are saved to data/output as .csv files.


 



