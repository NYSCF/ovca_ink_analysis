# Analysis of ovarian cancer organoid and IPSC-derived Natural Killer cell co-cultures

Tumor heterogeneity drives resistance to therapeutics in cancers including ovarian cancer. In this repository, we use Bayesian modelling to investigate resistant and susceptible tumor populations from co-cultures of ovarian cancer organoids and IPSC-derived Natural Kill cells (iNKs). 

## The apoptosis and organoid size datasets

* `organoid_sizes.csv.gz`
* `organoid_apoptosis_pd01.csv.gz`
* `organoid_apoptosis_pd02.csv.gz`

## Bayesian modelling and organoid size analysis 

* `fit_betamix_and_plot_results.R`
* `model_sigmoid_apoptotic_by_time.R`
* `ovarian_size_correlation_analysis.R`

## Generating figures of results

All plots resulting from the analysis can be reproduced by executing the scripts in this repository.

```
cd scripts/
Rscript fit_betamix_and_plot_results.R
Rscript model_sigmoid_apoptotic_by_time.R
Rscript ovarian_size_correlation_analysis.R
```

Data and fitted models are read from `data` but code for model fitting is still included (though commented out) in the scripts. Figures are written to `figures/`.

## Citation

Please cite the study related to this work when using the data or code in this repository. Publication of a bioRxiv preprint is pending.
