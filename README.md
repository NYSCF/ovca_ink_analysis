# Analysis of apoptosis in co-cultures of ovarian cancer organoids and IPSC-derived Natural Killer cells

Tumor heterogeneity drives resistance to therapeutics in cancers including ovarian cancer. In this repository, we use Bayesian modelling to investigate resistant and susceptible tumor populations from co-cultures of patient-derived ovarian cancer organoids and IPSC-derived Natural Kill cells (iNKs). For further details please see the preprint cited below. 
## The apoptosis and organoid size datasets

This repository contains three comma-separated tables of processed imaging data in the `data/` directory.

* `organoid_sizes.csv.gz`: This table contains one measurement of an organoid per row with columns describing experiment ID `expt`, `treatment` of either control or iNK co-culture, measurement `timepoint` in 15min increments encoded as 1-39, surface `area` of all organoids in squared micrometres, apoptotic surface area of all organoids `apop_area` in squared micrometres, `roundness` of organoids defined as `sqrt(4*pi*area)/perimetre)`, and `diametre in micrometres`
* `organoid_apoptosis_pd01.csv.gz`: This wide format table contains the percentage of apoptotic surface area in treatments and controls across different time points for patient PDO1. Each value is the measurement of a single organoid and each column contains independent measurements from the other columns. 
* `organoid_apoptosis_pd02.csv.gz`: This wide format table contains the percentage of apoptotic surface area in treatments and controls across different time points for patient PDO2 in the same format as above.

## Bayesian modelling and organoid size analysis 

To decipher differences between resistant and susceptible cancer organoids, we fit models to the apoptosis and organoid size data using a Bayesian modelling framework. Three custom scripts were used in the analysis: 

* `fit_betamix_and_plot_results.R`: Read in organoid apoptosis tables and fit a Beta mixture model to the observations at the experimental end point. 
* `model_sigmoid_apoptotic_by_time.R`: Read in organoid apoptosis tables and fit multiple sigmoid models to the observed time series, distinguishing resistant and susceptible organoids.
* `ovarian_size_correlation_analysis.R`: Read in organoid size table and visualise the correlation of organoid size and apoptosis, also fitting a Gaussian mixture model to observed organoid sizes at the experimental endpoint. 

All plots resulting from the analysis can be reproduced by executing the scripts in this repository.

```
cd scripts/
Rscript fit_betamix_and_plot_results.R
Rscript model_sigmoid_apoptotic_by_time.R
Rscript ovarian_size_correlation_analysis.R
```

Pre-fitted models are read from `data/models` but code for model fitting is still included (though commented out) in the scripts. Figures are written to `figures/`.

## Citation

When using the data or code in this repository, please cite:

> Marisa Mercadante, Armin Scheben, Jacob Estrada, Jan Savas-Carstens, William Sullivan, Nicholas Housel, Tatiana Volpari, Jax Hebner, Maria Sapar, Tom Rusielewicz, Frederick J. Monsma Jr., Stefan Semrau, Yinan Wang, Laura A. Martin. 
> A patient-derived ovarian cancer organoid platform to study susceptibility to natural killer cells. *bioRxiv*, 2025.03.06.641285, doi: [https://doi.org/10.1101/2025.03.06.641285](https://doi.org/10.1101/2025.03.06.641285)
