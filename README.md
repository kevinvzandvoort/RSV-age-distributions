# **Age distribution of respiratory syncytial virus disease in children \<5 years of age in low-income and middle-income countries: a systematic review and meta-analysis. Analysis code.**

This repository contains the analysis code used in our manuscript:

Sarwat Mahmud, Kevin van Zandvoort, Ling Guo, You Li, Harish Nair, Daniel R Feikin, Erin Sparrow, Fahmida Chowdhury, Cheryl Cohen, Ghassan Dbaibo, Angela Gentile, Christopher Gill, Aubrey Gordon, Katherine Horton, Cao Qing, Kirill Stolyarov, Andrew D Clark, and RSV Age Study Collaborators. *Age distribution of respiratory syncytial virus disease in children \<5 years of age in low-income and middle-income countries: a systematic review and meta-analysis.*

We are not able to make the individual study-level data used to fit age-distributions to publicly available. However, we here provide the Stan model code used to fit age-distributions to the data.

### **Models**

All models are included in the `./models/` folder, as implemented in Stan.

These are:

-   `rsv_flat_gamma.stan`, `rsv_flat_lognormal.stan`, `rsv_flat_loglogistic.stan`, `rsv_flat_burr.stan`: these are non-hierarchical (flat) models to fit a single age-distribution to a single study for a single RSV outcome. Used to compare the fit of different probability distributions.

-   `rsv_hier1_burr.stan`: this is a single-level hierarchical model in which study-specific parameters are assumed to be drawn from a global distribution (centered on hyperparameters). Used to fit age-distributions to multiple studies for a single RSV outcome, accounting for between-study variation.

-   `rsv_hier2_burr.stan`: this is a multilevel hierarchical model in which study-specific parameters are drawn from covariate-level distributions (e.g. WHO region, climate, under five mortality, and income). Each covariate-level distribution is in turn centered around global hyperparameters.

The Stan models are fitted using the `{rstan}` *R* package. The corresponding R-scripts are:

`rsv_fit_flat.R`: for fitting non-hierarchical (flat) models

`rsv_fit_hier1.R`: for fitting single-level hierarchical models

`rsv_fit_hier2.R`: for fitting multilevel hierarchical models with one coveriate
