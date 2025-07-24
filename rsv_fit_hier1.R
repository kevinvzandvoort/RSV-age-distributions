#' Run single-level hierarchical model
#' - one model per outcome
#' - study-level random effect
#' - covariate-level random effect

#' load packages and set Stan options to run in parallel
pacman::p_load(rstan, data.table, readxl, magrittr, posterior, bayesplot)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#' read data
#source("./process_dataset.R") #Nb. data omitted in repository
source("./distributions.R")

if(!dir.exists("./model_fits")) dir.create("model_fits")

#' outcomes in dataset
outcome_names = c("COM" = "Community cases",
                  "OPD" = "Outpatients",
                  "ERM" = "Emergency room visits",
                  "XTR" = "Inpatients (exc ICU)",
                  "IPD" = "Inpatients",
                  "ICU" = "Intensive care unit",
                  "DTH" = "Facility deaths",
                  "DTC" = "Community deaths")

for(DATASET_OUTCOME in names(outcome_names)){
  message(sprintf("Fitting DATASET_OUTCOME %s", DATASET_OUTCOME))
  
  #' let community death age-distribution be informed by facility death
  if(DATASET_OUTCOME == "DTC"){
    dataset_subset = rbind(dataset_subset,
                           data[dataset %in% characteristics[presentation_code == "DTH", dataset_code]])
  }
  dataset_subset %>% setorder(dataset, outcome, age_high, value)
  
  #' unique datasets for this outcome
  D = dataset_subset[, length(unique(dataset))]
  
  #' number of datapoints per dataset
  Ndat = dataset_subset[!is.na(value), .N, by="dataset"] %>% .[, N]
  
  #' total datapoints
  Ntotal = sum(Ndat)
  
  #' all observed number of cases
  obs = dataset_subset[!is.na(value), value]
  
  #' get upper and lower ages for every datapoint
  age_high = dataset_subset[!is.na(value), age_high]
  age_low = dataset_subset[!is.na(value), age_low]
  
  #' get max and min age for every dataset
  max_ages = dataset_subset[!is.na(value), .(maxage = max(age_high)), by="dataset"] %>% .[, maxage]
  min_ages = dataset_subset[!is.na(value), .(minage = min(age_low)), by="dataset"] %>% .[, minage]
  
  #' data to pass to STAN
  rsv_data = list(D = D,
                 Ntotal = Ntotal,
                 Ndat = Ndat,
                 max_ages = max_ages,
                 min_ages = min_ages,
                 obs = round(obs), #Nb stochastic rounding doesn't work with NUTS
                 age_low = age_low,
                 age_high = age_high)
  
  #' fit model
  #' set init_r to 0.5 so as to not sample any extreme initial values to minimize any convergence issues
  fit = stan(file = "./models/rsv_hier1_burr.stan",
             data = rsv_data,
             sample_file = sprintf("./model_fits/out_%s.RDS", DATASET_OUTCOME),
             chains = 5,
             iter = 2000,
             init = "random", init_r = 0.5)
  
  saveRDS(list(fit = fit, rsv_data = rsv_data),
          sprintf("./model_fits/fit_%s.RDS", DATASET_OUTCOME))
  
  #' asses traceplots if ran in interactive session
  #mcmc_trace(fit, pars = c("mu_k", "sigma_k", "mu_c", "sigma_c", "mu_lambda", "sigma_lambda", "k[1]", "k[2]", "lambda[1]", "lambda[2]", "c[1]", "c[2]"))
}
