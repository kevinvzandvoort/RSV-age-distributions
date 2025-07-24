#' Run two-level hierarchical model
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

outcome_names = c("COM" = "Community cases",
                  "OPD" = "Outpatients",
                  "ERM" = "Emergency room visits",
                  "XTR" = "Inpatients (exc ICU)",
                  "IPD" = "Inpatients",
                  "ICU" = "Intensive care unit",
                  "DTH" = "Facility deaths",
                  "DTC" = "Community deaths")
colnames(characteristics)

#' covariates
variables = c("region", "climate", "income", "u5mr_grouped")

#' analysis 
DATASET_OUTCOME = names(outcome_names)[5] #Nb this is only done for Inpatients (outcome with most datapoints)
variable = variables[1]
for(variable in variables){
  message(sprintf("Doing DATASET_OUTCOME %s - %s", DATASET_OUTCOME, variable))
  
  #' ensure data is ordered correctly by dataset and outcome
  dataset_subset = data[outcome == DATASET_OUTCOME]
  
  #' let community death age-distribution be informed by facility death
  if(DATASET_OUTCOME == "DTC"){
    dataset_subset = rbind(dataset_subset,
                           data[dataset %in% characteristics[presentation_code == "DTH" & complete_u5 == "yes", dataset_code]])
  }
  dataset_subset %>% setorder(dataset, outcome, age_high, value)
  
  variable_values = dataset_subset[, .(dataset = unique(dataset))] %>%
    merge(characteristics[, c("dataset_code", variable), with = FALSE] %>%
            setNames(c("dataset_code", "variable")),
          by.x = "dataset", by.y = "dataset_code") %>%
    .[, dataset := factor(dataset, levels(dataset_subset$dataset))] %>%
    .[order(dataset), factor(variable)]
  
  variable_labels = levels(variable_values)
  variable_values = variable_values %>% as.numeric()
  
  #unique datasets for this outcome
  D = dataset_subset[, length(unique(dataset))]
  
  #unique variable values
  V = length(variable_labels)
  
  #number of datapoints per dataset
  Ndat = dataset_subset[!is.na(value), .N, by="dataset"] %>% .[, N]
  
  #variable values per dataset
  var_val = variable_values
  
  #total datapoints
  Ntotal = sum(Ndat)
  
  #all observed number of cases
  obs = dataset_subset[!is.na(value), value]
  
  #get upper and lower ages for every datapoint
  age_high = dataset_subset[!is.na(value), age_high]
  age_low = dataset_subset[!is.na(value), age_low]
  
  #get max and min age for every dataset
  max_ages = dataset_subset[!is.na(value), .(maxage = max(age_high)), by="dataset"] %>% .[, maxage]
  min_ages = dataset_subset[!is.na(value), .(minage = min(age_low)), by="dataset"] %>% .[, minage]
  
  #data to pass to STAN
  rsv_data = list(D = D,
                  V = V,
                  Ntotal = Ntotal,
                  Ndat = Ndat,
                  var_val = var_val,
                  max_ages = max_ages,
                  min_ages = min_ages,
                  obs = round(obs), #Nb stochastic rounding doesn't work with NUTS
                  age_low = age_low,
                  age_high = age_high)
  
  additional_data = list(variable_values = variable_labels,
                         dataset_values = dataset_subset[, unique(dataset)])
  
  #fit model
  #' set init_r to 0.5 to not sample any extreme initial values to prevent flat traceplots
  fit = stan(file = "./models/rsv_hier2_burr3.stan",
             data = rsv_data,
             chains = 5,
             sample_file = sprintf("./model_fits/rsv_%s_by_%s.RDS", DATASET_OUTCOME, variable),
             iter = 2000, 
             init = "random", init_r = 0.5)
  
  saveRDS(list(fit = fit, rsv_data = rsv_data, additional_data = additional_data),
          sprintf("./model_fits/pooled_fit_%s_by_%s.RDS", DATASET_OUTCOME, variable))
  
  #' asses traceplots if ran in interactive session
  #trace_plot = mcmc_trace(fit, pars = c("mu_k", "sigma_k", "mu_c", "sigma_c", "mu_lambda", "sigma_lambda", "k[1]", "k[2]", "lambda[1]", "lambda[2]", "c[1]", "c[2]"))
}
