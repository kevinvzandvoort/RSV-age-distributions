#' Run flat (non-hierarchical) model
#' - one model per study per outcome

pacman::p_load(rstan, data.table, readxl, magrittr, posterior, bayesplot, actuar, ggh4x)

#' STAN options - run in parallel
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#' get dataset
#source("./process_dataset.R") #Nb. data omitted in repository
source("./distributions.R")

#' Datasets and outcomes to run flat models for
outcomes = list(
  "COHEN" = c("ICU", "OPD"),
  "DBAIBO" = c("ERM", "ICU", "OPD"),
  "QING" = c("ERM", "ICU", "OPD"))

#' evaluate ages at
eval_ages = seq(0, 5*52, 0.1)

#' get stan implementation of models (singlelevel)
mod_burr = stan_model("./models/rsv_flat_burr.stan")
mod_gamma = stan_model("./models/rsv_flat_gamma.stan")
mod_llogis = stan_model("./models/rsv_flat_loglogistic.stan")
mod_lnorm = stan_model("./models/rsv_flat_lognormal.stan")

out = lapply(names(outcomes), function(DATASET_CONTACT){
  lapply(outcomes[[DATASET_CONTACT]], function(DATASET_OUTCOME){
    message(sprintf("Doing DATASET CONTACT: %s; DATASET OUTCOME %s",
                    DATASET_CONTACT, DATASET_OUTCOME))
    
    llvals = data.table(outcome = DATASET_OUTCOME,
                        distribution = c("Burr", "Gamma", "Log Normal", "Log-logistic"),
                        value = NA_real_, ll = NA_real_)
    
    dataset_subset = data[outcome == DATASET_OUTCOME & type == DATASET_TYPE & dataset_contact == DATASET_CONTACT]
    dataset_subset %>% setorder(dataset, outcome, age_high, value)
    
    dataset = dataset_subset
    
    study_ids = dataset[, unique(dataset)]
    dataset_value_cdf = dataset %>% copy() %>%
      .[!is.na(value)] %>%
      .[, study_id := as.numeric(factor(dataset, study_ids))] %>%
      .[, .(age=age_high,
            y = cumsum(value)/sum(value),
            type = "cdf"), by=c("dataset", "study_id")]
    
    dataset_value_pdf = dataset %>% copy() %>%
      .[!is.na(value)] %>%
      .[, study_id := as.numeric(factor(dataset, study_ids))] %>%
      .[, .(age=age_high,
            y = (value)/sum(value),
            type = "pdf"), by=c("dataset", "study_id")]
    dataset_value = rbind(dataset_value_cdf, dataset_value_pdf)
    
    Ndat = dataset_subset[!is.na(value), .N, by="dataset"] %>% .[, N]
    obs = dataset_subset[!is.na(value), value]
    age_high = dataset_subset[!is.na(value), age_high]
    age_low = dataset_subset[!is.na(value), age_low]
    max_age = dataset_subset[!is.na(value), .(maxage = max(age_high)), by="dataset"] %>% .[, maxage]
    min_age = dataset_subset[!is.na(value), .(minage = min(age_low)), by="dataset"] %>% .[, minage]
    
    rsv_dat = list(Ndat = Ndat,
                   max_age = max_age,
                   min_age = min_age,
                   obs = round(obs), #stochastic rounding doesnt work with NUTS
                   age_low = age_low,
                   age_high = age_high)
    
    rcode = 1
    i = 1
    while(rcode != 0){
      if(i %% 10 == 0) message(sprintf("In Burr, i: %s", i)) 
      optfit_burr = optimizing(mod_burr, data = rsv_dat)
      rcode = optfit_burr["return_code"]
      i = i+1
    }
    llvals[distribution == "Burr", c("value", "ll") := .(optfit_burr$value,
                                                         optfit_burr$theta_tilde %>% (function(x) sum(x[grepl(pattern = "log_lik", colnames(x))])))]
    
    y_pdf = cburr_truncated(x = eval_ages, optfit_burr$par["k"], optfit_burr$par["c"], optfit_burr$par["lambda"], max=5*52)
    y_cdf = pburr_truncated(x = eval_ages, optfit_burr$par["k"], optfit_burr$par["c"], optfit_burr$par["lambda"], max=5*52)
    y_cdf_compare = pburr_truncated(x = age_high, optfit_burr$par["k"], optfit_burr$par["c"], optfit_burr$par["lambda"], max=5*52)
    out_burr = data.table(age = eval_ages, y = y_pdf, type = "pdf") %>%
      rbind(data.table(age = eval_ages, y = y_cdf, type = "cdf")) %>%
      rbind(data.table(age = age_high, y = y_cdf_compare, type = "cdf_compare")) %>%
      .[, distribution := "Burr"] %>% .[]
    
    rcode = 1
    i = 1
    while(rcode != 0){
      if(i %% 10 == 0) message(sprintf("In Gamma, i: %s", i)) 
      optfit_gamma = optimizing(mod_gamma, data = rsv_dat)
      rcode = optfit_gamma["return_code"]
      i = i+1
    }
    llvals[distribution == "Gamma", c("value", "ll") := .(optfit_gamma$value,
                                                          optfit_gamma$theta_tilde %>% (function(x) sum(x[grepl(pattern = "log_lik", colnames(x))])))]
    
    y_pdf = cgamma_truncated(x = eval_ages, optfit_gamma$par["alpha"], optfit_gamma$par["beta"], max=5*52)
    y_cdf = pgamma_truncated(x = eval_ages, optfit_gamma$par["alpha"], optfit_gamma$par["beta"], max=5*52)
    y_cdf_compare = pgamma_truncated(x = age_high, optfit_gamma$par["alpha"], optfit_gamma$par["beta"], max=5*52)
    out_gamma = data.table(age = eval_ages, y = y_pdf, type = "pdf") %>%
      rbind(data.table(age = eval_ages, y = y_cdf, type = "cdf")) %>%
      rbind(data.table(age = age_high, y = y_cdf_compare, type = "cdf_compare")) %>%
      .[, distribution := "Gamma"] %>% .[]
    
    rcode = 1
    i = 1
    while(rcode != 0){
      if(i %% 10 == 0) message(sprintf("In LLogis, i: %s", i)) 
      optfit_loglogistic = optimizing(mod_llogis, data = rsv_dat)
      rcode = optfit_loglogistic["return_code"]
      i = i+1
    }
    llvals[distribution == "Log-logistic", c("value", "ll") := .(optfit_loglogistic$value,
                                                                 optfit_loglogistic$theta_tilde %>% (function(x) sum(x[grepl(pattern = "log_lik", colnames(x))])))]
    
    y_pdf = cllogis_truncated(x = eval_ages, optfit_loglogistic$par["alpha"], optfit_loglogistic$par["beta"], max=5*52)
    y_cdf = pllogis_truncated(x = eval_ages, optfit_loglogistic$par["alpha"], optfit_loglogistic$par["beta"], max=5*52)
    y_cdf_compare = pllogis_truncated(x = age_high, optfit_loglogistic$par["alpha"], optfit_loglogistic$par["beta"], max=5*52)
    out_llogis = data.table(age = eval_ages, y = y_pdf, type = "pdf") %>%
      rbind(data.table(age = eval_ages, y = y_cdf, type = "cdf")) %>%
      rbind(data.table(age = age_high, y = y_cdf_compare, type = "cdf_compare")) %>%
      .[, distribution := "Log-logistic"] %>% .[]
    
    rcode = 1
    i = 1
    while(rcode != 0){
      if(i %% 10 == 0) message(sprintf("In lnorm, i: %s", i)) 
      optfit_lognormal = optimizing(mod_lnorm, data = rsv_dat)
      rcode = optfit_lognormal["return_code"]
      i = i+1
    }
    llvals[distribution == "Log Normal", c("value", "ll") := .(optfit_lognormal$value,
                                                               optfit_lognormal$theta_tilde %>% (function(x) sum(x[grepl(pattern = "log_lik", colnames(x))])))]
    
    y_pdf = clnorm_truncated(x = eval_ages, optfit_lognormal$par["mu"], optfit_lognormal$par["sigma"], max=5*52)
    y_cdf = plnorm_truncated(x = eval_ages, optfit_lognormal$par["mu"], optfit_lognormal$par["sigma"], max=5*52)
    y_cdf_compare = plnorm_truncated(x = age_high, optfit_lognormal$par["mu"], optfit_lognormal$par["sigma"], max=5*52)
    out_lnorm = data.table(age = eval_ages, y = y_pdf, type = "pdf") %>%
      rbind(data.table(age = eval_ages, y = y_cdf, type = "cdf")) %>%
      rbind(data.table(age = age_high, y = y_cdf_compare, type = "cdf_compare")) %>%
      .[, distribution := "Log Normal"] %>% .[]
    
    out = rbind(out_lnorm, out_llogis, out_gamma, out_burr)
    out[, c("dataset", "outcome") := .(DATASET_CONTACT, DATASET_OUTCOME)]
    llvals[, c("dataset", "outcome") := .(DATASET_CONTACT, DATASET_OUTCOME)]
    return(list(out=out, llvals=llvals))
  })
})

llvals = out %>% lapply(function(x) lapply(x, "[[", "llvals") %>% rbindlist()) %>% rbindlist()
llvals_target = llvals = llvals %>% dcast(dataset+outcome~distribution, value.var="ll")

fwrite(llvals, "ll_values.csv")
out = out %>% lapply(function(x) lapply(x, "[[", "out") %>% rbindlist()) %>% rbindlist()

dataset_values = lapply(names(outcomes), function(DATASET_CONTACT){
  lapply(outcomes[[DATASET_CONTACT]], function(DATASET_OUTCOME){
    message(sprintf("Doing DATASET CONTACT: %s; DATASET OUTCOME %s", DATASET_CONTACT, DATASET_OUTCOME))
    dataset_subset = data[outcome == DATASET_OUTCOME & type == DATASET_TYPE & dataset_contact == DATASET_CONTACT]
    dataset_subset %>% setorder(dataset, outcome, age_high, value)
    dataset = dataset_subset
    study_ids = dataset[, unique(dataset)]
    dataset_value_cdf = dataset %>% copy() %>%
      .[!is.na(value)] %>%
      .[, study_id := as.numeric(factor(dataset, study_ids))] %>%
      .[, .(age=age_high,
            y = cumsum(value)/sum(value),
            type = "cdf"), by=c("dataset", "study_id")]
    
    dataset_value_pdf = dataset_value_cdf %>% copy %>%
      .[, .(y=diff(y)/diff(age), type = "pdf"), by=c("dataset", "study_id")]
    dataset_value_pdf[, "age"] = dataset_value_cdf[-.N, age] + dataset_value_cdf[, diff(age)]
    
    dataset_value = rbind(dataset_value_cdf, dataset_value_pdf) %>%
      .[, c("dataset", "outcome") := .(DATASET_CONTACT, DATASET_OUTCOME)]
    return(dataset_value)
  }) %>% rbindlist()
}) %>% rbindlist()