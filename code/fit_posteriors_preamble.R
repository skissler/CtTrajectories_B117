library(tidyverse)
library(lazyeval)
library(rstan) 
library(shinystan) 
library(purrr)
library(data.table)
options(mc.cores=parallel::detectCores())
source('code/utils.R')
source("code/set_global_pars.R")

# Store the number of people we've kept: 
n_indiv <- length(unique(ct_dat_refined$PersonID))

# Define a pared-down dataset for passing to Stan: 
indiv_data <- ct_dat_refined %>% 
	select(PersonID, PersonIDClean, TestDateIndex, CtT1, B117Status) %>%
	rename(id=PersonID) %>%
	rename(id_clean=PersonIDClean) %>% 
	rename(t=TestDateIndex) %>%
	rename(y=CtT1) %>%
	rename(b117=B117Status) %>%
	mutate(b117=case_when(b117=="Yes"~1,TRUE~0)) %>%
	trim_negatives(global_pars)

# Trim to 6 b117 and 6 non-b117 (comment for a full run):
# indiv_data <- indiv_data %>% 
# 	split(.$b117) %>%
# 	map(~ group_by(., id)) %>% 
# 	map(~ sample_n_groups(., 10)) %>% 
# 	bind_rows() %>%
# 	rename(PersonID=id) %>%
# 	select(-id_clean) %>%
# 	clean_person_id %>%
# 	rename(id=PersonID) %>% 
# 	rename(id_clean=PersonIDClean) %>%
# 	select(id, id_clean, t, y, b117) %>% 
# 	ungroup() 
# n_indiv <- length(unique(indiv_data$id))

b117 <- indiv_data %>%
	group_by(id) %>%
	slice(1) %>%
	select(id, b117) %>%
	arrange(id) %>%
	pull(b117)

# Useful dataframe for mapping official ID to Stan ID:
id_map <- indiv_data %>% 
	group_by(id) %>%
	summarise(id_clean=first(id_clean)) %>% 
	select(id, id_clean) %>%
	mutate(id_clean=as.character(id_clean))

# Useful dataframe for mapping official ID to symptoms
b117_map <- indiv_data %>% 
	group_by(id) %>%
	summarise(b117=first(b117))

prior_pars <- list(
	b117=b117,
	tpsd=2,
	dpmean_prior=global_pars[["lod"]]/2,
	dpsd_prior=global_pars[["lod"]]/6,
	wpmax=14,
	wpmean_prior=14/2,
	wpsd_prior=14/6,
	wrmax=30,
	wrmean_prior=30/2,
	wrsd_prior=30/6,
	sigma_max=10,
	sigma_prior_scale=5,
	lambda=0.01,
	fpmean=1/log(10)  # so that 90% of mass is <1 and 99% is <2
	)