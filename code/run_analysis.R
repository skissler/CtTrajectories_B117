library(tidyverse) 
library(scales)

source('code/utils.R')
source("code/set_global_pars.R")

ct_dat_refined <- read_csv("data/ct_dat_refined.csv")

# , col_types=list(
# 	col_integer(),
# 	col_factor(),
# 	col_factor(),
# 	col_double(),
# 	col_double(),
# 	col_integer()))

source("code/fit_posteriors_preamble.R")
source("code/fit_posteriors.R")
source("code/make_figures.R")
source("code/summarise_dists.R")

savedir <- paste0("figures/")
source("code/save_figures.R")

save(ct_fit, file="output/ct_fit.RData")