library(tidyverse) 
library(scales)

source('code/utils.R')
source("code/set_global_pars.R")

B117inds <- c(1368,1371,1374,1375,3229,4399,4447)

for(indexZ in 1:length(B117inds)){

	ct_dat_refined <- read_csv("data/ct_dat_refined.csv")

	ct_dat_refined <- ct_dat_refined %>% 
		filter(PersonID!=B117inds[indexZ]) %>%
		select(-PersonIDClean) %>%
		clean_person_id

	source("code/fit_posteriors_preamble.R")
	source("code/fit_posteriors.R")

	source("code/make_figures.R")
	source("code/summarise_dists.R")
	savedir <- paste0("figures/no",B117inds[indexZ],"/")
	source("code/save_figures.R")
	write.csv(dist_summary, file=paste0("output/dist_summary_no",B117inds[indexZ],".csv"), row.names=FALSE)
	save(ct_fit, file=paste0("output/ct_fit_no",B117inds[indexZ],".RData"))

}

