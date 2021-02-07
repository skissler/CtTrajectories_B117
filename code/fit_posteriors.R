ct_model <- stan_model("code/fit_posteriors_b117.stan") 

fit_startq <- Sys.time()
ct_fit <- sampling(ct_model, 
	data=list(
		N=nrow(indiv_data), 
		n_id=length(unique(indiv_data$id_clean)),
		lod=global_pars[["lod"]], 
		id=indiv_data$id_clean,
		b117=as.list(prior_pars)$b117,
		t=indiv_data$t, 
		y=indiv_data$y, 
		tpsd=as.list(prior_pars)$tpsd,
		dpmean_prior=as.list(prior_pars)$dpmean_prior,
		dpsd_prior=as.list(prior_pars)$dpsd_prior,
		wpmax=as.list(prior_pars)$wpmax,
		wpmean_prior=as.list(prior_pars)$wpmean_prior,
		wpsd_prior=as.list(prior_pars)$wpsd_prior,
		wrmax=as.list(prior_pars)$wrmax,
		wrmean_prior=as.list(prior_pars)$wrmean_prior,
		wrsd_prior=as.list(prior_pars)$wrsd_prior,
		sigma_max=as.list(prior_pars)$sigma_max,
		sigma_prior_scale=as.list(prior_pars)$sigma_prior_scale,
		lambda=as.list(prior_pars)$lambda,
		fpmean=as.list(prior_pars)$fpmean), 
	iter=1000, chains=4)
# , control = list(adapt_delta=0.85)
# , control = list(adapt_delta=0.99)
# control = list(adapt_delta=0.95, max_treedepth=15)
fit_endq <- Sys.time()
print(paste0("Fit time: ",difftime(fit_endq, fit_startq, units="min")," mins"))

# launch_shinystan_nonblocking(ct_fit)

params <- rstan::extract(ct_fit)
indiv_params_df <- make_indiv_params_df(params, c("tp","dp","wp","wr"), n_indiv) %>% 
	rename(id_clean=id) %>% 
	left_join(id_map, by="id_clean") %>%
	left_join(b117_map, by="id")

shared_params_df <- make_shared_params_df(params, c("dpmeanW","wpmeanW","wrmeanW","dpmeanB","wpmeanB","wrmeanB","dpsd","wpsd","wrsd")) 

params_df <- indiv_params_df %>% 
	left_join(shared_params_df, by="iteration") %>% 
	select(-iteration) 

