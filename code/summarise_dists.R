dist_summary <- shared_params_df %>% 
	summarise(
		peak.ct.WT_mean=mean(global_pars[["lod"]]-dpmeanW),
		peak.ct.WT_lwr95=quantile(global_pars[["lod"]]-dpmeanW,0.05),
		peak.ct.WT_upr95=quantile(global_pars[["lod"]]-dpmeanW,0.95),
		proliferation.time.WT_mean=mean(wpmeanW),
		proliferation.time.WT_lwr95=quantile(wpmeanW,0.05),
		proliferation.time.WT_upr95=quantile(wpmeanW,0.95),
		clearance.time.WT_mean=mean(wrmeanW),
		clearance.time.WT_lwr95=quantile(wrmeanW,0.05),
		clearance.time.WT_upr95=quantile(wrmeanW,0.95),
		total.duration.WT_mean=mean(wpmeanW+wrmeanW),
		total.duration.WT_lwr95=quantile(wpmeanW+wrmeanW,0.05),
		total.duration.WT_upr95=quantile(wpmeanW+wrmeanW,0.95),
		peak.ct.B117_mean=mean(global_pars[["lod"]]-dpmeanB),
		peak.ct.B117_lwr95=quantile(global_pars[["lod"]]-dpmeanB,0.05),
		peak.ct.B117_upr95=quantile(global_pars[["lod"]]-dpmeanB,0.95),
		proliferation.time.B117_mean=mean(wpmeanB),
		proliferation.time.B117_lwr95=quantile(wpmeanB,0.05),
		proliferation.time.B117_upr95=quantile(wpmeanB,0.95),
		clearance.time.B117_mean=mean(wrmeanB),
		clearance.time.B117_lwr95=quantile(wrmeanB,0.05),
		clearance.time.B117_upr95=quantile(wrmeanB,0.95),
		total.duration.B117_mean=mean(wpmeanB+wrmeanB),
		total.duration.B117_lwr95=quantile(wpmeanB+wrmeanB,0.05),
		total.duration.B117_upr95=quantile(wpmeanB+wrmeanB,0.95)
		) %>%
	pivot_longer(everything()) %>%
	separate(name, c("parameter", "statistic"), sep="_") %>%
	pivot_wider(names_from=statistic, values_from=value) %>%
	arrange(parameter)


