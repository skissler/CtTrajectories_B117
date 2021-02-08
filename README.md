# Data and code associated with "Densely sampled viral trajectories suggest longer duration of acute infection with B.1.1.7 variant relative to non-B.1.1.7 SARS-CoV-2"
S.M. Kissler`*`, J.R. Fauver`*`, C. Mack`*`, C. Tai, M. Breban, A. Watkins, R.M. Samant, D.J. Anderson, D.D. Ho, N. Grubaugh`+`, Y.H. Grad`+`

`*` denotes equal contribution

`+` denotes co-senior author

Correspondence: skissler@hsph.harvard.edu

__run_analysis.R__ is the main analysis file. It calls all other functions in the proper order. It calls data from the data/ directory, saves key figures into the figures/ directory, and saves the fitted distributions into the output/ directory. 

__ct_dat_refined.csv__ is the main data file. It contains a person ID, a diagnosis of novel/persistent infection, the B.1.1.7 status (Yes/No), the Test Date Index (days from that person's minimum recorded Ct value), and a clean Person ID column running from 1 to 65.

__fit_posteriors.R__ runs the MCMC sampling scheme to fit the viral trajectory parameters.

__fit_posteriors_b117.R__ contains the Stan program for fitting the viral trajectory parameters.  

__fit_posteriors_preamble.R__ defines key data structures for fitting the viral trajectories.

__make_figures.R__ generates the figures from the manuscript. 

__save_figures.R__ writes the figures to .pdf. The directory for saving the figures is speified in __run_analysis.R__ (default is the figures/ directory). 

__set_global_pars.R__ defines parameters needed throughout the code (currently just the PCR limit of detection). 

__summarise_dists.R__ conveniently summarizes the mean posterior viral trajectory distributions. 

__utils.R__ contains key functions needed for the analysis. 


All code is available under the GNU General Public License, version 3, included in this repository under the following terms: 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

