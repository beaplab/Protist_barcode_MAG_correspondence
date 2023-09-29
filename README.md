# SMAG-V9OTU_crspd
Online Supplement to "On the relationship between protist metabarcoding and protist metagenome-assembled genomes"

# R_scripts

## calculation :
  + sim_calc.R - calculate pairwise correspondence between the features in simulated datasets
  + sim_data_gen.R - generate simulated datasets that mimic SMAG or V9 data
  + taxogroup_distr.R - generate summary statistics on SMAG and V9 taxonomy data
  + Evil_Shuffle.R - shuffling of SMAG and V9 datasets with subsequent correspondence calculation
  + zero_abundance_SMAGs_list.R , station_intersect_counts.R - minor summary stats on SMAG and V9 abundance datasets
  + script_summary_propr_calc.R , script_summary_shuffle_stats.R , script_summary_sim_stats.R - generate summary tables from correspondence estimates of SMAG and V9 data, shuffled SMAG and V9 data and simulated data, correspondingly

## downstream : 
  + breakpoint_detection.Rmd - generate breakpoint_summary_rho_cor.tsv
	+ scenario_assignation_for_candidate_pairs.Rmd - assign scenarios to the manually curated list of SMAG-V9 OTU pairs

### 	 plots : 
 + 	 scripts used to produce plots
		 
# data
## data_inputs :
 + taxa_counts.tsv - input file with SMAG and V9 taxonomic subset sizes
 + avg_abund_V9.tsv , V9_var_vector.tsv - average abundance counts and variance of abundance counts for V9 data, used by sim_data_gen.R as input parameters for metaSPARSim simulation
 + TAGs-18S-V4_NAME-PROJ-MAT-READ-CAB_nico.list - custom list of sample names for SMAG and V9 data, used to match V9 samples with SMAG samples
 + Online_data_available_elsewhere - list with download and reference links previously published and used in this study available online 
 + final_taxa_ext_genus.csv - taxonomy dictionary of SMAGs and V9

## calc_outputs :
 + all_scores.zip - correspondence for all SMAG - V9 OTU pairs in all taxonomic subset for all metrics
 + rand_all_scores_relabund_inv_ext_shuff_scratch.zip , all_scores_relabund_simulation_matr.zip - example outputs produced by Evil_Shuffle.R and sim_calc.R, correspondingly (Note they are not exhaustive, but just an instance of the output produced!)
 + to_plot_\*_shuffle1.tsv - summary stats used to generate plots
###  summary_calcs: 
 +   summary_propr_calc_stats\*.tsv , summary_shuffle_stats\*.tsv , summary_sim_stats_\*.tsv - summary of outputs produced by script_summary_propr_calc.R , script_summary_shuffle_stats.R script_summary_sim_stats.R, correspondingly
 +   se_summary_shuffle_stats\*.tsv , se_summary_sim_stats\*.tsv - - summary of se of estimates from the multiple replicates of outputs produced by script_summary_shuffle_stats.R script_summary_sim_stats.R, correspondingly
  
##  downstream_outputs :
 +  breakpoint_summary_rho_cor.tsv - results of automatic filtering by criteria applied for breakpoint detection
 + Candidate_pairs-for_R_script_scatterplot_corrected_Pairs_automated.csv - manually curated list of SMAG-V9 OTU pairs
 + Scenarios_auto_assigned_candidate_pairs.csv - manually curated list of SMAG-V9 OTU pairs with automatically assigned scenarios
 + blastSAGs.tsv - blast output for reference genome- SAG control
 + plots.zip - various plots, along with miscellaneous data for plot generation.

# manuscript_figs_n_code
 + svg of main figures
 + source LaTeX data for the manuscript pdf
 
   
  
