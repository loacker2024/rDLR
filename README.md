install the R package ./rDLR_final_code/rDLR_0.1.0.tar.gz

go to the folder ./rDLR_final_code/simulation/

# to generate the illustration plot 

run the file ./output_fig_illustration_manuscript_2024.R

input: ./illustration_fig_dat.rds
output: ./illustration_simulated_coefficients_2024.pdf 

# to generate the tables in simulation 

run the file ./rDLR_final_code/simulation/output_tab_manuscript_2024.R

input: ./rds/*.rds
output: 
./tab/tab_dev.cv_pvec_20_ar1vec_0.7_rhovec_0.35_.rds 
./tab/tab_misclass.cv_pvec_20_ar1vec_0.7_rhovec_0.35_.rds 


# to generate the figure in simulation 

run the file ./rDLR_final_code/simulation/output_fig_manuscript_2024.R

input: ./rds/*.rds
output: ./fig/p_s5_num_segs25_homo_ARone0.7_rho0.35_p20_huge.pdf 

--- 

then, go to the folder ./rDLR_final_code/real_case/

# to generate the table and figure in real case study

run the file ./2024_real_case_hk_code.R 

input: ./real_case_hk_full_2024.csv
output: ./plot_coefs_rDLR_hk.pdf and ./tab_case_study_hk_2024.RData




