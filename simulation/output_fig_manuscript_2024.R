# setwd("/Volumes/PS2000/ubuntu/rDLR_final_code/simulation")

source("./functions_fig.R")

fig.dir = "fig"
file.name="p"
ar1=0.7
num_p = 20

out1 = foo(make.pdf = FALSE, type = 1, p=num_p, s=5, num_segs = 5, ar1=ar1, rho = 0.35, case = "p20", spatial="homo", type_beta_true=1, legend.pos="none", type_plot = 5)
out2 = foo(make.pdf = FALSE, type = 1, p=num_p, s=5, num_segs = 5, ar1=ar1, rho = 0.35, case = "p20", spatial="homo", type_beta_true=2, legend.pos="none", type_plot = 5)
out3 = foo(make.pdf = FALSE, type = 1, p=num_p, s=5, num_segs = 5, ar1=ar1, rho = 0.35, case = "p20", spatial="homo", type_beta_true=3, legend.pos="none", type_plot = 5)
out4 = foo(make.pdf = FALSE, type = 1, p=num_p, s=5, num_segs = 5, ar1=ar1, rho = 0.35, case = "p20", spatial="homo", type_beta_true=4, legend.pos="bottom", type_plot = 4)

plots1 = grid.arrange(out1, out2, out3, out4, nrow = 4, ncol=1) 
ggsave(sprintf("%s/%s_s%.0f_num_segs%.0f_%s_ARone%.1f_rho%.2f_p%i_huge.pdf",fig.dir,file.name,5,25,"homo",ar1,0.35,num_p), plot = plots1,
       height=12, width=9, device="pdf")


