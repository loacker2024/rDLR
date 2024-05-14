library(ggplot2)
library(rDLR)
###############################################################
#folder = "/Volumes/PS2000/ubuntu/rDLR_final_code/"

source(paste0(folder,"real_case/functions_reg.funs.R"))
source(paste0(folder,"common.R"))
source(paste0(folder,"cvlasso_functions.R"))
source(paste0(folder,"cvmars_functions.R"))
source(paste0(folder,"cvspline_functions.R"))
source(paste0(folder,"spline_functions.R"))
source(paste0(folder,"cvquad_functions.R"))
source(paste0(folder,"metrics.R"))

# load csv
start_time = Sys.time()
print(paste0("start time is ", start_time))
hk_full = read.csv(file = paste0(folder,"real_case/real_case_hk_full_2024.csv"))
head(hk_full)

##############################################################################################
# run methodsr
##############################################################################################
out_cvlasso_hk_misclass = reg.funs[[1]](x=hk_full[,-1], y = hk_full[,1],type.measure="class")
out_cvlasso_hk = reg.funs[[1]](x=hk_full[,-1], y = hk_full[,1])

out_VCM_fuse_grp_hk = reg.funs[[3]](x=hk_full[,-1], y = hk_full[,1], nlambda = 18, verbose = FALSE)

tmp_path_lambda_fuse = round(exp(seq(log(1), log(0.00005),
                                     length.out = 30)), digits = 6)
out_VCM_fuse_hk = reg.funs[[4]](x=hk_full[,-1], y = hk_full[,1], path_lambda_fuse = tmp_path_lambda_fuse, verbose = FALSE)

out_myspline_hk = reg.funs[[6]](x=hk_full[,-1], y = hk_full[,1])
out_myquad_hk = reg.funs[[7]](x=hk_full[,-1], y = hk_full[,1])

##############################################################################################
# plot of coefficient
##############################################################################################
# rDLR
coefs = coef(out_VCM_fuse_grp_hk)
head(coefs)

gp = gridExtra::grid.arrange(plot_coef(dat = coefs, 2, ylim = c(-1.5, 1.5))+ theme(plot.margin=margin(5.5,5.5,-20,5.5, "pt")),
                             plot_coef(dat = coefs, 3, ylim = c(-1.5, 1.5))+ theme(plot.margin=margin(5.5,5.5,-20,5.5, "pt")),
                             plot_coef(dat = coefs, 4, ylim = c(-1.5, 1.5))+ theme(plot.margin=margin(5.5,5.5,-20,5.5, "pt")),
                             plot_coef(dat = coefs, 5, ylim = c(-1.5, 1.5))+ theme(plot.margin=margin(5.5,5.5,-20,5.5, "pt")),
                             plot_coef(dat = coefs, 6, ylim = c(-1.5, 1.5))+ theme(plot.margin=margin(5.5,5.5,-20,5.5, "pt")),
                             plot_coef(dat = coefs, 7, ylim = c(-1.5, 1.5))+ theme(plot.margin=margin(5.5,5.5,-20,5.5, "pt")),
                             plot_coef(dat = coefs, 8, ylim = c(-1.5, 1.5))+ theme(plot.margin=margin(5.5,5.5,-20,5.5, "pt")),
                             plot_coef(dat = coefs, 9, ylim = c(-1.5, 1.5))+ theme(plot.margin=margin(5.5,5.5,-20,5.5, "pt")),
                             plot_coef(dat = coefs, 10, ylim = c(-1.5, 1.5))+ theme(plot.margin=margin(5.5,5.5,-20,5.5, "pt")),
                             plot_coef(dat = coefs, 11, ylim = c(-1.5, 1.5))+ theme(plot.margin=margin(5.5,5.5,-20,5.5, "pt")),
                             plot_coef(dat = coefs, 12, ylim = c(-1.5, 1.5), xlab = 'Time'),
                             plot_coef(dat = coefs, 13, ylim = c(-1.5, 1.5), xlab = 'Time'),
                             plot_coef(dat = coefs, 14, xlab = 'Time', ylim = c(-1.5, 1.5)),
                             plot_coef(dat = coefs, 15, xlab = 'Time', ylim = c(-1.5, 1.5)),
                             #plot_coef(dat = coefs, 16, xlab = 'Time'),
                             plot_coef(dat = hk_full, 2, ylim = c(-1.5, 1.5), xlab = 'Time', geom_type = "point"),
                             nrow = 3, ncol=5)

fig.dir = folder 
file.name = "real_case/plot_coefs_rDLR_hk"
h = 7
w = 12
ggsave(sprintf("%s%s.pdf",fig.dir,file.name), plot = gp,
       height=h, width=w, device="pdf")

##############################################################################################
# table performance
##############################################################################################
dev_cvmin_hk = rbind(get_cvm(reg.obj=out_VCM_fuse_grp_hk, name_method="VCM_fuse_grp"),
                     get_cvm(reg.obj=out_VCM_fuse_hk, name_method="VCM_fuse" ),
                     get_cvm(reg.obj=out_myspline_hk, name_method="myspline"),
                     get_cvm(reg.obj=out_myquad_hk, name_method="myquad"),
                     c(out_cvlasso_hk$cvm_min, out_cvlasso_hk$cvsd[which.min(out_cvlasso_hk$cvm)],
                       out_cvlasso_hk_misclass$cvm_min,out_cvlasso_hk_misclass$cvsd[which.min(out_cvlasso_hk_misclass$cvm)])
)

rownames(dev_cvmin_hk) = c("rDLR", "FUSE","VCM1","VCM2","LASSO")
dev_cvmin_hk


size = c(get_model_size(object = out_VCM_fuse_grp_hk, digit_nzero = 3),
         get_model_size(object = out_VCM_fuse_hk, digit_nzero = 3),
         get_model_size(object = out_myspline_hk, digit_nzero = 3),
         get_model_size(object = out_myquad_hk, digit_nzero = 3),
         sum(as.numeric(coef(out_cvlasso_hk)!=0))-1
)
dev_cvmin_hk = cbind(dev_cvmin_hk,size)
dev_cvmin_hk

##############################################################################################
# save results
##############################################################################################
elapsed_time = Sys.time() - start_time
print(elapsed_time)

folder
save(dev_cvmin_hk, file = paste0(folder,"real_case/tab_case_study_hk_2024.RData"))

#


