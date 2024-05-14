#setwd("/Volumes/PS2000/ubuntu/rDLR_final_code/simulation")

source("./functions_tab.R")

what_vec = c("dev.cv","misclass.cv")
rho_vec = c(0.35) 
ar1_vec = c(0.7) 
p_vec = c(20) 
num_segs_vec = c(5)
type_beta_true_vec = c(1,2,3,4)
s=5

for (what in what_vec) {
  for (rho in rho_vec) {
    for (ar1 in ar1_vec) {
      tab=c()
      for (p in p_vec) {
        for (num_segs in num_segs_vec) {
          tab1 = c()
          for (type_beta_true in type_beta_true_vec) {

            if (FALSE) {
              p = 20
              ar1 = 0.7
              rho=0.35
              s=5
              num_segs = 5
              spatial = "homo"
              what = "dev.cv"
              type_beta_true = 1
            }
            
            tmp_tab = get_tab(p, s, num_segs, ar1, type_beta_true, rho, what=what)
            #print(tmp_tab)
            tab = rbind(tab, tmp_tab)
          }

        }

      } # end for-loop p_vec
      tab
      saved_file = paste0("./tab/", "tab_", what, "_pvec_", paste0(p_vec, collapse = "_"),"_ar1vec_", paste0(ar1, collapse = "_"),"_rhovec_", paste0(rho, collapse = "_"),"_",".rds")
      print(saved_file)
      saveRDS(tab, saved_file)
    }
  } # end for-loop rho_vec
}

dim(tab)
head(tab)



