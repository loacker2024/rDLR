
reg.funs = list()
reg.funs[["cvlasso"]] = function(x,y,type.measure="deviance") {
  out = sim_cvlasso(x,y,nlambda=50, intercept=TRUE, family = "binomial",type.measure=type.measure, nfolds=5)
  return(out)
}

reg.funs[["Mars"]] = function(x, y) {
  out = sim_cvMARS (x, y, numFolds = 5, plot = FALSE, trace = 0)
  return (out)
}

reg.funs[["VCM_fuse_grp"]] = function(x, y, intercept=FALSE, nlambda = 14, verbose=FALSE, digit_nzero = 3) {

  options(scipen = 999)
  # # both the default path_lambda_fuse and path_lambda_gr

  out = cv.vcm (x=x, response=y, family = "logit", fused_lasso = TRUE, group_lasso=TRUE
                , path_lambda_fuse = NULL, path_lambda_gr = NULL, nlambda_fuse = nlambda
                , intercept = intercept, standardize = TRUE
                , numFolds = 5, foldid = NULL, warm_start = FALSE
                , digit_nzero = digit_nzero
                , cv.type.measure = c("deviance")
                , type.measure = c("deviance") # class is the misclassification error
                , ABSTOL = 1e-4, RELTOL = 1e-2, num_iterADMM = 1e3
                , maxiter = 1, eps = 1e-16, eta = 2, cst_mu = 10, verbose=verbose)
  return (out)
}
reg.funs[["VCM_fuse"]] = function(x, y, intercept=FALSE, nlambda = 30, path_lambda_fuse = NULL, verbose=FALSE, digit_nzero = 3) {

  #Only the default path_lambda_fuse, and path_lambda_gr =0
  path_lambda_gr = 0

  out = cv.vcm (x=x, response=y, family = "logit", fused_lasso = TRUE, group_lasso=TRUE
                , path_lambda_fuse = path_lambda_fuse, path_lambda_gr = path_lambda_gr, nlambda_fuse = nlambda
                , intercept = intercept, standardize = TRUE
                , numFolds = 5, foldid = NULL, warm_start = FALSE
                , digit_nzero = digit_nzero
                , cv.type.measure = c("deviance")
                , type.measure = c("deviance") # class is the misclassification error
                , ABSTOL = 1e-4, RELTOL = 1e-2, num_iterADMM = 1e3
                , maxiter = 1, eps = 1e-16, eta = 2, cst_mu = 10, verbose=verbose)
  return (out)
}

reg.funs[["method2"]] = function(x, y, nlambda = 12, verbose=FALSE, digit_nzero=3, niter=10^5, tauStart=5) {
  out = 0
  return (out)
}


reg.funs[["basis_myspline"]] = function(t=NULL, x, y, intercept = FALSE, type_basis = "smoothspline", verbose=FALSE
                                      ) {

  out = sim_cv.myspline (t, x, y, intercept, type_basis, num_df_smspline=5
                         , plot = FALSE,
                         family = "logit", numFolds = 5, foldid = NULL, cv.type.measure = "deviance"
                         , structured_cv = TRUE
                         , trace = 0, verbose = verbose)
  return (out)
}


reg.funs[["basis_myquad"]] = function(t=NULL, x, y, intercept = FALSE, type_basis = "quad", verbose=FALSE
) {

  out = sim_cv.myquad (t, x, y, intercept, type_basis
                         , plot = FALSE,
                         family = "logit", numFolds = 5, foldid = NULL, cv.type.measure = "deviance"
                         , structured_cv = TRUE
                         , trace = 0, verbose = verbose)
  return (out)
}



get_model_size = function(object, digit_nzero=3){
  list = (apply(round(coef(object), digit_nzero), MARGIN=2, unique))
  count = 0
  for (k in 1:length(list)) {
    if (length(list[[k]]) == 1) {
      if (list[[k]] != 0) {
        count = count+1
      }
    } else {
      count = count+1
    }
  }
  out = count
  return (out)
}


###########################################
# function
# get_cvm from proposed method
# obtain the DEV and MER from the 5-fold CV
###########################################
get_cvm = function(reg.obj, name_method) {

  dev_cv_tmp = apply(reg.obj$metric_dev_raw, 2, mean, na.rm = TRUE)
  if ("mylongfused" %in% strsplit(name_method,"_")[[1]]) {
    num_obs_tmp = nrow(reg.obj$vcm.object[[1]]$x)
  } else if ("mygam" %in% strsplit(name_method,"_")[[1]]) {
    num_obs_tmp = length(reg.obj$t)
  } else if ("myspline" %in% strsplit(name_method,"_")[[1]]) {
    # names(reg.obj$myspline.object)
    num_obs_tmp = nrow(reg.obj$myspline.object$x)
  } else if ("myquad" %in% strsplit(name_method,"_")[[1]]) {
    # names(reg.obj$myspline.object)
    num_obs_tmp = nrow(reg.obj$myspline.object$x)
  } else {
    num_obs_tmp = ifelse( "grp" %in% strsplit(name_method,"_")[[1]], nrow(reg.obj$vcm.object[[1]]$z2_mat), nrow(reg.obj$vcm.object$z2_mat) )
  }
  # num_obs_tmp = ifelse( "grp" %in% strsplit(name_method,"_")[[1]], nrow(reg.obj$vcm.object[[1]]$z2_mat), nrow(reg.obj$vcm.object$z2_mat) )
  misclass_cv_tmp = apply(reg.obj$metric_misclass_raw, 2, mean, na.rm = TRUE)

  # sd
  sd_dev_cv_tmp = apply(reg.obj$metric_dev_raw, 2, sd, na.rm = TRUE)
  sd_misclass_cv_tmp = apply(reg.obj$metric_misclass_raw, 2, sd, na.rm = TRUE)

  if (("myspline" %in% strsplit(name_method,"_")[[1]]) | ("myquad" %in% strsplit(name_method,"_")[[1]])) {
    dev.cv = dev_cv_tmp
    dev.cv.sd = sd_dev_cv_tmp
    misclass.cv = misclass_cv_tmp
    misclass.cv.sd = sd_misclass_cv_tmp
  } else{
    if (sum(as.numeric(unlist(reg.obj$nzero) < num_obs_tmp))>0) {
      which_cvmin = dev_cv_tmp %in% min(dev_cv_tmp[unlist(reg.obj$nzero) < num_obs_tmp]) # obtain cvmin, after filtering the large nonzero cases
      dev.cv = dev_cv_tmp[which_cvmin][1] # choose [1] is because there could be ties
      dev.cv.sd = sd_dev_cv_tmp[which_cvmin][1]

      which_cvmin = misclass_cv_tmp %in% min(misclass_cv_tmp[unlist(reg.obj$nzero) < num_obs_tmp]) # obtain cvmin, after filtering the large nonzero cases
      misclass.cv = misclass_cv_tmp[which_cvmin][1] # choose [1] is because there could be ties
      misclass.cv.sd = sd_misclass_cv_tmp[which_cvmin][1]
    }
  }
  out = c(dev.cv, dev.cv.sd, misclass.cv, misclass.cv.sd )
  names(out) = c("dev.cv", "dev.cv.sd", "misclass.cv", "misclass.cv.sd" )
  return(out)
}




plot_coef = function(dat, j, xlab = '', geom_type = "point", ylim = NULL){
  coefs = dat
  if (colnames(coefs)[1] != "id") {
    coefs = cbind(seq(1,nrow(coefs)), coefs)
    colnames(coefs)[1] = "id"
  }
  # j is the index o column in dat, corresponding to the variable column
  #j = 2
  if (xlab !=""){
    xlab = "Time"
  }
  ylab = colnames(coefs)[j]
  min_y = min(coefs[,j])
  max_y = max(coefs[,j])
  if (is.null(ylim)) {
    ylim_left = ifelse(min_y < 0, min_y * 1.2, min_y * 0.8)
    ylim_right = ifelse(max_y < 0, max_y * 0.8, max_y * 1.2)
  } else {
    ylim_left = ylim[1]
    ylim_right = ylim[2]
  }
  # print(paste0("ylim_right is ", ylim_right, sep = " "))
  myy = coefs[,j]
  if (geom_type == "line") {
    geom_my = geom_line
  } else {
    geom_my = geom_point
  }
  gp1 = ggplot(data = as.data.frame(coefs) , aes(id)) +
    geom_my(aes(y = myy, colour = "X2"), lwd = 2) +
    xlab(xlab) + ylab(ylab) + theme_bw() +ylim(c(ylim_left, ylim_right)) +
    #theme(legend.pos=c(0.6,0.90), legend.direction = "horizontal") +
    theme(legend.pos="none") +
    theme(axis.text.x = element_text(face="bold", color="#993333",
                                     size=18, angle=45),
          axis.text.y = element_text(face="bold", color="#993333",
                                     size=18, angle=45),
          axis.title=element_text(size=18,face="bold") ) +
    geom_hline(yintercept=0, linetype="dashed", color = "grey") +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  return (gp1)
}


