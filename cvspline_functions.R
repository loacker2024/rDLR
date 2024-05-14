#' sim_cv.myspline
sim_cv.myspline = function(t=NULL, x, y, intercept = TRUE, num_df_smspline=7
                           #, id_sig_vars = NULL, knots = 10
                           , plot = FALSE, type_basis = "smoothspline",
                           family = "logit", numFolds = 5, foldid = NULL, cv.type.measure = "deviance"
                        , structured_cv = TRUE
                        , trace = 0, verbose = TRUE){

  if (class(x) != "matrix") x = as.matrix(x)
  if (is.null(t)) {
    message("the argument t in the functional coefficient is not provided. Using the default 1:nrow(x).")
    t = 1:nrow(x)
  }
  flag_break = FALSE

  # 5-fold cv
  if (is.null(foldid)) {
    if (structured_cv == TRUE) {
      nrow(x)
      foldid = c(0,rep(seq(1,numFolds),nrow(x)-2)[seq(1,nrow(x)-2)],0) # GAM can do both structured_cv or random_cv, use this to be consistent with proposed DLR
    } else {
      foldid = sample(rep(seq(numFolds), length = nrow(x))) # the functional coefficient is a function of t
    }
  } else {
    numFolds = max(foldid)
  }
  if (numFolds < 3)
    stop("numFolds must be bigger than 3; numFolds=5 is recommended")

  count_unique = function(x){
    length(unique(x))
  }
  while(1){ # adjust the correct num_df_smspline
    if (num_df_smspline < 3) {
      num_df_smspline = 3
      break
    }
    if (num_df_smspline > min(apply(x, MARGIN = 2, count_unique))){
      break
    }
    num_df_smspline = num_df_smspline - 1
  }
  df_smspline_vec = round(seq(3,(2+num_df_smspline),length.out=num_df_smspline))
  metric_misclass_raw_mat = metric_dev_raw_mat = c()

  for (kk in 1:length(df_smspline_vec)){ # this for-loop can be used to tune other hyper-parameters
    df_tmp = df_smspline_vec[kk]
    metric_list = metric_misclass_list = metric_dev_list = pred_test = outlist = as.list(seq(numFolds))

    for (i in seq(numFolds)) {
      which = foldid == i
      not_which = !which
      if (is.matrix(y)) {
        y_sub = y[not_which, ]
      } else {
        y_sub = y[not_which]
      }
      x_sub = x[not_which, ] # training set
      t_sub = t[not_which]

      # can do random cv, as long as t is given as a covariate for the splines
      possibleError = tryCatch({
        # outlist_tmp = sim_mygam(t=t_sub, x=x_sub,
        #                         y=y_sub, id_sig_vars = id_sig_vars,
        #                         knots = knots_tmp, plot = plot)
        outlist_tmp = sim_myspline (x=x_sub, response = y_sub, tim = t_sub, intercept = intercept, df_smspline = df_tmp,
                                    type_basis = type_basis,
                                    family = "logit", verbose = F)
      }, error = function(err) {
        if (verbose) {
          cat(paste("    Oops! The for-loop for the knots selection breaks, see message below",
                    "below; ...\n"))
          cat("    ***** Error message *****\n")
          #cat(paste0("    When n is ",n," and the number of knots is ", knots_tmp, ",\n"))
          cat(sprintf("    %s\n",err$message))
          cat("    *** End error message ***\n")
          err
          # saveRDS(reg.obj, file=paste0("./",names(reg.funs)[j],"_sim.n",n,".p",p,".s",s, ".segs", num_segs,".",spatial, ".ARone",ar1,sprintf(".rho%0.2f", rho),".nrep",i,"_reg_obj.rds"))
        }

      })
      # https://stackoverflow.com/questions/8093914/use-trycatch-skip-to-next-value-of-loop-upon-error
      if("error" %in% class(possibleError)) {
        flag_break = TRUE
        break
      }

      newx_sub = data.frame(t[which], x[which, ]) # x[which, ] is test subset
      colnames(newx_sub)[1] = c("t") # include the time t as the first column in the newx_sub

      coef_test = coef.myspline (outlist_tmp, tim_test = newx_sub[,which(colnames(newx_sub) %in% "t")])
      pred_test = predict(object = outlist_tmp, newx = newx_sub)

      # contributions_terms = predict( outlist_tmp, newdata = newx_sub, type = "terms") # contributions from each terms in gam
      
      # compute dev.cv and misclass.cv
      if (is.matrix(y)) {
        y_sub_test = y[which, ]
      } else {
        y_sub_test = y[which]
      }

      metric_dev_list[[i]] = sapply(data.frame(pred_test), metrics_dev, response=y_sub_test)
      metric_misclass_list[[i]] = sapply(data.frame(pred_test), metrics, response=y_sub_test, type = "misclass")

    } # end for-loop i

    if (flag_break) break
    metric_dev_raw = do.call(rbind, metric_dev_list)
    metric_misclass_raw = do.call(rbind, metric_misclass_list)

    metric_dev_raw_mat = cbind(metric_dev_raw_mat, metric_dev_raw)
    metric_misclass_raw_mat = cbind(metric_misclass_raw_mat, metric_misclass_raw)

  } # end for-loop kk
  if (is.null(metric_dev_raw_mat) ) {
    out = list(df_smspline_vec=NA, cvm=NA, cvsd=NA, cvup=NA, cvlo=NA, cvmin=NA, which_cvmin=NA
               , knots_cvmin=NA, nzero=NA, metric_dev_raw=NA, metric_misclass_raw=NA, mygam.object=NA,
               message="Model has more coefficients than data even at knots = 3.")
    class(out) = c("cv.mygam")
  } else {
    if (flag_break) {
      knots_end = kk-1
    } else {
      knots_end = kk
    }
    if (!is.matrix(metric_dev_raw_mat)) {
      metric_dev_raw_mat = matrix(metric_dev_raw_mat, ncol=1)
    }
    if (!is.matrix(metric_misclass_raw_mat)) {
      metric_misclass_raw_mat = matrix(metric_misclass_raw_mat, ncol=1)
    }
    colnames(metric_dev_raw_mat) = colnames(metric_misclass_raw_mat) = paste("df", df_smspline_vec, sep="=") #[1:knots_end] #seq(1,knots_end) # cause I removed the knots parameter vector

    if (cv.type.measure == "deviance") {
      metric_list = metric_dev_raw_mat #[,which_cvmin]
    } else if (cv.type.measure == "misclass") {
      metric_list = metric_misclass_raw_mat #[,which_cvmin]
    }

    cvraw = metric_list #do.call(rbind, metric_list)
    cvm  <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- apply(cvraw, 2, sd, na.rm = TRUE)
    cvup = lapply(1:length(cvm), function(m) cvm[[m]] + cvsd[[m]])
    cvlo = lapply(1:length(cvm), function(m) cvm[[m]] - cvsd[[m]])

    cvup = do.call(c, cvup)
    cvlo = do.call(c, cvlo)

    # obtain the min dev and misclass among cv under different criterion
    metric_dev_raw = metric_dev_raw_mat[,which.min(apply(metric_dev_raw_mat, MARGIN=2, mean))]
    metric_misclass_raw = metric_misclass_raw_mat[,which.min(apply(metric_misclass_raw_mat, MARGIN=2, mean))]
    cvmin = min(cvm)
    which_cvmin = which.min(cvm)
    df_cvmin = df_smspline_vec[which_cvmin]

    myspline.object = sim_myspline(x=x, response=y, intercept = intercept, df_smspline = df_cvmin, type_basis = type_basis, family = family, verbose = verbose)

    if (class(metric_dev_raw) == "numeric") metric_dev_raw = matrix(metric_dev_raw, ncol=1)
    if (class(metric_misclass_raw) == "numeric") metric_misclass_raw = matrix(metric_misclass_raw, ncol=1)

    cvm_dev = min(apply(metric_dev_raw_mat, 2, mean, na.rm = TRUE))
    cvm_misclass = min(apply(metric_misclass_raw_mat, 2, mean, na.rm = TRUE))

    out = enlist(cvm, cvsd, cvup, cvlo, cvmin, which_cvmin
                 , df_cvmin 
                 , cvm_misclass, cvm_dev
                 , metric_dev_raw, metric_misclass_raw
                 , myspline.object)

    class(out) = c("cv.myspline", class(myspline.object))

  }

  return (out)
}

#' coef function for cv.myspline object.
coef.cv.myspline = function(object, type = "matrix"){

  obj = object$myspline.object
  out = coef(object = obj,  type=type)
  return(out)
}

#' predict function for cv.myspline object.
predict.cv.myspline = function(object, newx = NULL){

  obj = object$myspline.object
  out = predict(obj, newx)
  return(out)
}


