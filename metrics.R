
metrics = function(response, fitted, type = "mse"){
  if (type == "mse"){
    out = mean((response - fitted)^2)
  }
  if (type == "auc") {
    # http://www.ritchieng.com/machine-learning-evaluate-classification-model/
    if (drop(length(unique(response))) == 2) {
      # https://stackoverflow.com/questions/4903092/calculate-auc-in-r
      # https://cran.r-project.org/web/packages/ROCR/ROCR.pdf
      # https://stackoverflow.com/questions/4903092/calculate-auc-in-r
      auroc <- function(response, fitted) {
        score = fitted
        bool = response
        n1 <- sum(!bool)
        n2 <- sum(bool)
        U  <- sum(rank(score)[!bool]) - n1 * (n1 + 1) / 2
        return(1 - U / n1 / n2)
      }
      auc = auroc(response, fitted)
    } else {
      print("in this fold, the class of test y is not equal to 2, so auc is set as NA ")
      auc = NA
    }
    out = auc
  }

  ## confusion matrix
  fitted = ifelse(fitted < 0.5, 0 , 1) # if fitted <= 0.5, then precision could be NaN b/c both TP, FP are zero
  tab = table(response, fitted)

  if (length(unique(fitted)) > 1) {
    TN = tab[1,1]
    FP = tab[1,2]
    FN = tab[2,1]
    TP = tab[2,2]
  } # end if

  if (length(unique(fitted)) == 1) {
    # when the fitting is only have unique value
    if (unique(fitted) == 0) {
    TN = tab[1,1]
    FP = 0
    FN = tab[2,1]
    TP = 0
    }
    if (unique(fitted) == 1) {
      TN = 0
      FP = tab[1,1]
      FN = 0
      TP = tab[2,1]
    }
  } # end if


  if (type == "misclass") {
  #  https://stat.ethz.ch/pipermail/r-help/2011-September/288885.html
    out = 1-sum(diag(tab))/sum(tab)
  }

  if (type == "TPR") {
    out = TP/(TP + FN)
  }
  if (type == "TNR") {
    out = TN/(TN + FP)
  }
  if (type == "specificity") {
    out = FP/(TN + FP)
  }

  if (type == "precision") {
    out = TP/(TP + FP)
  }
  if (type == "f1") {
    out = 2*TP/(2*TP + FP + FN)
  }

  return(out)
}


metrics_coef = function(true_coef, esti_coef, esti_coef_base = NULL, type = "PM"){
  if (type == "PM"){
    # Dynamic quality process model in consideration of equipment degradation
    if (class(esti_coef) != "numeric") esti_coef = as.vector(esti_coef) # convert matrix to a vector by stacking column wise
    if (class(true_coef) != "numeric") true_coef = as.vector(true_coef) # convert matrix to a vector by stacking column wise
    PM = sum(abs(esti_coef - true_coef))/sum(abs(true_coef))
    out = PM
  }

  if (type == "REE") {
    # true_coef is matrix
    # esti_coef is matrix, from the proposed method
    # esti_coef_base is matrix, from the baseline method
    out = sum(abs(esti_coef - true_coef))/sum(abs(esti_coef_base - true_coef)) #* 100
  }

  if (type == "corr_zero") {
    id_0_true_coef = which(true_coef == 0,arr.ind = T)
    id_0_esti_coef = which(esti_coef == 0,arr.ind = T)

    out1 = match(data.frame(t(id_0_true_coef)), data.frame(t(id_0_esti_coef)))

    id_n0_true_coef = which(true_coef != 0,arr.ind = T)
    id_n0_esti_coef = which(esti_coef != 0,arr.ind = T)

    out2 = match(data.frame(t(id_n0_true_coef)), data.frame(t(id_n0_esti_coef))) #!is.na means correctly indentized zero or non-zero

    out = (sum(!is.na(out1))+sum(!is.na(out2)))/(nrow(esti_coef) * ncol(esti_coef))

  }
  return(out)
}


metrics_dev = function(response, fitted, Dev2 = FALSE, eps = .Machine$double.eps) {
  # https://stats.stackexchange.com/questions/108995/interpreting-residual-and-null-deviance-in-glm-r
  if (class(response) != "matrix") response = as.matrix(response, ncol=1)
  if (class(fitted) != "matrix") fitted = as.matrix(fitted, ncol=1)

  p_temp = fitted
  prob = ifelse(p_temp==0, eps,
                ifelse(p_temp ==1, p_temp - eps,p_temp))
  nll_model = -base::colSums( log(prob) * response + (1-response) * log(1-prob) )


  # null deviance
  p_temp = rep(mean(response), times = nrow(response))
  prob = ifelse(p_temp==0, eps,
                ifelse(p_temp ==1, p_temp - eps,p_temp))
if (class(prob) != "matrix") prob = matrix(prob,ncol=1)
  nll_null = -base::colSums( log(prob) * response + (1-response) * log(1-prob) )

  # saturated deviance
  p_temp = as.vector(response) # response is matrix
  prob = ifelse(p_temp==0, eps,
                ifelse(p_temp ==1, p_temp - eps,p_temp))
if (class(prob) != "matrix") prob = matrix(prob,ncol=1)
  nll_sat = -base::colSums( log(prob) * response + (1-response) * log(1-prob) )

  # fraction of deviance explained
  Dev_null = 2*(nll_null - nll_sat)
  Dev_model = 2*(nll_model - nll_sat)
  Dev2 = 1 - Dev_model/Dev_null

  if (Dev2 == TRUE) {
	out = Dev2
} else {
	out = Dev_model/nrow(response) # normalize
}
  return(out)
}
