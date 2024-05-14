library(gamsel)

sim_myspline = function(x, response, tim = NULL, intercept = TRUE, df_smspline = 5, plot = FALSE, type_basis = "smoothspline",
                        family = "logit", verbose = FALSE,
                        tol_inner = 10^{-3}, # weighted backfitting alg
                        tol_outer = 10^(-5), # general local scoring alg
                        maxiter_k = 200, # converged in outer loop k quickly
                        maxiter_t = 50,
                        degree_truncate = 10
                        ){
  # x is a matrix
  # y is a vector
  if (class(response) == "matrix") {
    if (ncol(response) == 1) {
      y = as.numeric(response)
    } else {
      stop("The dimension of response is not correct. Plz check the response.")
    }
  } else {
    y = response
  }

  nobs = nrow(x)
  p = ncol(x)
  if (is.null(tim)) {
    tim = seq(1,nobs)
  }
  # obtain the smoothing-spline basis in terms of time
  # this does not consider the intercept, since the first column of sm_spline and first value of d is for the linear transformation
  if (type_basis == "smoothspline"){
    sm_spline = gamsel::basis.gen(x=tim, df = df_smspline, degree = degree_truncate) # degree is for the approximated basis mat of size n x degree
    # this is the same for each x_j because it is based on the time 1, 2, ..., nobs
    sm_spline_basis = sm_spline
    
    lambda_omega = attr(sm_spline_basis,"parms")$d
  } else if (type_basis == "quad") {
    sm_spline = stats::poly(x = tim, degree = 2) # help(poly), no intercept column
    sm_spline_basis = sm_spline
    lambda_omega = rep(0, 2) # naive quad basis does not have the penalty omega matrix
    }


  ybar = mean(y)
  if (intercept) {
    s0 = rep(compute_logit(ybar),nobs)
  } else {
    s0 = rep(0,nobs)
  }
  smat = matrix(0, ncol = p, nrow = nobs)
  coef_s = smat # coef_s stores the functional coefficients
  tmp_dev = 0
  dev_vec = rep(NA, maxiter_k)

  # based on the algorithm in the paper
  for (k in 1:maxiter_k) {
    # k = 1
    tmp_dev_old = tmp_dev
    if (k == 1) {
      tmp_eta = s0
    } else {
      tmp_eta = s0 + rowSums(fitted_s)
    }
    tmp_p = plogis(tmp_eta, lower.tail=TRUE) # compute_sigmoid(tmp_eta)
    tmp_w_vec = tmp_p * (1- tmp_p)
    tmp_z = tmp_eta + 1/tmp_w_vec * (y - tmp_p)
    residual = matrix(0, ncol = p, nrow = nobs)
    #if (k==1) fitted_s_old = fitted_s = smat  # init
    criterion_vec = rep(NA, maxiter_t)
    
    # below is the weighted backfitting algorithm
    for (t in 1:maxiter_t){
      # j = 1
      if (t == 1) {
        fitted_s_old = fitted_s = smat  #init as 0 for a new Weight
      } else {
        fitted_s_old = fitted_s
      }
      for  (j in 1:p){
        if (j==1) {
          residual[,1] = tmp_z - s0  - rowSums(check_mat_type(fitted_s[,c(2:p)]))
        } else if (j==p) {
          residual[,j] = tmp_z - s0 - rowSums(check_mat_type(fitted_s[,c(1:j-1)]))
        } else {
          residual[,j] = tmp_z - s0 - rowSums(check_mat_type(fitted_s[,seq(from = 1, to = j-1)])) - rowSums(check_mat_type(fitted_s_old[,seq(from=j+1, to=p)]))
        }
        xtWx = x[,j] * tmp_w_vec * x[,j]
        diag_xtWx = diag(xtWx)
        lhs = crossprod(sm_spline_basis, diag_xtWx) %*% sm_spline_basis + diag(lambda_omega)
        rhs = crossprod(sm_spline_basis, x[,j] * tmp_w_vec * residual[,j])
        rhs
        tmp_theta_single_s = solve(lhs, rhs)

        coef_s[,j] = sm_spline_basis %*% tmp_theta_single_s
        fitted_s[,j] = x[,j] * coef_s[,j]

      } # end for-loop j
      
      # check convergence
      if (t != 1) {
        numerator = sum(apply(X = fitted_s_old - fitted_s, MARGIN=2, FUN = compute_wxx, w = tmp_w_vec))  # (fitted_s_old - fitted_s)^2*tmp_w_vec
        denominator = sum(apply(X = fitted_s_old, MARGIN=2, FUN = compute_wxx, w = tmp_w_vec)) + sum(tmp_w_vec) #(fitted_s_old)^2*tmp_w_vec + tmp_w_vec
        criterion_vec[t] = sum(numerator)/sum(denominator)
        if (criterion_vec[t] < tol_inner) {
          # extract fitted_s and coef_s
          break
        }
      }
    } # for-loop t

    # check the deviance convergence
    tmp_prob = plogis(rowSums(fitted_s) + s0, lower.tail=TRUE) #  compute_sigmoid(rowSums(fitted_s) + s0)

    tmp_dev = dev_vec[k] = metrics_dev (response = y, fitted = tmp_prob)
    if (k != 1) {
      if (abs(dev_vec[k] - tmp_dev_old)/tmp_dev_old < tol_outer){
        break
      }
    }
  } # for-loop k

  auc = metrics (response = response, fitted = tmp_prob, type = "auc")
  misclass = metrics (response = response, fitted = tmp_prob, type = "misclass")

  coef_final = cbind(s0, coef_s)
  colnames(coef_final) = c("intercept", colnames(x))
  out = list(fitted = tmp_prob, coef = coef_final, tim = tim, dev = tmp_dev, auc = auc, misclass = misclass, numiter_outer = k,intercept=intercept, x = x, response = response )
  class(out) = "myspline"

  if (verbose){
    par(mfrow=c(2,2))
    for (k in 1:ncol(out$coef)){
      plot(out$coef[,k], main = paste0("column is ", k))
    }
    par(mfrow=c(1,1))
    plot(round(out$fitted), ylim = c(-0.1, 1.1), col = "blue", pch = 2); points(out$response, col = "red")
  }

  return(out)


}

track_tim = function(tim, ref_tim, ref_coeff) {
  # tim is a single integer
  # ref_tim is a numeric vector
  # ref_coef is a matrix

  tim_left = tim - 1
  tim_right = tim + 1
  if (tim %in% ref_tim){
    out = ref_coeff[which(ref_tim %in% tim),]
  } else if ((tim_left %in% ref_tim) & (tim_right %in% ref_tim)) {
    out = (ref_coeff[which(ref_tim %in% tim_left),] + ref_coeff[which(ref_tim %in% tim_right),])/2
  } else if ((tim_left %in% ref_tim) & !(tim_right %in% ref_tim)){
    out = ref_coeff[which(ref_tim %in% tim_left),]
  } else if (!(tim_left %in% ref_tim) & (tim_right %in% ref_tim)){
    out = ref_coeff[which(ref_tim %in% tim_right),]
  }
  return (out) # numeric vector
}


coef.myspline = function(object, tim_test = NULL, type = "matrix"){
  #  type = "vector" is used to be consistent with other methods

  if (is.null(tim_test)) {
    out = object$coef
  } else {
    out = t(sapply(tim_test, track_tim, ref_tim=object$tim, ref_coeff=object$coef))
  }

  if (type == "vector") {
    # such that the beta is in a vector format, for calculation in the simulation study
    out = as.vector(out)
  }
  return(out)
}


predict.myspline= function(object,newx =NULL, coef_test = NULL, tim_test = NULL, type = "probability"){

  if (is.null(newx)) {
    newx = object$x
  }
  if ("t" %in% colnames(newx)) {
    idx = which(colnames(newx) %in% "t")
    tim_test = newx[,idx]
    newx = newx[,-idx]
  }

  if (is.vector(newx)) {
    newx = as.matrix(newx, ncol=1)
  }
  if (!all(newx[,1] == 1)) {
    newx = cbind(1,newx) # intercept is always counted in gam
  }

  fitted = c()
  coeff = coef(object, tim_test)

  # row-by-row computation
  for (j in 1:nrow(newx)){
    tmp = sum(coeff[j,] * newx[j,]) #xbeta for the smoothing spline representations
    fitted[j] = plogis(tmp, lower.tail=TRUE)
  }
  out = fitted
  if (type == "response") {
    out = round(out)
  }
  return (out)
}



compute_logit = function(x){
  return (log(x/(1-x)))
}


compute_sigmoid = function(x){
  exp_x = exp(x)
  return (exp_x/(1+exp_x))
}


compute_wxx= function(x, w) {
  # x, w are both vectors
  return(x*x*w)
}


check_mat_type = function(x){
  #if (class(x) == "numeric"){
  if ("numeric" %in% class(x)) {
    out = as.matrix(x)
  } else {
    out = x
  }
  return (out)
}
