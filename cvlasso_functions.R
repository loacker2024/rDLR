#' reference
#' https://github.com/ryantibs/best-subset/blob/master/bestsubset/R/lasso.R
sim_cvlasso = function(x, y, alpha=1, nrelax=1, nlambda=50,
                 lambda.min.ratio=ifelse(nrow(x)<ncol(x),0.01,0.0001), family = "binomial",
                 type.measure="deviance", nfolds=5,
                 lambda=NULL, intercept=TRUE, standardize=TRUE) {

  # Check for glmnet package
  if (!require("glmnet",quietly=TRUE)) {
    stop("Package glmnet not installed (required here)!")
  }

  # Set up data
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)

  # Set dfmax manually
  dfmax = p
  if (nrelax > 1 && n < (p-intercept)) dfmax = n-intercept
  # TH: this sometimes still returns too many!

  # Reset nlambda if a specific lambda sequence is passed
  if (!is.null(lambda)) nlambda = length(lambda)

  # Run glmnet
  obj = cv.glmnet(x, y, alpha=alpha, nlambda=nlambda, dfmax=dfmax,family = family,
               lambda.min.ratio=lambda.min.ratio, lambda=lambda,
               type.measure=type.measure, nfolds=nfolds,
               intercept=intercept, standardize=standardize)

  # Append a few things to the returned object
  obj$nrelax = nrelax
  obj$nlambda = nlambda
  obj$intercept = intercept
  obj$x = x; obj$y = y

  obj$cvm_min = obj$cvm[obj$lambda %in% obj$lambda.min]

  class(obj) = "cv.lasso"
  return(obj)
}

#' coef function for cv.lasso object.
coef.cv.lasso = function(object, s=NULL) {

  class(object) = "cv.glmnet"
  if (is.null(s)) s = "lambda.min"
  out = glmnet::coef.cv.glmnet(object,s=s)
  out = as.matrix(out)
  return(out)
}

#' predict function for cv.lasso object.
predict.cv.lasso = function(object, newx, s=NULL, label = FALSE) {

  if (missing(newx)) newx = object$x
  if (object$intercept) {
    if (!all(newx[,1] == 1)) {
      newx = cbind(rep(1,nrow(newx)),newx)
    }
  }
  # https://stackoverflow.com/questions/39184991/why-is-predict-glmnet-not-predicting-probabilities
  xbeta = newx %*% coef.cv.lasso(object,s) # for the coef, it is a matrix, each column corresponds to a lambda value
  out = exp(xbeta)/(1+exp(xbeta)) # prob
  if (label) out = ifelse(out <= 0.5, 0 , 1)
  return(out)
}


#' Coef function for lasso object.
coef.lasso = function(object, s=NULL, gamma=NULL) {
  beta.lasso = coef.lasso.from.glmnet(object,s)
  if (object$nrelax == 1 && is.null(gamma)) {
    if (object$intercept) return(beta.lasso)
    else return(beta.lasso[-1,])
  }
  if (is.null(gamma)) gamma = seq(1,0,length=object$nrelax)

  beta.ls = coef.ls(beta.lasso,object$x,object$y)
  beta.left = matrix(apply(beta.lasso,2,function(b){b%o%gamma}),
                     nrow=nrow(beta.lasso))
  beta.right = matrix(apply(beta.ls,2,function(b){b%o%(1-gamma)}),
                      nrow=nrow(beta.lasso))
  beta.mat = beta.left + beta.right

  if (object$intercept) return(beta.mat)
  else return(beta.mat[-1,])
}

coef.lasso.from.glmnet = function(object, s=NULL) {
  class(object) = "glmnet"
  if (length(object$lambda)==object$nlambda) {
    return(glmnet::coef.glmnet(object,s=s))
  }
  else {
    min.lam = min(object$lambda)
    max.lam = max(object$lambda)
    svec = exp(seq(log(max.lam),log(min.lam),length=object$nlambda))
    return(glmnet::coef.glmnet(object,s=svec))
    ## RJT TODO: should we used exact=TRUE above? Requires additional
    ## arguments to match the initial call to glmnet(), kind of clunky
    ## TH: use glmnet.control(fdev=0) at beginning of session
    ## Still needed though for cases when df exceeds p (can happen with
    ## glmnet, and bad for relaxed lasso)
  }
}

coef.ls = function(beta, x, y) {
  n = nrow(x); p = ncol(x)
  apply(beta, 2, function(b) {
    act.set = which(b[-1] != 0)
    intercept = b[1]!=0
    if (length(act.set)==0) return(c(b[1],rep(0,p)))
    if (length(act.set)>(n-intercept)) {
      # Take any n-intercept elements (which ones dont matter)
      act.set = act.set[seq(n-intercept)]
    }
    b.new = rep(0,p+1)
    if (intercept) b.new[c(1,1+act.set)] = lsfit(x[,act.set],y)$coef
    else b.new[1+act.set] = lsfit(x[,act.set],y,int=FALSE)$coef
    return(b.new)
  })
}

#' Predict function for lasso object.
predict.lasso = function(object, newx, s=NULL, label = FALSE) {
  if (missing(newx)) newx = object$x
  if (object$intercept) {
    if (!all(newx[,1] == 1)) {
      newx = cbind(rep(1,nrow(newx)),newx)
    }
  }
  # https://stackoverflow.com/questions/39184991/why-is-predict-glmnet-not-predicting-probabilities
  xbeta = newx %*% coef.lasso(object,s) # for the coef, it is a matrix, each column corresponds to a lambda value
  out = exp(xbeta)/(1+exp(xbeta)) # prob
  if (label) out = ifelse(out <= 0.5, 0 , 1)
  return(out)
}


# obtain the coef_glmnet as the refernce
# note: this coef_glmnet at one lambda is duplicated and augmented so that all the points have a coef
# this is for the admm proposal
coef_aug_one_lambda = function(coeff, x) {
  if (nrow(matrix(coeff)) == ncol(x)){
    coef_vec = as.vector(coeff)
    coef = t(replicate(nrow(x), coef_vec)) # duplicate into a matrix
  }
  if (nrow(matrix(coeff)) == (1+ncol(x))) { # because of the intercept, in admm, no intercept currently to compare, 04102018
    coef_vec = matrix(coeff)[-1,] # coef_glmnet is now a vector
    coef = t(replicate(nrow(x), coef_vec)) # duplicate into a matrix
  }

  out = coef
  return (out)
}

# this coef_aug function deal with a matrix of coefficient and each column correspond to a lambda
coef_aug = function(coeff, x) {

  if (class(coeff) != "matrix") coeff = as.matrix (coeff)
  if (min(dim(coeff)) == 1) {
    # indicate only one lambda
    out = coef_aug_one_lambda (coeff, x)
  } else {
    storage = vector(mode="list",length=ncol(coeff))
    for (j in 1:ncol(coeff)) {
      tmp = coef_aug_one_lambda (coeff[,j], x)
      storage[[j]] = as.vector(tmp)
    }
    out = do.call(cbind, storage)
  }
  return (out)
}

