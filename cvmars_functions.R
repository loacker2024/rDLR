library(earth)
#'  # fit ab addutuve MARS data
#' “Multivariate Adaptive Regression Splines”

# MARS works by splitting input variables into multiple basis functions and then fitting a linear regression model to those basis functions. The basis
# functions used by MARS come in pairs: f(x) = {x − t if x > t, 0 otherwise} and g(x)
# = {t − x if x < t, 0 otherwise}. These functions are piecewise linear functions. The
# value t is called a knot.

sim_cvMARS = function(x, y, numFolds = 5, plot = FALSE, trace = 0){
 # intercept is always in the formula
  if (class(x) != "matrix") x = as.matrix(x)

  fit <- earth(x = x, y = y
               #,dat
               ,trace = trace       #Trace = 1: shows the details of computation
               ,ncross= 1
               ,nfold = numFolds
               ,pmethod="backward"
               ,nprune= nrow(x)*0.5
               , glm = list(family = binomial, control = list(maxit = 1000))
  )

  if (plot == TRUE) plotmo(fit)

  out = fit

  out$x = x; out$y = y
  out$intercept = TRUE # verify by fit$prune.terms
  out$numFolds = numFolds

  # https://cran.r-project.org/web/packages/earth/earth.pdf
  out$cvm_dev = fit$cv.deviance.tab[nrow(fit$cv.deviance.tab), ncol(fit$cv.deviance.tab)]
  out$cvm_misclass = 1- fit$cv.class.rate.tab[nrow(fit$cv.class.rate.tab), ncol(fit$cv.class.rate.tab)] #the fraction ofclasses correctly predicted

  class(out) = c("myearth", "earth")
  # https://cran.r-project.org/web/packages/earth/earth.pdf

  # ref
  # https://cran.r-project.org/web/packages/earth/earth.pdf
  # https://github.com/fclesio/learning-space/blob/master/R/mars-regression-study-draft.R
  # http://www.milbo.org/doc/earth-notes.pdf
  return(out)
}

coef.myearch = function(object){
 # 
}


predict.myearth = function(object, newx, label = FALSE){

  if (missing(newx)) newx = object$x

  class(object) = "earth"
  out = predict(object, newdata = newx,type = "response")
  if (label) out = ifelse(out <= 0.5, 0 , 1)
  return(out)
}


