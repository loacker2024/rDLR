get_tab = function(p, s, num_segs, ar1, type_beta_true, rho, what="misclass.cv", spatial="homo" ){
  file.list = system(paste0("ls rds/sim.n*.p",p,".s",s,".seg",num_segs,".",spatial,".ARone",ar1,".type_beta_true",type_beta_true, sprintf(".rho%0.2f", rho),".rds"),intern=TRUE) # add .ARone
  file.list

  if (what == "misclass.cv") {
    method.names = c("cvlasso", "Mars","VCM_fuse_grp", "VCM_fuse", "VCM_mylongfused", "basis_myspline", "basis_myquad")
  } else if (what == "pm" || what == "cicr" || what == "F1" || what == "fpos" || what == "fneg" || what == "nzs" || what == "tneg"  || what == "tpos" || what == "size"  ) {
    method.names = c("cvlasso", "VCM_fuse_grp", "VCM_fuse", "VCM_mylongfused", "basis_myspline", "basis_myquad")
  } else {
    method.names = c("cvlasso", "Mars","VCM_fuse_grp", "VCM_fuse", "VCM_mylongfused", "basis_myspline", "basis_myquad")
  }

  method.nums =  seq(1, length(method.names));
  N = length(method.nums)

  yvec = ybar = num_segs.vec = rho.vec = p.vec = n.vec = c()
  n_size = c(100,200,300,400,500)
  for (i in 1:length(file.list)) {
    sim.obj = readRDS(file.list[i]) # names(sim.obj)
    
    num_segs.vec = c(num_segs.vec,rep(5,N))
    rho.vec = c(rho.vec,rep(rho,N))
    p.vec = c(p.vec,rep(p,N))
    n.vec = c(n.vec,rep(n_size[i],N))

    z = sim.obj[[what]]

    res = z[method.names]
    res

    # For prop, F  and nonzero we ignore any request for a relative metric
    if (what=="prop" || what=="F1" || what=="nzs" || what=="fpos" || what=="tneg" ||
        what=="mer" || what=="err.val" ||
        what=="dev.cv" || what=="misclass.cv" || what == "err.train" || what == "size" ||
        what=="cicr" || what=="pm" || what=="tpos" || what =="fneg" || what == "runtime") {
      if (1) { #(what == "runtime") {
        # for all metrics, use this to gather results
        z.ave = z.std = z.med = z.mad = vector(mode="list",length=N)
        z = res
        for (j in 1:N) {
          # z[[j]] = matrix(z[[j]], nrow=nrep) # Just in case it is not a matrix
          z.ave[[j]] = apply(z[[j]], 2, mean, na.rm=TRUE) # colMeans(z[[j]], na.rm=TRUE)
          z.std[[j]] = apply(z[[j]], 2, sd, na.rm=TRUE) 
          z.med[[j]] = apply(z[[j]], 2, median, na.rm=TRUE)
          z.mad[[j]] = apply(z[[j]], 2, mad, na.rm=TRUE) # / sqrt(colSums(!is.na(z[[j]])))
        }
        yvec = c(yvec, z.ave[method.nums])
        ybar = c(ybar, z.std[method.nums])
      } else {
        yvec = c(yvec,res[[paste0("z.",substr(tuning,1,3),".",type)]][method.nums])
        ybar = c(ybar,res[[paste0("z.",substr(tuning,1,3),".",
                                  ifelse(type=="ave","std","mad"))]][method.nums])
      } # end if what == runtime condition
    }

  } # end for-i loop

  yvec = unlist(yvec)
  ybar = unlist(ybar)
  xvec = n.vec

  num_segs.vec = factor(num_segs.vec)
  rho.vec = factor(rho.vec)
  n.vec = factor(n.vec)
  levels(num_segs.vec) = paste("Number of Segments", levels(num_segs.vec))
  levels(rho.vec) = paste("Correlation", levels(rho.vec))

  dat = data.frame(x=xvec, y=yvec, se=ybar,
                   num_segs=num_segs.vec, rho=rho.vec, p=p.vec,
                   Method=factor(rep(method.names, length=length(xvec))),
                   typeBeta=type_beta_true, temporally=spatial, typeAR=ar1
  )
  colnames(dat)[colnames(dat)%in%"y"] = what
  colnames(dat)[colnames(dat)%in%"se"] = paste0(what,'.se')

  dat = base::transform(dat, Method=ifelse(Method=="cvlasso", "LASSO",
                                           ifelse(Method=="Mars", "MARS",
                                                  ifelse(Method=="VCM_fuse_grp", "rDLR",
                                                         ifelse(Method=="VCM_fuse", "FUSE",
                                                                ifelse(Method=="VCM_mygam", "method1",
                                                                       ifelse(Method=="VCM_mylongfused", "method2",
                                                                              ifelse(Method=="basis_myspline", "VCM1",
                                                                                     ifelse(Method=="basis_myquad", "VCM2","others")))))))))
  out = dplyr::filter(dat, Method %in% c("rDLR", "FUSE", "VCM1", "VCM2", "LASSO", "MARS"))
  out
}

