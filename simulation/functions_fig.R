library(ggplot2)
library(gridExtra)
library(grid)

# https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right","top","none")) {

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)),
                     "top" = arrangeGrob(do.call(arrangeGrob, gl),
                                         legend,
                                         ncol = 1,
                                         heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "none" = arrangeGrob(do.call(arrangeGrob, gl),
                                          legend,
                                          ncol = 1,
                                          heights = unit.c(unit(1, "npc") - lheight, lheight))
  )

  grid.newpage()
  grid.draw(combined)

  # return gtable invisibly
  invisible(combined)
}

plot.from.file = function(file.list,
                          row=c("num_segs","rho","n"), col=c("rho","num_segs","n"),
                          method.nums=NULL, method.names=NULL,
                          what=c("mer","pm","cicr","F1","fpos","fneg","runtime","nzs","tpos","err.train","dev.cv","misclass.cv"
                                 ,"tneg","err.val","size"), rel.to=NULL,
                          tuning=c("validation","oracle"), type=c("ave","med"),
                          plot_title = TRUE,
                          std=TRUE, lwd=1, point_size=3, pch=19, main=NULL, ylim=NULL,
                          legend.pos=c("bottom","right","top","left","none"),
                          make.pdf=FALSE, fig.dir=".", file.name="sim",
                          nobackground = TRUE,
                          n=NULL, p=NULL, s=NULL, num_segs=NULL, rho = NULL,ar1=NULL,
                          w=8, h=10) {

  row = match.arg(row)
  col = match.arg(col)
  if (row==col) stop("row and col must be different")

  what = match.arg(what)
  tuning = match.arg(tuning)
  type = match.arg(type)
  legend.pos = match.arg(legend.pos)

  if (what == "misclass.cv") {
    method.names = c("cvlasso", "Mars","VCM_fuse_grp", "VCM_fuse", "VCM_mylongfused", "basis_myspline", "basis_myquad")
  } else if (what == "pm" || what == "cicr" || what == "F1" || what == "fpos" || what == "fneg" || what == "nzs" || what == "tneg"  || what == "tpos" || what == "size"  ) {
    method.names = c("cvlasso", "VCM_fuse_grp", "VCM_fuse", "basis_myspline", "basis_myquad")
  } else {
    method.names = c("cvlasso", "Mars","VCM_fuse_grp", "VCM_fuse", "VCM_mylongfused", "basis_myspline", "basis_myquad")
  }

  method.nums =  seq(1, length(method.names));
  N = length(method.nums)

  # Set the base number and name
  if (is.null(rel.to)) {
    base.num = 0
    base.name = ifelse(what=="error","Bayes","null model")
  } else {
    base.num = which(method.nums==rel.to)
    base.name = tolower(method.names[base.num])
  }

  ylab = switch(what,
                error=paste0("Relative test error (to ",base.name,")"),
                risk=paste0("Relative risk (to ",base.name,")"),
                prop="Proportion of variance explained",
                F1="F1 classification of nonzeros",
                tpos = "True Positive",
                fpos = "False Positive",
                fneg = "False Negative",
                tneg = "True Negative",
                mer = "Misclassification error rate",
                pm = "Performance measurement",
                cicr = "Correctly identified coefficient rate",
                err.val = "Deviance (validation set)",
                dev.cv = "dev.cv", err.train = "err.train",
                misclass.cv = "misclass.cv",
                runtime = "Runtime/s",
                nzs="Number of nonzeros",
                size="Size"
  )
  ylab
  # Collect the y-variable from the file list
  yvec = ybar = num_segs.vec = rho.vec = p.vec = n.vec = c()
  for (i in 1:length(file.list)) {
    sim.obj = readRDS(file.list[i]) 

    n_tmp = c(100,200,300,400,500)[i]
    p_tmp = p
    num_segs_tmp = num_segs
    rho_tmp = rho
  
    n.vec = c(n.vec,rep(n_tmp,N))
    p.vec = c(p.vec,rep(p_tmp,N))
    num_segs.vec = c(num_segs.vec,rep(num_segs_tmp,N))
    rho.vec = c(rho.vec,rep(rho_tmp,N))

    z = sim.obj[[switch(what,
                        error="err.test", err.train="err.train", dev.cv="dev.cv", misclass.cv="misclass.cv",
                        risk="risk",
                        prop="prop",
                        F1="F1",
                        tpos = "tpos",
                        fpos = "fpos",
                        fneg = "fneg",
                        runtime = "runtime", size="size",
                        tneg = "tneg", mer="mer",pm="pm",cicr="cicr",err.val="err.val",
                        nzs="nzs")]]
      
    if (what == "runtime") {
      res = sim.obj$runtime 
    } else {
      res = z[method.names]
    }

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
    #
    else {
      # First build the relative metric
      met = res[[paste0("z.",substr(tuning,1,3))]][method.nums]
      if (base.num == 0 && what=="error") denom = sim.obj$sigma^2
      else if (base.num == 0 && what=="risk") denom = sim.obj$risk.null
      else denom = met[[base.num]]
      z.rel = lapply(met, function(v) v / denom)
      # Now aggregate the relative metric
      res2 = tune.and.aggregate(sim.obj, z.rel, tune=FALSE)
      yvec = c(yvec,unlist(res2[[paste0("z.",type)]]))
      ybar = c(ybar,unlist(res2[[paste0("z.",ifelse(type=="ave",
                                                    "std","mad"))]]))
    } # end if condition
  } # end for-i loop
  if (TRUE) { 
    yvec = unlist(yvec)
    ybar = unlist(ybar)
  }
  # Set the x-variable and x-label
  xvec = n.vec 
  xlab = "Number of observations"


  if (is.null(ylim)) ylim = range(yvec-ybar, yvec+ybar, na.rm=TRUE)
  num_segs.vec = factor(num_segs.vec)
  rho.vec = factor(rho.vec)
  n.vec = factor(n.vec)
  levels(num_segs.vec) = paste("Number of Segments", levels(num_segs.vec))
  levels(rho.vec) = paste("Correlation", levels(rho.vec))

  dat = data.frame(x=xvec, y=yvec, se=ybar,
                   num_segs=num_segs.vec, rho=rho.vec, p=p.vec,
                   Method=factor(rep(method.names, length=length(xvec))))

  # update the legend using the formal names
  dat = base::transform(dat, Method=ifelse(Method=="cvlasso", "LASSO",
                                           ifelse(Method=="Mars", "MARS",
                                                  ifelse(Method=="VCM_fuse_grp", "rDLR",
                                                         ifelse(Method=="VCM_fuse", "FUSE", 
                                                                ifelse(Method=="VCM_mygam", "VCM",
                                                                       ifelse(Method=="VCM_mylongfused", "method2",
                                                                              ifelse(Method=="basis_myspline", "VCM1",
                                                                                     ifelse(Method=="basis_myquad", "VCM2","others")))))))))

  gp = ggplot(dat, aes(x=x,y=y,color=Method)) +
    xlab(xlab) + ylab(ylab) + coord_cartesian(ylim=ylim) +
    geom_line(aes(linetype=Method), lwd=lwd) + geom_point(aes(pch=Method),size=point_size) +
    scale_shape_manual(values=c(1,2,3,4,5,8,11))+ #https://stackoverflow.com/questions/16813278/cycling-through-point-shapes-when-more-than-6-factor-levels
    #facet_grid(formula(paste(row,"~",col))) +
    theme_bw() + theme(legend.pos=legend.pos)
  if (nobackground) gp = gp +  theme(panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  if (!nobackground) gp = gp + facet_grid(formula(paste(row,"~",col)))
  # https://felixfan.github.io/ggplot2-remove-grid-background-margin/
  if (!("snr" %in% c(row,col))) {
    # If SNR is being plotted on the x-axis in each plot, then define special
    # x-axis ticks and put the x-axis on a log scale
    snr.breaks = round(exp(seq(from=min(log(xvec)),
                               to=max(log(xvec)),length=4)),2)
    # log-x-scale
    # gp = gp + scale_x_continuous(trans="log", breaks=snr.breaks)
  }
  if (plot_title) gp = gp + ggtitle(paste0("p=",unique(dat$p))) + theme(plot.title = element_text(hjust = 0.5))
  if (std) gp = gp + geom_errorbar(aes(ymin=y-se,ymax=y+se), width=0.02)

  if (!is.null(main)) gp = gp + ggtitle(main)
  if (!is.null(ylim)) gp = gp + coord_cartesian(ylim=ylim)

  if (make.pdf==TRUE) {
    ggsave(sprintf("%s/%s.pdf",fig.dir,file.name),
           height=h, width=w, device="pdf")
  } else if (make.pdf=="table") {
    return(dat)
  } else {
    gp
  }
}



pdf_fig_cases = function(n=NULL, p, s, num_segs = 5, rho = 0.35, case = "p20", spatial="homo", ar1=0
                         ,type_beta_true=1
                         , make.pdf = TRUE, fig.dir = "fig", file.name = paste(case,"",sep='')
                         , method.names = NULL, method.nums =  NULL,
                         plot_title = TRUE, legend.pos="right", type_plot = 1
                         ,nobackground = TRUE
){
  
  file.list = system(paste0("ls rds/sim.n*.p",p,".s",s,".seg",num_segs,".",spatial,".ARone",ar1,".type_beta_true",type_beta_true,".rho",rho,".rds"),intern=TRUE)
  file.list
  
  short.list = grep(gsub("\\.","\\\\.",paste0("*seg",num_segs,".",spatial,".ARone",ar1,".type_beta_true",type_beta_true,".rho",rho)),
                    file.list, val=TRUE)
  print(short.list)
  
  method.nums = method.names = NULL

  plot0 = plot.from.file(short.list, what="dev.cv", rel.to=NULL, tuning="val",
                         method.nums=method.nums, method.names=method.names,
                         legend.pos=legend.pos, make.pdf=make.pdf, fig.dir="fig",
                         plot_title = plot_title,
                         file.name=paste(case,".dev",sep=''), h=4, w=5.5
                         , n=n, p=p, s=s, num_segs=num_segs, rho = rho,ar1=ar1
                         ,nobackground = nobackground
  )

  plot1 = plot.from.file(short.list, what="misclass.cv", rel.to=NULL, tuning="val",
                         method.nums=method.nums, method.names=method.names,
                         legend.pos="right", make.pdf=make.pdf, fig.dir="fig",
                         plot_title = plot_title,
                         file.name=paste(case,".misclass",sep=''), h=4, w=5.5
                         , n=n, p=p, s=s, num_segs=num_segs, rho = rho,ar1=ar1
                         ,nobackground = nobackground
  )

  plot2 = plot.from.file(short.list, what="pm", rel.to=NULL, tuning="val",
                         method.nums=method.nums, method.names=method.names,
                         legend.pos=legend.pos, make.pdf=make.pdf, fig.dir="fig",
                         plot_title = plot_title,
                         file.name=paste(case,".pm",sep=''), h=4, w=4
                         , n=n, p=p, s=s, num_segs=num_segs, rho = rho,ar1=ar1
                         ,nobackground = nobackground
  )

  plot3 = plot.from.file(short.list, what="cicr", rel.to=NULL, tuning="val",
                         method.nums=method.nums, method.names=method.names,
                         legend.pos=legend.pos, make.pdf=make.pdf, fig.dir="fig",
                         plot_title = plot_title,
                         file.name=paste(case,".cicr",sep=''), h=4, w=4
                         , n=n, p=p, s=s, num_segs=num_segs, rho = rho,ar1=ar1
                         ,nobackground = nobackground
  )

  plot4 = plot.from.file(short.list, what="F1", rel.to=NULL, tuning="val",
                         method.nums=method.nums, method.names=method.names,
                         legend.pos=legend.pos, make.pdf=make.pdf, fig.dir="fig",
                         plot_title = plot_title,
                         file.name=paste(case,".F1",sep=''), h=4, w=4
                         , n=n, p=p, s=s, num_segs=num_segs, rho = rho,ar1=ar1
                         ,nobackground = nobackground
  )

  plot5 = plot.from.file(short.list, what="nzs", rel.to=NULL, tuning="val",
                         method.nums=method.nums, method.names=method.names,
                         legend.pos=legend.pos, make.pdf=make.pdf, fig.dir="fig",
                         plot_title = plot_title,
                         file.name=paste(case,".nonzero",sep=''), h=4, w=5.5
                         , n=n, p=p, s=s, num_segs=num_segs, rho = rho,ar1=ar1
                         ,nobackground = nobackground
  )

  plot6 = plot.from.file(short.list, what="fneg", rel.to=NULL, tuning="val",
                         method.nums=method.nums, method.names=method.names,
                         legend.pos=legend.pos, make.pdf=make.pdf, fig.dir="fig",
                         plot_title = plot_title,
                         file.name=paste(case,".fneg",sep=''), h=4, w=5.5
                         , n=n, p=p, s=s, num_segs=num_segs, rho = rho,ar1=ar1
                         ,nobackground = nobackground)

  plot7 = plot.from.file(short.list, what="fpos", rel.to=NULL, tuning="val",
                         method.nums=method.nums, method.names=method.names,
                         legend.pos=legend.pos, make.pdf=make.pdf, fig.dir="fig",
                         plot_title = plot_title,
                         file.name=paste(case,".fpos",sep=''), h=4, w=5.5
                         , n=n, p=p, s=s, num_segs=num_segs, rho = rho,ar1=ar1
                         ,nobackground = nobackground)

  plot8 = plot.from.file(short.list, what="tpos", rel.to=NULL, tuning="val",
                         method.nums=method.nums, method.names=method.names,
                         legend.pos=legend.pos, make.pdf=make.pdf, fig.dir="fig",
                         plot_title = plot_title,
                         file.name=paste(case,".tpos",sep=''), h=4, w=5.5
                         , n=n, p=p, s=s, num_segs=num_segs, rho = rho,ar1=ar1
                         ,nobackground = nobackground)

  plot9 = plot.from.file(short.list, what="tneg", rel.to=NULL, tuning="val",
                         method.nums=method.nums, method.names=method.names,
                         legend.pos=legend.pos, make.pdf=make.pdf, fig.dir="fig",
                         plot_title = plot_title,
                         file.name=paste(case,".tneg",sep=''), h=4, w=5.5
                         , n=n, p=p, s=s, num_segs=num_segs, rho = rho,ar1=ar1
                         ,nobackground = nobackground)

  plot10 = plot.from.file(short.list, what="err.train", rel.to=NULL, tuning="val",
                          method.nums=method.nums, method.names=method.names,
                          legend.pos=legend.pos, make.pdf=make.pdf, fig.dir="fig",
                          plot_title = plot_title,
                          file.name=paste(case,".err.train",sep=''), h=4, w=5.5
                          , n=n, p=p, s=s, num_segs=num_segs, rho = rho,ar1=ar1
                          ,nobackground = nobackground)

  if (make.pdf == FALSE) {
    # make combined pdf plots
    #plots1 = grid_arrange_shared_legend(plot0, plot10, ncol = 2, nrow = 1, position = "right")
    #plots1 = grid.arrange(plots1, plot1, nrow = 2, ncol=1, widths=c(2)) # b/c cv.dev and cv.class have different legends
    # https://stackoverflow.com/questions/40877386/grid-arrange-ggplot2-plots-by-columns-instead-of-by-row-using-lists
    #plots2 =  grid_arrange_shared_legend(plot6, plot7, plot8, plot9, ncol=2, nrow=2, position = "right")
    plots3 =  grid_arrange_shared_legend(plot2, plot3, plot4, plot5, ncol=2, nrow=2, position = "right")
    out = plots3 
  }
  if (type_plot == 1) { # another version of making plots to use in the paper
    plots3 = grid_arrange_shared_legend(plot2, plot3, plot4, ncol=3, nrow=1, position = "right")
    # plots1 = grid_arrange_shared_legend(plot0, plot10, ncol = 2, nrow = 1, position = "right")
    plots2 = grid.arrange(plot0, plot1, nrow = 1, ncol=2) # b/c cv.dev and cv.class have different legends
    plots1 = grid.arrange(plots2, plots3, nrow = 2, ncol=1) # b/c cv.dev and cv.class have different legends

    ggsave(sprintf("%s/%s_s%.0f_num_segs%.0f_%s_ARone%.1f_rho%.2f_large.pdf",fig.dir,file.name,s,num_segs,spatial,ar1,rho), plot = plots1,
           height=6, width=9, device="pdf")

    out = plots1

  } else if (type_plot == 2) {
    plots3 = grid_arrange_shared_legend(plot2, plot5, plot3, plot4, ncol=4, nrow=1, position = legend.pos)
    out = plots3
  } else if (type_plot == 3) { # no shared legend
    plots3 = grid.arrange(plot2, plot5, plot3, plot4, ncol=4, nrow=1)
    out = plots3
  } else if (type_plot == 4) { # for the format used in ppt/poster
    plots4 = grid_arrange_shared_legend(plot2, plot4, plot3, plot5, ncol=4, nrow=1, position = "bottom")
    out = plots4
  }else if (type_plot == 5) { # for the format used in ppt/poster, no shared legend
    plots5 = grid.arrange(plot2, plot4, plot3, plot5, ncol=4, nrow=1)
    out = plots5
  }
  return(out)

}


foo = function(make.pdf=TRUE, type = 1, n=NULL, p, s, num_segs = 5, rho = 0.35, case = "p20", spatial="homo",ar1=0
               , type_beta_true=1
               , fig.dir = "fig", file.name = paste(case,"",sep='')
               , method.names = NULL, method.nums =  NULL,legend.pos="right", type_plot = 1, nobackground = TRUE
) {

  if (type == 1) {
    pdf_fig_cases (n=NULL, p=p, s=s, num_segs=num_segs, rho=rho, case=case, spatial=spatial,ar1=ar1, type_beta_true=type_beta_true
                   , make.pdf=make.pdf, fig.dir=fig.dir,legend.pos=legend.pos, type_plot = type_plot, nobackground=nobackground
    )
  } 
}


