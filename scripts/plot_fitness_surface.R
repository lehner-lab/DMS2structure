#################### plot A-B-AB surfaces #################### 
#subfunction for call_epistasis class of functions
plot_fitness_surface = function(double_data,F_fit_loess_model,List,dataset_dir,prefix) {
  
  #range of data to plot, omitting long tails in single mutant space
  xyrange = double_data[is.fitness == T & is.reads0==T,quantile(fitness1,c(0.005,0.995),na.rm=T)]
  xy = seq(xyrange[1],
           xyrange[2],
           abs(diff(xyrange))/25)
  
  Fd_pred = predict(List$F_median_fit,data.frame(fitness1=rep(xy,length(xy)),fitness2=rep(xy,each=length(xy))))
  Fd_pred2 =  predict(F_fit_loess_model,data.frame(fitness1=rep(xy,length(xy)),fitness2=rep(xy,each=length(xy))))
  
  Fd_pred_q05 = predict(List$F_lower_fit,data.frame(fitness1=rep(xy,length(xy)),fitness2=rep(xy,each=length(xy))))
  Fd_pred_q95 = predict(List$F_upper_fit,data.frame(fitness1=rep(xy,length(xy)),fitness2=rep(xy,each=length(xy))))
  
  xyz = double_data[is.fitness == T & is.reads0==T,between(fitness1,xyrange[1],xyrange[2]) & between(fitness2,xyrange[1],xyrange[2]),
                    .(x=fitness1,y=fitness2,z=fitness,below_q05=fitness < F_fit_lower,above_q95 = fitness>F_fit_upper)]
  
  #number of points to plot, 10k is sufficient, otherwise PDF becomes very large
  r = sample(x = nrow(xyz),size = min(c(10000,nrow(xyz))))
  x=xyz$x
  y=xyz$y
  z=xyz$z
  
  #axis limits, adjust
  xlim_plot = ylim_plot = c(xyrange[1] - 0.1*diff(xyrange),xyrange[2] + 0.1*diff(xyrange))
  zlim_plot = quantile(xyz$z,probs = c(0.005,0.995),na.rm = T)
  #plot angles
  theta_plot = c(15,55)
  phi_plot = 15 
  for (idx in 1:length(theta_plot)) {
    #### upper and lower surface with points
    pdf(paste0(dataset_dir, "results/epistasis/",prefix,"epistasis_surface",idx,".pdf"), useDingbats=FALSE)
    a=persp(xy,xy,matrix(Fd_pred2+Fd_pred_q05,nrow=length(xy),ncol=length(xy)),
            xlab="single mutant fitness 1",ylab="single mutant fitness 2",zlab="double mutant fitness",
            xlim = xlim_plot,ylim = ylim_plot,zlim = zlim_plot, 
            theta = theta_plot[idx], phi = phi_plot, 
            col=NA, nticks=5,ticktype="detailed",expand=0.75, box=TRUE)
    b=trans3d(xyz[intersect(r,which(below_q05==T))]$x,
              xyz[intersect(r,which(below_q05==T))]$y,
              xyz[intersect(r,which(below_q05==T))]$z,a)
    points(b$x,b$y,col=rgb(1,0.1,0.1),pch=16,cex=0.75)
    par(new=TRUE)
    a=persp(xy,xy,matrix(Fd_pred2+Fd_pred_q05,nrow=length(xy),ncol=length(xy)),
            xlab="single mutant fitness 1",ylab="single mutant fitness 2",zlab="double mutant fitness",
            xlim = xlim_plot,ylim = ylim_plot,zlim = zlim_plot,
            theta = theta_plot[idx], phi = phi_plot,
            col=NA, nticks=5,ticktype="detailed",expand=0.75, box=TRUE)
    par(new=TRUE)
    b=trans3d(xyz[intersect(r,which(below_q05==F & above_q95==F))]$x,
              xyz[intersect(r,which(below_q05==F & above_q95==F))]$y,
              xyz[intersect(r,which(below_q05==F & above_q95==F))]$z,a)
    points(b$x,b$y,col=rgb(1,0.5,0.5),pch=16,cex=0.75)
    par(new=TRUE)
    a2=persp(xy,xy,matrix(Fd_pred2+Fd_pred_q95,nrow=length(xy),ncol=length(xy)),
             xlab="single mutant fitness 1",ylab="single mutant fitness 2",zlab="double mutant fitness",
             xlim = xlim_plot,ylim = ylim_plot,zlim = zlim_plot,
             theta = theta_plot[idx], phi = phi_plot,
             col=NA, nticks=5,ticktype="detailed",expand=0.75, box=TRUE)
    par(new=TRUE)
    b=trans3d(xyz[intersect(r,which(above_q95==T))]$x,
              xyz[intersect(r,which(above_q95==T))]$y,
              xyz[intersect(r,which(above_q95==T))]$z,a)
    points(b$x,b$y,col=rgb(0.5,0.9,0.1),pch=16,cex=0.75)
    dev.off()
  }
}