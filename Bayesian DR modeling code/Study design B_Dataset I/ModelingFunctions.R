##create/organize data for STAN to utilize in modeling
make_stan_dat <- function(dat,chemnum, parmcol, parmname="", parmnamenorm="",
                          quants = c(0.01,0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975,0.99)) {
  dat$cell.f<-factor(dat$Cell.line)
  dat$y <- dat[,parmcol];
  chemname<-ChemKey[as.character(chemnum),1]
  dat.tmp<-subset(dat,(Chemical.Number == chemnum | is.na(Chemical.Number)))
  dat.tmp<-subset(dat.tmp,!is.na(y))
  ###
  ### For all other chemicals, need to do additional subset for same plate-controls
  ###
  fileprefix<-paste("C", chemnum,"_","standat_v7_do_stan_fit_",parmnamenorm,sep="")    ##No ACO
  cat("##############",fileprefix,"\n")
  stan_dat <- list(
    scale_factor = median(dat$y,na.rm=TRUE), # dat not dat.tmp
    Ni = length(unique(dat$cell.f)) + 1, # dat not dat.tmp
    Nj = dim(dat.tmp)[1],
    x = dat.tmp$Chemical.Concentration,
    ys = dat.tmp$y/median(dat$y,na.rm=TRUE),
    cell = as.numeric(dat.tmp$cell.f),
    quants = quants,
    Nquants = length(quants)
  )
  
  with(stan_dat, {
    stan_rdump(names(stan_dat),file=paste(chempath, "/", fileprefix,"_dat.R",sep=""))
  })
  list(stan_dat=stan_dat,fileprefix=fileprefix,chemname=chemname,parmcol=parmcol,
       parmname=parmname,parmnamenorm=parmnamenorm)
}
#specify the model to use based on direction of data and specify iteration

# iteration number specified in chemmap

# increased maximum treedepth from default 10 to 15

# increased adapt_alpha from default 0.8 to 0.99

# seed changed to 669661

do_stan_fit <- function(stan_dat, fileprefix="",direction = 1,
                        verbose=FALSE,seed=669661,iter=4000,
                        ...) {
  if (direction == 1) {
    #need to change paths accordingly
    stanmodel <- 'conc_resp_v7_up.stan' 
  } else if (direction == -1) {
    stanmodel <- 'conc_resp_v7_dn.stan' 
  } else if (direction == 0){
    stanmodel <- 'conc_resp_v8_zero.stan' 
  }
  
  time.start<-proc.time()
  stan_fit <- stan(file=stanmodel,data=stan_dat,verbose=verbose,seed=seed,iter=iter, save_warmup=FALSE, chains = 4,
                   control = list(max_treedepth = 15, adapt_delta = 0.99),
                   sample_file=paste(chempath,"/",fileprefix,"_samples",sep=""),...);
  time.end<-proc.time()
  print(time.end-time.start)
  stan_samples<-extract(stan_fit)
  list(stan_dat=stan_dat,stan_fit=stan_fit,stan_samples=stan_samples,direction=direction)
}

get_stan_fit_csv <- function(stan_dat,fileprefix,direction) {
  stan_fit <- read_stan_csv(dir(path = chempath, pattern=paste(fileprefix,"_samples_","[[:digit:]]",sep=""), 
                                full.names = TRUE))
  #  stan_fit[[1]] <- NULL # get rid of "energy"
  stan_samples<-extract(stan_fit) 
  #  stan_samples[[1]]<-NULL # get rid of "energy"
  list(stan_dat=stan_dat,stan_fit=stan_fit,stan_samples=stan_samples,direction=direction)
}

plot_stan_fit_rhat <- function(stan_fit,stan_samples,fileprefix,dopdf=FALSE) {    
  if (dopdf) pdf(paste(fileprefix,"rhat","pdf",sep="."),height=10.5,width=8)
  print(plot(stan_fit,plotfun="stan_rhat",par=names(stan_samples)))
  if (dopdf) dev.off()
}

plot_stan_fit_parms <- function(stan_fit,fileprefix, dopdf=FALSE) { 
  if (dopdf) pdf(paste(fileprefix,"parms","pdf",sep="."),height=10.5,width=8)
  .print(stan_plot(stan_fit,
                  par=c("m_y0","m_x0","m_Emax","m_n","sd_y0","sd_x0","sd_Emax","sd_n","sigma_y"))+ 
          ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="z_y0")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="z_x0")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="z_Emax")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="z_n")+ ggtitle(fileprefix))
  if (dopdf) dev.off()
}

plot_stan_fit_pop <- function(stan_fit,fileprefix, dopdf=FALSE) { ##same
  if (dopdf) pdf(paste(fileprefix,"pop","pdf",sep="."),height=10.5,width=8)
  print(stan_plot(stan_fit,par="y0")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="x0")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="Emax")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="n")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="ec99")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="ec95")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="ec50")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="ec10")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="ec05")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="ec025")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="ec01")+ ggtitle(fileprefix))
  
  if (dopdf) dev.off()
}

plot_stan_fit_iter <- function(stan_samples,stan_dat,direction,chemname,parmnamenorm, 
                               fileprefix="",dopdf=FALSE,niter=50) {
  if (dopdf) pdf(paste(fileprefix,"iter","pdf",sep="."),height=10.5,width=8)
  par(mfrow=c(9,3),mar=c(0,2,0.5,1),oma=c(4,3,3,1))
  x50 <- stan_samples$x0^stan_samples$n*stan_samples$Emax
  xmin=max(1e-6,min(0.001,10^floor(log10(quantile(x50,0.1,na.rm=TRUE))),na.rm=TRUE))
  yrange=range(stan_dat$scale_factor*stan_dat$ys)
  for (i in 1:stan_dat$Ni) {
    xtmp <- stan_dat$x[stan_dat$cell == i]+xmin
    ytmp <- stan_dat$scale_factor*stan_dat$ys[stan_dat$cell == i]
    yrange=range(ytmp,na.rm=TRUE)
    xx<- 10^(((10*(ceiling(log10(xmin)))):30)/10)
    lastiter<-dim(stan_samples$y0)[1]
    for (iter in 0:(niter-1)) {
      y0<-stan_samples$y0[lastiter-iter,i]
      x0<-stan_samples$x0[lastiter-iter,i]
      Emax<-stan_samples$Emax[lastiter-iter,i]
      n<-stan_samples$n[lastiter-iter,i]
      if (iter == 0) {
        plot(xx,y0*stan_dat$scale_factor*(1 + direction * (xx / x0)^n / (1 + (xx / x0)^n/Emax)),
             ylim=yrange, 
             #c(yrange[1]-0.25*(yrange[2]-yrange[1]),yrange[2]),
             type="l",log="x",col=i,xlab="",ylab="",axes=FALSE,lwd=0.5);
        axis(2,line=0)
        axis(1,at=10^(log10(xmin):2),lab=c("C",log10(xmin*10):2),line=0,outer=TRUE)
        box()
      } else {
        lines(xx,y0*stan_dat$scale_factor*(1 + direction * (xx / x0)^n / (1 + (xx / x0)^n/Emax)),col=i,lwd=0.5)
      }
    }
    xec99<-10^seq(max(log10(xmin),quantile(log10(stan_samples$ec99[,i]),prob=c(0.05))),
                  quantile(log10(stan_samples$ec99[,i]),prob=c(0.95)),
                  length.out=100)
    yec99<-stan_dat$scale_factor*rep(median(stan_samples$y0[,i])*(1+0.1*direction),100)
    points(xec99,yec99,pch=15)
    points(xec99,yec99,pch=15,col="grey",cex=0.5)
    points(xtmp,ytmp,pch=21,bg="white")
    text(xmin,yrange[2],lab=i,adj=c(0,1))
    
    xec95<-10^seq(max(log10(xmin),quantile(log10(stan_samples$ec95[,i]),prob=c(0.05))),
                  quantile(log10(stan_samples$ec95[,i]),prob=c(0.95)),
                  length.out=100)
    yec95<-stan_dat$scale_factor*rep(median(stan_samples$y0[,i])*(1+0.1*direction),100)
    points(xec95,yec95,pch=15)
    points(xec95,yec95,pch=15,col="grey",cex=0.5)
    points(xtmp,ytmp,pch=21,bg="white")
    text(xmin,yrange[2],lab=i,adj=c(0,1))
    
    xec50<-10^seq(max(log10(xmin),quantile(log10(stan_samples$ec50[,i]),prob=c(0.05))),
                  quantile(log10(stan_samples$ec50[,i]),prob=c(0.95)),
                  length.out=100)
    yec50<-stan_dat$scale_factor*rep(median(stan_samples$y0[,i])*(1+0.1*direction),100)
    points(xec50,yec50,pch=15)
    points(xec50,yec50,pch=15,col="grey",cex=0.5)
    points(xtmp,ytmp,pch=21,bg="white")
    text(xmin,yrange[2],lab=i,adj=c(0,1))
    
    xec10<-10^seq(max(log10(xmin),quantile(log10(stan_samples$ec10[,i]),prob=c(0.05))),
                  quantile(log10(stan_samples$ec10[,i]),prob=c(0.95)),
                  length.out=100)
    yec10<-stan_dat$scale_factor*rep(median(stan_samples$y0[,i])*(1+0.1*direction),100)
    points(xec10,yec10,pch=15)
    points(xec10,yec10,pch=15,col="grey",cex=0.5)
    points(xtmp,ytmp,pch=21,bg="white")
    text(xmin,yrange[2],lab=i,adj=c(0,1))
    
    xec05<-10^seq(max(log10(xmin),quantile(log10(stan_samples$ec05[,i]),prob=c(0.05))),
                  quantile(log10(stan_samples$ec05[,i]),prob=c(0.95)),
                  length.out=100)
    yec05<-stan_dat$scale_factor*rep(median(stan_samples$y0[,i])*(1+0.1*direction),100)
    points(xec05,yec05,pch=15)
    points(xec05,yec05,pch=15,col="grey",cex=0.5)
    points(xtmp,ytmp,pch=21,bg="white")
    text(xmin,yrange[2],lab=i,adj=c(0,1))
    
    xec025<-10^seq(max(log10(xmin),quantile(log10(stan_samples$ec025[,i]),prob=c(0.05))),
                   quantile(log10(stan_samples$ec025[,i]),prob=c(0.95)),
                   length.out=100)
    yec025<-stan_dat$scale_factor*rep(median(stan_samples$y0[,i])*(1+0.1*direction),100)
    points(xec025,yec025,pch=15)
    points(xec025,yec025,pch=15,col="grey",cex=0.5)
    points(xtmp,ytmp,pch=21,bg="white")
    text(xmin,yrange[2],lab=i,adj=c(0,1))
    
    xec01<-10^seq(max(log10(xmin),quantile(log10(stan_samples$ec01[,i]),prob=c(0.05))),
                  quantile(log10(stan_samples$ec01[,i]),prob=c(0.95)),
                  length.out=100)
    yec01<-stan_dat$scale_factor*rep(median(stan_samples$y0[,i])*(1+0.1*direction),100)
    points(xec01,yec01,pch=15)
    points(xec01,yec01,pch=15,col="grey",cex=0.5)
    points(xtmp,ytmp,pch=21,bg="white")
    text(xmin,yrange[2],lab=i,adj=c(0,1))
  }
  
  mtext(paste(chemname,parmnamenorm),line=1,outer=TRUE)
  mtext("Log10 Concentration",side=1,line=2.2,outer=TRUE)
  mtext(parmnamenorm,side=2,line=1,outer=TRUE)
  if (dopdf) dev.off()
}