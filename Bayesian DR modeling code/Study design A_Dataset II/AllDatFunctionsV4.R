library(reshape2)

## DR Functions ##

DR_up<-function(x=1, y0=1, x0=1, Emax=1, n=1){
  y0*(1+((x/x0)^n)/(1+((x/x0)^n)*(1/Emax))) 
}

DR_dn<-function(x=1, y0=1, x0=1, Emax=1, n=1){
  y0*(1-((x/x0)^n)/(1+((x/x0)^n)/(Emax))) 
}

DR_zero<-function(x=1, y0=1, x0=1, n=1){
  y0*(1-((x/x0)^n)/(1+((x/x0)^n)/1)) 
}

DR_up_Frac<- function(x=1, x0=1, Emax=1, n=1){
  (x/x0)^n/(1+((x/x0)^n)*(1/Emax))
}

## function to remove outliers from data based on quartiles

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

## Function to read in data by chemical name and number ##

ChemDatRead<- function(chemical.name = chemical.name, chemical.num = chemical.num, endpoint = endpoint){
  path.paste<- paste(getwd(), "/", chemical.name, sep = "")
  pattern.R<- paste("C",chemical.num, "_standat_v7_do_stan_fit_", endpoint, "_dat.R",sep="")
  pattern.csv<- paste("C", chemical.num,"_standat_v7_do_stan_fit_", endpoint, "_samples_*",sep="")
  file.sources=list.files(path=path.paste,pattern=pattern.R)
  sapply(paste(path.paste, "/", file.sources,sep=""),source,.GlobalEnv)
  files = list.files(path=path.paste , pattern= pattern.csv)
  return(files)
}

SamplesDatRead <- function(chemical.name = chemical.name, chemical.num = chemical.num, endpoint = endpoint){
  path.paste<- paste(getwd(), "/", chemical.name, sep = "")
  pattern.R<- paste("C",chemical.num, "_standat_v7_do_stan_fit_", endpoint, "_dat.R",sep="")
  pattern.csv<- paste(path.paste,"/","C", chemical.num,"_standat_v7_do_stan_fit_", endpoint, "_MCMC_Samples.csv",sep="")
  file.sources=list.files(path=path.paste,pattern=pattern.R)
  sapply(paste(path.paste, "/", file.sources,sep=""),source,.GlobalEnv)
  MCMC_Samples <- read.csv(pattern.csv)
  return(MCMC_Samples)
}

## Function to associate endpoint with proper direction ##
ep.compass<- function(endpoint = endpoint, return.dir = FALSE){
  endpoint.directions<-data.frame(endpoint.name = c("Peak_Width_up","Peak_Width_dn", "Peak_Frequency_up", "BPM_up", "Peak_Frequency_dn", "BPM_dn", "Peak_Decay_Rise_Ratio_up", "Peak_Amplitude_zero", "Peak_Frequency_zero", "Total_Cells_dn", "Total_Cells_zero", "Peak_Spacing_CV_up"),
                                  endoint.direction = c(1,-1,1,1,-1,-1,1,0,0, -1, 0, 1))
  row.num<- match(endpoint, endpoint.directions$endpoint.name)
  dir<- endpoint.directions[row.num, 2]
  if (dir == 1){
    DR_func<- DR_up
  } else if (dir == -1) {
    DR_func<- DR_dn
  } else if (dir == 0){
    DR_func <- DR_zero
  }
  if(return.dir == FALSE){
  return(DR_func)
  } else if (return.dir == TRUE){
    return(dir)
  }
}

## Function to return TDVF catagory based on chems median TDVF in comparison to 33 percentila and 66 percentile ##
TDVFCat<- function(x = x){
  if (x < TDVF33){
    category<- "Narrow"
  } else if (x > TDVF33 & x < TDVF66){
    category<- "Mid"
  } else if ( x > TDVF66) {
    category<- "Wide"
  }
  return(category)
}

## Function to calc quantiles and apply TDVF cat function ##
TDVFQuantCat<- function(df, dfall){
  TDVF33<- quantile(dfall[,3], 0.33)
  TDVF66<- quantile(dfall[,3], 0.66)
  df.temp<- as.data.frame(df[,3])
  TDVF.Category<- apply(df.temp, MARGIN = 1, FUN = TDVFCat)
  TDVF.Category<- cbind(df, TDVF.Category)
  return(TDVF.Category)
  
}

##Function to extract parms

ParmsExtract<- function(files = files){
  ParmsList<-list()
  y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, 186:227], (eval(parse(text=files[2])))[2:2001, 186:227], (eval(parse(text=files[3])))[2:2001, 186:227], (eval(parse(text=files[4])))[2:2001, 186:227]))
  colnames(y0)<- (1:42)
  rownames(y0)<- c(1:8000)
  ParmsList$y0<- y0
  x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, 228:269], (eval(parse(text=files[2])))[2:2001, 228:269], (eval(parse(text=files[3])))[2:2001, 228:269], (eval(parse(text=files[4])))[2:2001, 228:269]))
  colnames(x0)<- paste("ID",1:42, sep="")
  rownames(x0)<- c(1:8000)
  ParmsList$x0<- x0
  Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, 270:311], (eval(parse(text=files[2])))[2:2001, 270:311], (eval(parse(text=files[3])))[2:2001, 270:311], (eval(parse(text=files[4])))[2:2001, 270:311]))
  colnames(Emax)<- paste("ID",1:42, sep="")
  rownames(Emax)<- c(1:8000)
  ParmsList$Emax<- Emax
  n<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, 312:353], (eval(parse(text=files[2])))[2:2001, 312:353], (eval(parse(text=files[3])))[2:2001, 312:353], (eval(parse(text=files[4])))[2:2001, 312:353]))
  colnames(n)<- paste("ID",1:42, sep="")
  rownames(n)<- c(1:8000)
  ParmsList$n<- n
  return(ParmsList)
}

ParmsExtractSamples <- function(samples = samples){
  dir <- ep.compass(endpoint = endpoint, return.dir = T)
  ParmsList <- list()
  y0<-samples[,which(colnames(samples) == "y0.1"):which(colnames(samples) == "y0.42")]
  colnames(y0) <- paste("ID",1:42, sep = "")
  rownames(y0)<-1:1000
  ParmsList$y0 <- y0
  x0<-samples[,which(colnames(samples) == "x0.1"):which(colnames(samples) == "x0.42")]
  colnames(x0) <- paste("ID",1:42, sep = "")
  rownames(x0)<-1:1000
  ParmsList$x0 <- x0
  if(dir != 0){
    Emax<-samples[,which(colnames(samples) == "Emax.1"):which(colnames(samples) == "Emax.42")]
    colnames(Emax) <- paste("ID",1:42, sep = "")
    rownames(Emax)<-1:1000
    ParmsList$Emax <- Emax
  }
  n<-samples[,which(colnames(samples) == "n.1"):which(colnames(samples) == "n.42")]
  colnames(n) <- paste("ID",1:42, sep = "")
  rownames(n)<-1:1000
  ParmsList$n <- n
  return(ParmsList)
}

## Function to Generate prediction frame ##
PredFrameGen<- function(files = files, returnparms = F ){
  concframe<-data.frame(x=c(1e-6,10^(seq(-3,2,0.05))))
  y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, 186:227], (eval(parse(text=files[2])))[2:2001, 186:227], (eval(parse(text=files[3])))[2:2001, 186:227], (eval(parse(text=files[4])))[2:2001, 186:227]))
  colnames(y0)<- (1:42)
  rownames(y0)<- c(1:8000)
  parmframe<- melt(as.matrix(y0))
  colnames(parmframe)<- c("Iter", "ID", "y0")
  x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, 228:269], (eval(parse(text=files[2])))[2:2001, 228:269], (eval(parse(text=files[3])))[2:2001, 228:269], (eval(parse(text=files[4])))[2:2001, 228:269]))
  colnames(x0)<- paste("ID",1:42, sep="")
  rownames(x0)<- c(1:8000)
  parmframe$x0<- melt(as.matrix(x0))$value
  Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, 270:311], (eval(parse(text=files[2])))[2:2001, 270:311], (eval(parse(text=files[3])))[2:2001, 270:311], (eval(parse(text=files[4])))[2:2001, 270:311]))
  colnames(Emax)<- paste("ID",1:42, sep="")
  rownames(Emax)<- c(1:8000)
  parmframe$Emax<- melt(as.matrix(Emax))$value
  n<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, 312:353], (eval(parse(text=files[2])))[2:2001, 312:353], (eval(parse(text=files[3])))[2:2001, 312:353], (eval(parse(text=files[4])))[2:2001, 312:353]))
  colnames(n)<- paste("ID",1:42, sep="")
  rownames(n)<- c(1:8000)
  parmframe$n<- melt(as.matrix(n))$value
  if(returnparms == T){
    return(parmframe)
  } else {
  ### Concentration-response
  Prediction_Frame<-data.frame()
  for (k in 1:nrow(concframe)){
    cat(k,"..")
    parmframeconc<- parmframe
    parmframeconc$concentration<-concframe$x[k]
    parmframeconc$Pred<- DR_func(x=parmframeconc$concentration, 
                                 y0=parmframeconc$y0, 
                                 x0=parmframeconc$x0, 
                                 Emax=parmframeconc$Emax, 
                                 n=parmframeconc$n)* scale_factor
    parmframeconc$concentration<-concframe$x[k]
    predframe_median<- aggregate(Pred~concentration+ID, data=parmframeconc, quantile, prob=0.5)
    predframe_975<- aggregate(Pred~concentration+ID, data=parmframeconc, quantile, prob=0.975)
    predframe_25<- aggregate(Pred~concentration+ID, data=parmframeconc, quantile, prob=0.025)
    Prediction_Frame<- rbind(Prediction_Frame, cbind(predframe_median, predframe_25$Pred,predframe_975$Pred))
  }
  names(Prediction_Frame)[2]<- "Individual"
  names(Prediction_Frame)[3]<- "p50"
  names(Prediction_Frame)[4]<- "p2.5"
  names(Prediction_Frame)[5]<- "p97.5"
  return(Prediction_Frame)
  }
}

PredFrameGenSamples <- function(samples = samples, returnparms = F){
  concframe<-data.frame(x=c(1e-6,10^(seq(-3,2,0.05))))
    y0<-samples[,which(colnames(samples) == "y0.1"):which(colnames(samples) == "y0.42")]
    colnames(y0)<- 1:42
    rownames(y0)<- c(1:1000)
    parmframe<- melt(as.matrix(y0))
    colnames(parmframe)<- c("Iter", "ID", "y0")
    x0<-samples[,which(colnames(samples) == "x0.1"):which(colnames(samples) == "x0.42")]
    colnames(x0)<- 1:42
    rownames(x0)<- c(1:1000)
    parmframe$x0<- melt(as.matrix(x0))$value  
    if(dir !=0){
    Emax<-samples[,which(colnames(samples) == "Emax.1"):which(colnames(samples) == "Emax.42")]
    colnames(Emax)<- 1:42
    rownames(Emax)<- c(1:1000)
    parmframe$Emax<- melt(as.matrix(Emax))$value
    }
    n<-samples[,which(colnames(samples) == "n.1"):which(colnames(samples) == "n.42")]
    colnames(n)<- 1:42
    rownames(n)<- c(1:1000)
    parmframe$n<- melt(as.matrix(n))$value
  if(returnparms == T){
    return(parmframe)
  } else {
    ### Concentration-response
    Prediction_Frame<-data.frame()
    for (k in 1:nrow(concframe)){
      cat(k,"..")
      parmframeconc<- parmframe
      parmframeconc$concentration<-concframe$x[k]
      if(dir != 0){
      parmframeconc$Pred<- DR_func(x=parmframeconc$concentration, 
                                   y0=parmframeconc$y0, 
                                   x0=parmframeconc$x0, 
                                   Emax=parmframeconc$Emax, 
                                   n=parmframeconc$n)* scale_factor
      } else if (dir == 0){
        parmframeconc$Pred<- DR_func(x=parmframeconc$concentration, 
                                     y0=parmframeconc$y0, 
                                     x0=parmframeconc$x0, 
                                     n=parmframeconc$n)* scale_factor  
      }
      parmframeconc$concentration<-concframe$x[k]
      predframe_median<- aggregate(Pred~concentration+ID, data=parmframeconc, quantile, prob=0.5)
      predframe_975<- aggregate(Pred~concentration+ID, data=parmframeconc, quantile, prob=0.975)
      predframe_25<- aggregate(Pred~concentration+ID, data=parmframeconc, quantile, prob=0.025)
      Prediction_Frame<- rbind(Prediction_Frame, cbind(predframe_median, predframe_25$Pred,predframe_975$Pred))
    }
    names(Prediction_Frame)[2]<- "Individual"
    names(Prediction_Frame)[3]<- "p50"
    names(Prediction_Frame)[4]<- "p2.5"
    names(Prediction_Frame)[5]<- "p97.5"
    return(Prediction_Frame)
  }
}

PredFrameGenSamples.for1434 <- function(samples = samples, returnparms = F){
  concframe<-data.frame(x=c(1e-6,10^(seq(-3,2,0.05))))
  y0<-samples[,which(colnames(samples) == "y0.1")]
  parmframe<- data.frame(Iter = 1:1000,
                         ID = 1,
                         y0 = samples[,which(colnames(samples) == "y0.1")])
  parmframe$x0<- samples[,which(colnames(samples) == "x0.1")]  
  if(dir !=0){
    parmframe$Emax<- samples[,which(colnames(samples) == "Emax.1")]
  }
  parmframe$n<- samples[,which(colnames(samples) == "n.1")]
  if(returnparms == T){
    return(parmframe)
  } else {
    ### Concentration-response
    Prediction_Frame<-data.frame()
    for (k in 1:nrow(concframe)){
      cat(k,"..")
      parmframeconc<- parmframe
      parmframeconc$concentration<-concframe$x[k]
      if(dir != 0){
        parmframeconc$Pred<- DR_func(x=parmframeconc$concentration, 
                                     y0=parmframeconc$y0, 
                                     x0=parmframeconc$x0, 
                                     Emax=parmframeconc$Emax, 
                                     n=parmframeconc$n)* scale_factor
      } else if (dir == 0){
        parmframeconc$Pred<- DR_func(x=parmframeconc$concentration, 
                                     y0=parmframeconc$y0, 
                                     x0=parmframeconc$x0, 
                                     n=parmframeconc$n)* scale_factor  
      }
      parmframeconc$concentration<-concframe$x[k]
      predframe_median<- aggregate(Pred~concentration+ID, data=parmframeconc, quantile, prob=0.5)
      predframe_975<- aggregate(Pred~concentration+ID, data=parmframeconc, quantile, prob=0.975)
      predframe_25<- aggregate(Pred~concentration+ID, data=parmframeconc, quantile, prob=0.025)
      Prediction_Frame<- rbind(Prediction_Frame, cbind(predframe_median, predframe_25$Pred,predframe_975$Pred))
    }
    names(Prediction_Frame)[2]<- "Individual"
    names(Prediction_Frame)[3]<- "p50"
    names(Prediction_Frame)[4]<- "p2.5"
    names(Prediction_Frame)[5]<- "p97.5"
    return(Prediction_Frame)
  }
}


## Four Tiered Filter ##
chemfilter<- function(files=files){
  DRDat<- data.frame(cell, x, ys)
  DRCount<- data.frame()
  for(i in 1:42) {
    DRtemp<- subset(DRDat, cell == i)
    DPCount<- (length(unique(DRtemp$x)) - 1)
    if(DPCount < 3){
      DP.3 = FALSE
    } else {
      DP.3 = TRUE
    }
    DPtemp<- data.frame(i, DP.3)
    DRCount<- rbind(DRCount, DPtemp)
  }
  DP3<- sum(DRCount$DP.3 == TRUE)
  if(DP3 < 20){
    filter.res = "Fail"
    return(filter.res)
  } else {
  sigma<-as.matrix(rbind((eval(parse(text=files[1])))[2001,188],(eval(parse(text=files[2])))[2001,188], (eval(parse(text=files[3])))[2001,188], (eval(parse(text=files[4])))[2001,188]))
  sigma<- median(sigma, na.rm = T)
  EC05_med<-as.matrix(t(rbind((eval(parse(text=files[1])))[2:2001,622],(eval(parse(text=files[2])))[2:2001,622], (eval(parse(text=files[3])))[2:2001,622], (eval(parse(text=files[4])))[2:2001,622])))
  ECmed<-median(EC05_med, na.rm = T)
  EC05.05<-quantile(EC05_med, 0.05, na.rm = T)
  EC05.95<- quantile(EC05_med, 0.95, na.rm = T)
  Point.frame<- data.frame(x, ys)
  colnames(Point.frame)<- c("Concentration", "Response")
  if (sigma >= 0.1){
    sigma.filter<- FALSE
  } else {
    sigma.filter<-  TRUE
  }
  if (sigma.filter == TRUE){
    if (ECmed < max(Point.frame$Concentration)){
      EC05.ratio<- EC05.95/EC05.05
      if (EC05.ratio < 100){
        filter.res<- "Pass"
      } else {
        filter.res<- "Fail"
      }
    } else {
      filter.res<- "Fail"
    }
  } else {
    
    filter.res<- "Fail"
    
  }
  return(filter.res)
  }
}

ds.samplefilter <- function(samples = samples){
  sigma<- median(samples$sigma_y, na.rm = T)
  if(endpoint == "Peak_Amplitude_dn" | endpoint == "Peak_Amplitude_zero" | endpoint == "Peak_Frequency_zero"){
    EC_matrix <- samples$ec95_median
  } else {
    EC_matrix <- samples$ec05_median
  }
  ECmed<-median(EC_matrix, na.rm = T)
  EC.05<-quantile(EC_matrix, 0.05, na.rm = T)
  EC.95<- quantile(EC_matrix, 0.95, na.rm = T)
  Point.frame<- data.frame(x, ys)
  colnames(Point.frame)<- c("Concentration", "Response")
  if (sigma <= 0.1 | endpoint == "Peak_Amplitude_dn"){
    sigma.filter<- TRUE
  } else {
    sigma.filter<-  FALSE
  }
  if (sigma.filter == TRUE){
    if (ECmed < max(Point.frame$Concentration)){
      EC.ratio<- EC.95/EC.05
      if (EC.ratio < 100){
        filter.res<- "Pass"
      } else {
        filter.res<- "Fail"
      }
    } else {
      filter.res<- "Fail"
    }
  } else {
      filter.res<- "Fail"
    }
    return(filter.res)
  }

samplefilter <- function(samples = samples){
  if(Rdat$MaxRhat >= 1.2){
    filter.res <- "Fail"
    return(filter.res)
  }
  DRDat<- data.frame(cell, x, ys)
  DRCount<- data.frame()
  for(i in 1:42) {
    DRtemp<- subset(DRDat, cell == i)
    DPCount<- (length(unique(DRtemp$x)) - 1)
    if(DPCount < 3){
      DP.3 = FALSE
    } else {
      DP.3 = TRUE
    }
    DPtemp<- data.frame(i, DP.3)
    DRCount<- rbind(DRCount, DPtemp)
  }
  DP3<- sum(DRCount$DP.3 == TRUE)
  if(DP3 < 20){
    filter.res = "Fail"
    return(filter.res)
  } else {
    sigma<- median(samples$sigma_y, na.rm = T)
    if(endpoint == "Peak_Frequency_zero"){
      EC_matrix <- samples$ec95_median
    } else if (endpoint == "Total_Cells_zero"){
      EC_matrix <- samples$ec10_median
      #cat("\n", "EC10", "\n")
    } else {
      EC_matrix <- samples$ec05_median
    }
    ECmed<-median(EC_matrix, na.rm = T)
    EC.05<-quantile(EC_matrix, 0.05, na.rm = T)
    EC.95<- quantile(EC_matrix, 0.95, na.rm = T)
    Point.frame<- data.frame(x, ys)
    colnames(Point.frame)<- c("Concentration", "Response")
    if (sigma <= 0.1 | endpoint == "Peak_Amplitude_dn"){
      sigma.filter<- TRUE
    } else {
      sigma.filter<-  FALSE
    }
    if (sigma.filter == TRUE){
      if (ECmed < max(Point.frame$Concentration)){
        EC.ratio<- EC.95/EC.05
        if (EC.ratio < 100){
          filter.res<- "Pass"
        } else {
          filter.res<- "Fail"
        }
      } else {
        filter.res<- "Fail"
      }
    } else {
      
      filter.res<- "Fail"
      
    }
    return(filter.res)
  }
}

samplefilter.logic <- function(samples = samples){
  if(Rdat$MaxRhat >= 1.2){
    filter.res <- F
    return(filter.res)
  }
  DRDat<- data.frame(cell, x, ys)
  DRCount<- data.frame()
  for(i in 1:42) {
    DRtemp<- subset(DRDat, cell == i)
    DPCount<- (length(unique(DRtemp$x)) - 1)
    if(DPCount < 3){
      DP.3 = FALSE
    } else {
      DP.3 = TRUE
    }
    DPtemp<- data.frame(i, DP.3)
    DRCount<- rbind(DRCount, DPtemp)
  }
  DP3<- sum(DRCount$DP.3 == TRUE)
  if(DP3 < 20){
    filter.res = F
    return(filter.res)
  } else {
    sigma<- median(samples$sigma_y, na.rm = T)
    if(endpoint == "Peak_Frequency_zero"){
      EC_matrix <- samples$ec95_median
    } else if (endpoint == "Total_Cells_zero"){
      EC_matrix <- samples$ec10_median
      #cat("\n", "EC10", "\n")
    } else {
      EC_matrix <- samples$ec05_median
    }
    ECmed<-median(EC_matrix, na.rm = T)
    EC.05<-quantile(EC_matrix, 0.05, na.rm = T)
    EC.95<- quantile(EC_matrix, 0.95, na.rm = T)
    Point.frame<- data.frame(x, ys)
    colnames(Point.frame)<- c("Concentration", "Response")
    if (sigma <= 0.1 | endpoint == "Peak_Amplitude_dn"){
      sigma.filter<- TRUE
    } else {
      sigma.filter<-  FALSE
    }
    if (sigma.filter == TRUE){
      if (ECmed < max(Point.frame$Concentration)){
        EC.ratio<- EC.95/EC.05
        if (EC.ratio < 100){
          filter.res<- T
        } else {
          filter.res<- F
        }
      } else {
        filter.res<- F
      }
    } else {
      
      filter.res<- F
      
    }
    return(filter.res)
  }
}

samplefilter.partial<- function(samples = samples){
  if(Rdat$MaxRhat >= 1.2){
      filter.res <- "Fail"
      return(filter.res)
    }
  sigma<- median(samples$sigma_y, na.rm = T)
  if (sigma <= 0.1 | endpoint == "Peak_Amplitude_dn"){
    sigma.filter<- TRUE
  } else {
    sigma.filter<-  FALSE
  }
  if(endpoint == "Peak_Frequency_zero"){
    EC_matrix <- samples$ec95_median
  } else if (endpoint == "Total_Cells_zero"){
    EC_matrix <- samples$ec10_median
    #cat("\n", "EC10", "\n")
  } else {
    EC_matrix <- samples$ec05_median
  }
  ECmed<-median(EC_matrix, na.rm = T)
  #EC.05<-quantile(EC_matrix, 0.05, na.rm = T)
  #EC.95<- quantile(EC_matrix, 0.95, na.rm = T)
  Point.frame<- data.frame(x, ys)
  colnames(Point.frame)<- c("Concentration", "Response")
  if (sigma.filter == TRUE){
    if (ECmed < max(Point.frame$Concentration)){
      filter.res <- "Pass"
    } else {
      filter.res <- "Fail"
    }
  } else {
    filter.res <- "Fail"
  }
  
  return(filter.res)
  
}

samplefilter.partial.logic<- function(samples = samples){
  if(Rdat$MaxRhat >= 1.2){
    filter.res <- F
    return(filter.res)
  }
  sigma<- median(samples$sigma_y, na.rm = T)
  if (sigma <= 0.1 | endpoint == "Peak_Amplitude_dn"){
    sigma.filter<- TRUE
  } else {
    sigma.filter<-  FALSE
  }
  if(endpoint == "Peak_Frequency_zero"){
    EC_matrix <- samples$ec95_median
  } else if (endpoint == "Total_Cells_zero"){
    EC_matrix <- samples$ec10_median
    #cat("\n", "EC10", "\n")
  } else {
    EC_matrix <- samples$ec05_median
  }
  ECmed<-median(EC_matrix, na.rm = T)
  #EC.05<-quantile(EC_matrix, 0.05, na.rm = T)
  #EC.95<- quantile(EC_matrix, 0.95, na.rm = T)
  Point.frame<- data.frame(x, ys)
  colnames(Point.frame)<- c("Concentration", "Response")
  if (sigma.filter == TRUE){
    if (ECmed < max(Point.frame$Concentration)){
      filter.res <- T
    } else {
      filter.res <- F
    }
  } else {
    filter.res <- F
  }
  
  return(filter.res)
  
}

samplefilter_1000chems <- function(samples = samples){
  Point.frame<- data.frame(x, ys)
  colnames(Point.frame)<- c("Concentration", "Response")
  sigma<- median(samples$sigma_y, na.rm = T)
  if(endpoint == "Peak_Frequency_zero"){
  EC_matrix <- samples$ec95_median
  } else if (endpoint == "Total_Cells_zero"){
  EC_matrix <- samples$ec10_median
    #cat("\n", "EC10", "\n")
  } else {
    EC_matrix <- samples$ec05_median
  }
  ECmed<-median(EC_matrix, na.rm = T)
  if (sigma <= 0.1 | endpoint == "Peak_Amplitude_dn"){
    sigma.filter<- TRUE
  } else {
    sigma.filter<- FALSE
  }
  if (ECmed < max(Point.frame$Concentration)){
        EC.filter<- TRUE
      } else {
        EC.filter<- FALSE
      }
  
  if(sigma.filter == T & EC.filter == T){
    
    filter.res <- "Pass"
    
  } else {
    
    filter.res <- "Fail"
    
  }
  
  return(filter.res)
}  

indiv.samplefilter <- function(samples = samples){
  if(Rdat$MaxRhat >= 1.2){
    conv <- "Fail"
  } else if (Rdat$MaxRhat < 1.2){
    conv <- "Pass"
  }
  DRDat<- data.frame(cell, x, ys)
  DRCount<- data.frame()
  for(i in 1:42) {
    DRtemp<- subset(DRDat, cell == i)
    DPCount<- (length(unique(DRtemp$x)) - 1)
    if(DPCount < 3){
      DP.3 = FALSE
    } else {
      DP.3 = TRUE
    }
    DPtemp<- data.frame(i, DP.3)
    DRCount<- rbind(DRCount, DPtemp)
  }
  DP3<- sum(DRCount$DP.3 == TRUE)
  if(DP3 < 20){
    filter.1 <- "Fail"
  } else {
    filter.1 <- "Pass"
  }
    sigma<- median(samples$sigma_y, na.rm = T)
    if(endpoint == "Peak_Frequency_zero"){
      EC_matrix <- samples$ec95_median
    } else if (endpoint == "Total_Cells_zero"){
      EC_matrix <- samples$ec10_median
    } else {
      EC_matrix <- samples$ec05_median
    }
    ECmed<-median(EC_matrix, na.rm = T)
    EC.05<-quantile(EC_matrix, 0.05, na.rm = T)
    EC.95<- quantile(EC_matrix, 0.95, na.rm = T)
    Point.frame<- data.frame(x, ys)
    colnames(Point.frame)<- c("Concentration", "Response")
    if (sigma <= 0.1 | endpoint == "Peak_Amplitude_dn"){
      filter.2 <- "Pass"
    } else {
      filter.2<-"Fail"
    }
    if (ECmed < max(Point.frame$Concentration)){
      filter.3 <- "Pass"
    } else {
      filter.3 <- "Fail"
    }
    EC.ratio<- EC.95/EC.05
        if (EC.ratio < 100){
          filter.4 <- "Pass"
        } else {
          filter.4<- "Fail"
        }
    if(conv == "Pass" & filter.1 == "Pass" & filter.2 == "Pass" & filter.3 == "Pass" & filter.4 == "Pass"){
      All.filter <- "Pass"
    } else {
      All.filter <- "Fail"
    }
    
    filter.return <- data.frame(convergence = conv, filter.1 = filter.1 , filter.2 = filter.2, filter.3 = filter.3, filter.4 = filter.4, All.filter = All.filter)
    return(filter.return)
}

indiv.samplefilter.partial <- function(samples = samples){
  if(Rdat$MaxRhat >= 1.2){
    conv <- "Fail"
  } else if (Rdat$MaxRhat < 1.2){
    conv <- "Pass"
  }
  
  sigma<- median(samples$sigma_y, na.rm = T)
  if(endpoint == "Peak_Frequency_zero"){
    EC_matrix <- samples$ec95_median
  } else if (endpoint == "Total_Cells_zero"){
    EC_matrix <- samples$ec10_median
  } else {
    EC_matrix <- samples$ec05_median
  }
  ECmed<-median(EC_matrix, na.rm = T)
  Point.frame<- data.frame(x, ys)
  colnames(Point.frame)<- c("Concentration", "Response")
  if (sigma <= 0.1 | endpoint == "Peak_Amplitude_dn"){
    filter.1 <- "Pass"
  } else {
    filter.1<-"Fail"
  }
  if (ECmed < max(Point.frame$Concentration)){
    filter.2 <- "Pass"
  } else {
    filter.2 <- "Fail"
  }
  if(conv == "Pass" & filter.1 == "Pass" & filter.2 == "Pass"){
    All.filter <- "Pass"
  } else {
    All.filter <- "Fail"
  }
  
  filter.return <- data.frame(convergence = conv, filter.1 = filter.1 , filter.2 = filter.2, All.filter = All.filter)
  return(filter.return)
}


indiv.samplefilter.1000chem <- function(samples = samples){
  sigma<- median(samples$sigma_y, na.rm = T)
  if(endpoint == "Peak_Frequency_zero"){
    EC_matrix <- samples$ec95_median
  } else if (endpoint == "Total_Cells_zero"){
    EC_matrix <- samples$ec10_median
  } else {
    EC_matrix <- samples$ec05_median
  }
  ECmed<-median(EC_matrix, na.rm = T)
  EC.05<-quantile(EC_matrix, 0.05, na.rm = T)
  EC.95<- quantile(EC_matrix, 0.95, na.rm = T)
  Point.frame<- data.frame(x, ys)
  colnames(Point.frame)<- c("Concentration", "Response")
  if (sigma <= 0.1 | endpoint == "Peak_Amplitude_dn"){
    filter.1 <- "Pass"
  } else {
    filter.1<-"Fail"
  }
  if (ECmed < max(Point.frame$Concentration)){
    filter.2 <- "Pass"
  } else {
    filter.2 <- "Fail"
  }
  #EC.ratio<- EC.95/EC.05
  #if (EC.ratio < 100){
    #filter.3 <- "Pass"
  #} else {
    #filter.3<- "Fail"
  #}
  if(filter.1 == "Pass" & filter.2 == "Pass"){
    All.filter <- "Pass"
  } else {
    All.filter <- "Fail"
  }
  filter.return <- data.frame(filter.1 = filter.1 , filter.2 = filter.2, All.filter = All.filter)
  return(filter.return)
}



##Generate Summary Statistics##
SumStatGen<- function(files = files, CI = 90){
  ##Define CI
  if(is.na(CI) == TRUE){
    warning("\n", "ERROR: No Confidence Interval Defined", "\n")
    stop()
  }
  if(CI == 90){
    low.CI = 0.05
    high.CI = 0.95
  } else if (CI == 95){
    low.CI = 0.025
    high.CI = 0.975
  } else {
    warning("\n", "ERROR: Select either a 90% CI or a 95% CI")
    stop()
  }
  ##Treatment Data Point Number
  DRDat<- data.frame(cell, x, ys)
  DPTot = nrow(subset(DRDat, x > 0))
  DRCount<- data.frame()
  for(i in 1:42) {
    DRtemp<- subset(DRDat, cell == i)
    DPCount<- (length(unique(DRtemp$x)) - 1)
    if(DPCount < 2){
      DP.2 = FALSE
      DP.3 = FALSE
    }
    if(DPCount >= 2){
      DP.2 = TRUE
    } else { DP.2 = FALSE}
    if(DPCount >= 3){
      DP.3 = TRUE
    } else { DP.3 = FALSE}
    DPtemp<- data.frame(i, DP.2, DP.3)
    DRCount<- rbind(DRCount, DPtemp)
  }
  DP2<- sum(DRCount$DP.2 == TRUE)
  DP3<- sum(DRCount$DP.3 == TRUE)
  #DPCount<- data.frame(DPTot = DPTot, Indiv2C = DP2, Indiv3C = DP3)
  DPCount<- data.frame(Indiv3C = DP3)
  ## Load in hyperparms
  m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001,c(8,12)],(eval(parse(text=files[2])))[2:2001, c(8,12)], (eval(parse(text=files[3])))[2:2001,c(8,12)], (eval(parse(text=files[4])))[2:2001,c(8,12)]))
  m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, c(9,13)], (eval(parse(text=files[2])))[2:2001, c(9,13)], (eval(parse(text=files[3])))[2:2001, c(9,13)], (eval(parse(text=files[4])))[2:2001, c(9,13)]))
  m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, c(10,14)], (eval(parse(text=files[2])))[2:2001,c(10,14)], (eval(parse(text=files[3])))[2:2001, c(10,14)], (eval(parse(text=files[4])))[2:2001, c(10,14)]))
  m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, c(11,15)], (eval(parse(text=files[2])))[2:2001, c(11,15)], (eval(parse(text=files[3])))[2:2001, c(11,15)], (eval(parse(text=files[4])))[2:2001, c(11,15)]))
  ##TDVFs
  TDVFframe<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, c(702,704)], (eval(parse(text=files[2])))[2:2001, c(702,704)], (eval(parse(text=files[3])))[2:2001, c(702,704)], (eval(parse(text=files[4])))[2:2001, c(702,704)]))
  TDVF01.Med<- median(TDVFframe[,1], na.rm = T)
  TDVF01.05<- quantile(TDVFframe[,1], 0.05,na.rm = T)
  TDVF01.95<- quantile(TDVFframe[,1], 0.95,na.rm = T)
  TDVF05.Med<- median(TDVFframe[,2], na.rm = T)
  TDVF05.05<- quantile(TDVFframe[,2], 0.05,na.rm = T)
  TDVF05.95<- quantile(TDVFframe[,2], 0.95,na.rm = T)
  TDVF<- data.frame(TDVF01.Med = TDVF01.Med, TDVF01.05 = TDVF01.05, TDVF01.95 = TDVF01.95, TDVF05.Med = TDVF05.Med, TDVF05.05 = TDVF05.05, TDVF05.95 = TDVF05.95)
  ##Calculate CIs
  m_x0.CI<- data.frame(signif(median(m_x0[,1],na.rm = T), digits = 3), signif(quantile(m_x0[,1], low.CI,na.rm = T), digits = 3), signif(quantile(m_x0[,1], high.CI,na.rm = T), digits = 3))
  m_x0.diff<- data.frame(m_x0.CI[,3] - m_x0.CI[,2])
  m_x0.SD<- data.frame(signif(median(m_x0[,2],na.rm = T), digits = 3))
  m_Emax.CI<- data.frame(signif(median(m_Emax[,1],na.rm = T), digits = 3), signif(quantile(m_Emax[,1], low.CI,na.rm = T), digits = 3), signif(quantile(m_Emax[,1], high.CI,na.rm = T), digits = 3))
  m_Emax.SD<- data.frame(signif(median(m_Emax[,2]), digits = 3))
  m_n.CI<- data.frame(signif(median(m_n[,1],na.rm = T), digits = 3), signif(quantile(m_n[,1], low.CI,na.rm = T), digits = 3), signif(quantile(m_n[,1], high.CI,na.rm = T), digits = 3))
  m_n.SD<- data.frame(signif(median(m_n[,2],na.rm = T), digits = 3))
  sigma<-as.matrix(rbind((eval(parse(text=files[1])))[2:2001,188],(eval(parse(text=files[2])))[2:2001,188], (eval(parse(text=files[3])))[2:2001,188], (eval(parse(text=files[4])))[2:2001,188]))
  sigma<- data.frame(signif(median(sigma,na.rm = T), digits = 3), signif(quantile(sigma, low.CI,na.rm = T), digits=3), signif(quantile(sigma, high.CI,na.rm = T), digits = 3))
  EC05_med<-as.matrix(t(rbind((eval(parse(text=files[1])))[2:2001,622],(eval(parse(text=files[2])))[2:2001,622], (eval(parse(text=files[3])))[2:2001,622], (eval(parse(text=files[4])))[2:2001,622])))
  ECmed<-median(EC05_med,na.rm = T)
  EC05.05<-quantile(EC05_med, 0.05,na.rm = T)
  EC05.95<- quantile(EC05_med, 0.95,na.rm = T)
  EC05.ratio<- EC05.95/EC05.05
  EC05.frame<- data.frame(EC05_med = ECmed, EC05.05 = EC05.05, EC05.95 = EC05.95, EC05.ratio = EC05.ratio)
  filter.res<- chemfilter(files = files)
  ##Format Confidence Frame
  Confidence.frame<- cbind(chemical.name, filter.res, DPCount, m_x0.CI, m_x0.diff, m_x0.SD, m_Emax.CI, m_Emax.SD, m_n.CI, m_n.SD, sigma, EC05.frame, TDVF )
  colnames(Confidence.frame)[c(1,2,4:19)]<- c("Chemical", "Filter", "x0.Med", "x0.5", "x0.95", "x0.diff" ,"x0.SD", "Emax.Med", "Emax.5", "Emax.95", "Emax.SD", "n.Med", "n.5", "n.95", "n.SD", "Sigma.Med", "Sigma.5", "Sigma.95")
  return(Confidence.frame)
}

SumStatGenSamples <- function(samples = samples, CI = 90){
  dir <- ep.compass(endpoint = endpoint, return.dir = T)
  ##Define CI
  if(is.na(CI) == TRUE){
    warning("\n", "ERROR: No Confidence Interval Defined", "\n")
    stop()
  }
  if(CI == 90){
    low.CI = 0.05
    high.CI = 0.95
  } else if (CI == 95){
    low.CI = 0.025
    high.CI = 0.975
  } else {
    warning("\n", "ERROR: Select either a 90% CI or a 95% CI")
    stop()
  }
  ##Convergence Check
  Max.Rhat <- Rdat$MaxRhat
  ##Treatment Data Point Number
  DRDat<- data.frame(cell, x, ys)
  DPTot = nrow(subset(DRDat, x > 0))
  DRCount<- data.frame()
  for(i in 1:42) {
    DRtemp<- subset(DRDat, cell == i)
    DPCount<- (length(unique(DRtemp$x)) - 1)
    if(DPCount < 2){
      DP.2 = FALSE
      DP.3 = FALSE
    }
    if(DPCount >= 2){
      DP.2 = TRUE
    } else { DP.2 = FALSE}
    if(DPCount >= 3){
      DP.3 = TRUE
    } else { DP.3 = FALSE}
    DPtemp<- data.frame(i, DP.2, DP.3)
    DRCount<- rbind(DRCount, DPtemp)
  }
  DP2<- sum(DRCount$DP.2 == TRUE)
  DP3<- sum(DRCount$DP.3 == TRUE)
  #DPCount<- data.frame(DPTot = DPTot, Indiv2C = DP2, Indiv3C = DP3)
  DPCount<- data.frame(Indiv3C = DP3)
  ## Load in hyperparms
  m_y0<-rbind(samples$m_y0, samples$sd_y0)
  m_x0<-rbind(samples$m_x0, samples$sd_x0)
  if(dir != 0){
  m_Emax<- rbind(samples$m_Emax, samples$sd_Emax)
  }
  m_n<-rbind(samples$m_n, samples$sd_n)
  ##TDVFs
  if(endpoint == "Peak_Amplitude_dn" | endpoint == "Peak_Amplitude_zero" | endpoint == "Peak_Frequency_zero"){
    TDVFframe<-samples$ec95_quant_ratios.3
  } else if(endpoint == "Total_Cells_zero"){
    TDVFframe <- samples$ec10_quant_ratios.3
  } else {
    TDVFframe<-samples$ec05_quant_ratios.3
  }
  TDVF05.Med<- median(TDVFframe, na.rm = T)
  TDVF05.05<- quantile(TDVFframe, 0.05,na.rm = T)
  TDVF05.95<- quantile(TDVFframe, 0.95,na.rm = T)
  TDVF<- data.frame(TDVF05.Med = TDVF05.Med, TDVF05.05 = TDVF05.05, TDVF05.95 = TDVF05.95)
  ##Calculate CIs
  m_x0.CI<- data.frame(signif(median(m_x0[,1],na.rm = T), digits = 3), signif(quantile(m_x0[,1], low.CI,na.rm = T), digits = 3), signif(quantile(m_x0[,1], high.CI,na.rm = T), digits = 3))
  m_x0.diff<- data.frame(m_x0.CI[,3] - m_x0.CI[,2])
  m_x0.SD<- data.frame(signif(median(m_x0[,2],na.rm = T), digits = 3))
  if(dir != 0){
  m_Emax.CI<- data.frame(signif(median(m_Emax[,1],na.rm = T), digits = 3), signif(quantile(m_Emax[,1], low.CI,na.rm = T), digits = 3), signif(quantile(m_Emax[,1], high.CI,na.rm = T), digits = 3))
  m_Emax.SD<- data.frame(signif(median(m_Emax[,2]), digits = 3))
  }
  m_n.CI<- data.frame(signif(median(m_n[,1],na.rm = T), digits = 3), signif(quantile(m_n[,1], low.CI,na.rm = T), digits = 3), signif(quantile(m_n[,1], high.CI,na.rm = T), digits = 3))
  m_n.SD<- data.frame(signif(median(m_n[,2],na.rm = T), digits = 3))
  sigma<- samples$sigma_y
  sigma<- data.frame(signif(median(sigma,na.rm = T), digits = 3), signif(quantile(sigma, low.CI,na.rm = T), digits=3), signif(quantile(sigma, high.CI,na.rm = T), digits = 3))
  if(endpoint == "Peak_Amplitude_dn" | endpoint == "Peak_Amplitude_zero" | endpoint == "Peak_Frequency_zero") {
    EC_matrix <- samples$ec95_median
    EC.name <- "EC95"
  } else if(endpoint == "Total_Cells_zero"){
    EC_matrix <- samples$ec10_median
    EC.name <- "EC10"
  } else {
    EC_matrix <- samples$ec05_median
    EC.name <- "EC05"
  }
  ECmed<-median(EC_matrix,na.rm = T)
  EC.05<-quantile(EC_matrix, 0.05,na.rm = T)
  EC.95<- quantile(EC_matrix, 0.95,na.rm = T)
  EC.ratio<- EC.95/EC.05
  EC.frame<- data.frame(ECmed, EC.05, EC.95, EC.ratio)
  colnames(EC.frame) <- paste(EC.name, c("_med", ".05", ".95", ".ratio"), sep = "")
  filter.res<- samplefilter(samples = samples)
  ##Format Confidence Frame
  if(dir != 0){
  Confidence.frame<- cbind(chemical.name, filter.res, Max.Rhat, DPCount, m_x0.CI, m_x0.diff, m_x0.SD, m_Emax.CI, m_Emax.SD, m_n.CI, m_n.SD, sigma, EC.frame, TDVF )
  colnames(Confidence.frame)[c(1,2,5:20)]<- c("Chemical", "Filter", "x0.Med", "x0.5", "x0.95", "x0.diff" ,"x0.SD", "Emax.Med", "Emax.5", "Emax.95", "Emax.SD", "n.Med", "n.5", "n.95", "n.SD", "Sigma.Med", "Sigma.5", "Sigma.95")
  } else {
    Confidence.frame<- cbind(chemical.name, filter.res, Max.Rhat, DPCount, m_x0.CI, m_x0.diff, m_x0.SD, m_n.CI, m_n.SD, sigma, EC.frame, TDVF )
    colnames(Confidence.frame)[c(1,2,5:16)]<- c("Chemical", "Filter", "x0.Med", "x0.5", "x0.95", "x0.diff" ,"x0.SD", "n.Med", "n.5", "n.95", "n.SD", "Sigma.Med", "Sigma.5", "Sigma.95")
  }
  return(Confidence.frame)
}


## Normalization of data with y0 ##

resp.y0norm<- function(files = files, name = chemical.name, num = chemical.num){
  ## Set up Point Frame
  Point.frame<- data.frame(cell,x,ys)
  colnames(Point.frame)<- c("Cell.Line", "Concentration", "Response")
  ## set up y0 frame ##
  y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, 186:227], (eval(parse(text=files[2])))[2:2001, 186:227], (eval(parse(text=files[3])))[2:2001, 186:227], (eval(parse(text=files[4])))[2:2001, 186:227]))
  colnames(y0)<- (1:42)
  rownames(y0)<- c(1:8000)
  y0<- sapply(y0, median)
  y0<- data.frame(Cell.Line = 1:42, y0_med = y0)
  ## merge point frame and y0 frame 
  Point.frame.y0<- merge(Point.frame, y0, by = "Cell.Line")
  Point.frame.y0$NormalizedResponse<- Point.frame.y0$Response/ Point.frame.y0$y0_med
  Point.frame.y0$Chemical.name<- name
  Point.frame.y0$Chemical.num<- num
  return(Point.frame.y0)
}

resp.y0norm.samples <- function(samples = samples, name = chemical.name, num = chemical.num){
  ## Set up Point Frame
  Point.frame<- data.frame(cell,x,ys)
  colnames(Point.frame)<- c("Cell.Line", "Concentration", "Response")
  ## set up y0 frame ##
  if(endpoint == "Peak_Frequency_zero" | endpoint == "Total_Cells_zero"){
    y0 <- samples[,142:183]  
  } else {
    y0 <- samples[,186:227]
  }
  colnames(y0)<- (1:42)
  rownames(y0)<- c(1:1000)
  y0<- sapply(y0, median)
  y0<- data.frame(Cell.Line = 1:42, y0_med = y0)
  ## merge point frame and y0 frame 
  Point.frame.y0<- merge(Point.frame, y0, by = "Cell.Line")
  Point.frame.y0$NormalizedResponse<- Point.frame.y0$Response/ Point.frame.y0$y0_med
  Point.frame.y0$Chemical.name<- name
  Point.frame.y0$Chemical.num<- num
  return(Point.frame.y0)
}

## Simulation Frame Setup ##

SimFrameGen<- function(endpoint = endpoint, files = files, PredType = "Pop", Scaling = FALSE, y0Norm = FALSE, returnparms = FALSE, prob.uppper = 0.95, prob.lower = 0.05){
    concframe<-data.frame(x=c(1e-6,10^(seq(-3,2,0.05))))
    ### Pop median prediction
    m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001,c(8,12)],(eval(parse(text=files[2])))[2:2001, c(8,12)], (eval(parse(text=files[3])))[2:2001,c(8,12)], (eval(parse(text=files[4])))[2:2001,c(8,12)]))
    ##### x0
    m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, c(9,13)], (eval(parse(text=files[2])))[2:2001, c(9,13)], (eval(parse(text=files[3])))[2:2001, c(9,13)], (eval(parse(text=files[4])))[2:2001, c(9,13)]))
    ##### Emax
    m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, c(10,14)], (eval(parse(text=files[2])))[2:2001,c(10,14)], (eval(parse(text=files[3])))[2:2001, c(10,14)], (eval(parse(text=files[4])))[2:2001, c(10,14)]))
    ##### n
    m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, c(11,15)], (eval(parse(text=files[2])))[2:2001, c(11,15)], (eval(parse(text=files[3])))[2:2001, c(11,15)], (eval(parse(text=files[4])))[2:2001, c(11,15)]))
    if(PredType == "Random"){
    dir<- ep.compass(endpoint, return.dir = T)
    ### Pop median prediction
    m_y0$z_score_sim<- rnorm(8000,0,1)
    simframe<- as.data.frame(exp(m_y0$m_y0+(m_y0$sd_y0*m_y0$z_score_sim)))
    rownames(simframe)<- c(1:8000)
    colnames(simframe)<- "sim_y0"
    ##### x0
    m_x0$z_score_sim<- rnorm(8000,0,1)
    simframe$sim_x0<-exp(m_x0$m_x0+(m_x0$sd_x0*m_x0$z_score_sim))
    ##### Emax
    m_Emax$z_score_sim<- rnorm(8000,0,1)
    if (dir == -1){
      simframe$sim_Emax<- inv.logit(m_Emax$m_Emax+(m_Emax$sd_Emax*m_Emax$z_score_sim))
    } else if (dir == 1){
      simframe$sim_Emax<- exp(m_Emax$m_Emax+(m_Emax$sd_Emax*m_Emax$z_score_sim))
    }
    ##### n
    m_n$z_score_sim<- rnorm(8000,0,1)
    simframe$sim_n<-exp(m_n$m_n+(m_n$sd_n*m_n$z_score_sim))
    ## set up simframeconc
    simframeconc<- data.frame()
    for (i in 1:nrow(concframe)) {
      simframe$concentration<-concframe$x[i]
      simframeconc<- rbind(simframeconc, simframe)
    }
    if(y0Norm == F){simframeconc$PredFold<- DR_func(x=simframeconc$concentration, y0=simframeconc$sim_y0,x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n)
    }else{simframeconc$PredFold<- DR_func(x=simframeconc$concentration, y0=1,x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n)}
    if(Scaling == TRUE){simframeconc$PredFold<- (simframeconc$PredFold)*scale_factor}
    ##### Quantiles
    simframe_median_fold<- aggregate(PredFold~concentration, data=simframeconc, quantile, prob=0.5)
    simframe_975_fold<- aggregate(PredFold~concentration, data=simframeconc, quantile, prob=0.95)
    simframe_25_fold<- aggregate(PredFold~concentration, data=simframeconc, quantile, prob=0.05)
    Simulation_Frame_Random<- simframe_median_fold
    Simulation_Frame_Random$p95<- simframe_975_fold$PredFold
    Simulation_Frame_Random$p05<- simframe_25_fold$PredFold
    names(Simulation_Frame_Random)[2]<- "p50"
    if(returnparms == FALSE){
    return(Simulation_Frame_Random)
    } else {return(simframe[,1:4])}
    } else if (PredType == "Pop"){
      dir<- ep.compass(endpoint, return.dir = T)
      
      #### Population median
      if(dir == 1){
      simframe_zero<-data.frame(y0=(exp(m_y0$m_y0)), x0=(exp(m_x0$m_x0)), Emax=(exp(m_Emax$m_Emax)), n=(exp(m_n$m_n)))
      } else if (dir == -1){
        simframe_zero<-data.frame(y0=(exp(m_y0$m_y0)), x0=(exp(m_x0$m_x0)), Emax=(inv.logit(m_Emax$m_Emax)), n=(exp(m_n$m_n)))
      }
      simframe_zero_conc<- data.frame() 
      for (i in 1:nrow(concframe)) {
        simframe_zero$concentration<-concframe$x[i]
        simframe_zero_conc<- rbind(simframe_zero_conc, simframe_zero)
      }
      #### Fold change - Population median
      if(y0Norm == F){simframe_zero_conc$PredFold<- DR_func(x=simframe_zero_conc$concentration, y0=simframe_zero_conc$y0, x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n) 
      } else {simframe_zero_conc$PredFold<- DR_func(x=simframe_zero_conc$concentration, y0=1, x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n)}
      if(Scaling == TRUE){simframe_zero_conc$PredFold<- (simframe_zero_conc$PredFold)*scale_factor}
      ##### Quantiles
      simframe_zero_median_fold<- (aggregate(PredFold~concentration, data=simframe_zero_conc, quantile, prob=0.5))
      simframe_zero_975_fold<- (aggregate(PredFold~concentration, data=simframe_zero_conc, quantile, prob= prob.upper))
      simframe_zero_25_fold<- (aggregate(PredFold~concentration, data=simframe_zero_conc, quantile, prob= prob.lower))
      Simulation_Frame_Pop<- simframe_zero_median_fold
      Simulation_Frame_Pop$pupper<- simframe_zero_975_fold$PredFold
      Simulation_Frame_Pop$plower<- simframe_zero_25_fold$PredFold
      names(Simulation_Frame_Pop)[2]<- "p50"
      if(returnparms == FALSE){
      return(Simulation_Frame_Pop)
      } else {return(simframe_zero[,1:4])}
    }
}

SimFrameGenSamples <- function(endpoint = endpoint, samples = samples, PredType = "Pop", Scaling = FALSE, y0Norm = FALSE, 
                               returnparms = FALSE, n.samples = 1000, GetFreeConc = F, 
                               prob.upper = 0.95, prob.lower = 0.05, returnconcframe = F){
  DR_func <- ep.compass(endpoint = endpoint)
  dir <- ep.compass(endpoint = endpoint, return.dir = T)
  concframe<-data.frame(x=c(1e-6,10^(seq(-3,2,0.05))))
  ### Pop median prediction
  m_y0<- data.frame(m_y0 = samples$m_y0, sd_y0 = samples$sd_y0)
  ##### x0
  m_x0<- data.frame(m_x0 = samples$m_x0, sd_x0 = samples$sd_x0)
  if(dir != 0){
  ##### Emax
  m_Emax<- data.frame(m_Emax = samples$m_Emax, sd_Emax = samples$sd_Emax)
  }
  ##### n
  m_n<- data.frame(m_n = samples$m_n, sd_n = samples$sd_n)
  if(PredType == "Random"){
    dir<- ep.compass(endpoint, return.dir = T)
    ### Pop median prediction
    m_y0$z_score_sim<- rnorm(n.samples,0,1)
    simframe<- as.data.frame(exp(m_y0$m_y0+(m_y0$sd_y0*m_y0$z_score_sim)))
    rownames(simframe)<- c(1:n.samples)
    colnames(simframe)<- "sim_y0"
    ##### x0
    m_x0$z_score_sim<- rnorm(n.samples,0,1)
    simframe$sim_x0<-exp(m_x0$m_x0+(m_x0$sd_x0*m_x0$z_score_sim))
    if(dir != 0){
      ##### Emax
      m_Emax$z_score_sim<- rnorm(n.samples,0,1)
      if (dir == -1){
        simframe$sim_Emax<- inv.logit(m_Emax$m_Emax+(m_Emax$sd_Emax*m_Emax$z_score_sim))
      } else if (dir == 1){
        simframe$sim_Emax<- exp(m_Emax$m_Emax+(m_Emax$sd_Emax*m_Emax$z_score_sim))
      }
    }
    ##### n
    m_n$z_score_sim<- rnorm(n.samples,0,1)
    simframe$sim_n<-exp(m_n$m_n+(m_n$sd_n*m_n$z_score_sim))
    ## set up simframeconc
    simframeconc<- data.frame()
    for (i in 1:nrow(concframe)) {
      simframe$concentration<-concframe$x[i]
      if(GetFreeConc == T){
        simframe$freeconcentration<-concframe$x[i]*InVivoFree$FracFreeMedia
      }
      simframeconc<- rbind(simframeconc, simframe)
    }
    if(dir != 0){
      if(y0Norm == F){simframeconc$PredFold<- DR_func(x=simframeconc$concentration, y0=simframeconc$sim_y0,x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n)
      }else{simframeconc$PredFold<- DR_func(x=simframeconc$concentration, y0=1,x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n)}
      if(Scaling == TRUE){simframeconc$PredFold<- (simframeconc$PredFold)*scale_factor}
    } else {
      if(y0Norm == F){simframeconc$PredFold<- DR_func(x=simframeconc$concentration, y0=simframeconc$sim_y0,x0=simframeconc$sim_x0, n=simframeconc$sim_n)
      }else{simframeconc$PredFold<- DR_func(x=simframeconc$concentration, y0=1,x0=simframeconc$sim_x0, n=simframeconc$sim_n)}
      if(Scaling == TRUE){simframeconc$PredFold<- (simframeconc$PredFold)*scale_factor}
    }
    if(returnconcframe == T){
      return(simframeconc)
    }
    ##### Quantiles
    simframe_median_fold<- aggregate(PredFold~concentration, data=simframeconc, quantile, prob=0.5)
    simframe_975_fold<- aggregate(PredFold~concentration, data=simframeconc, quantile, prob= prob.upper)
    simframe_25_fold<- aggregate(PredFold~concentration, data=simframeconc, quantile, prob= prob.lower)
    Simulation_Frame_Random<- simframe_median_fold
    if(GetFreeConc == T){
      Simulation_Frame_Random$freeconcentration <- Simulation_Frame_Random$concentration*InVivoFree$FracFreeMedia
    }
    Simulation_Frame_Random$pupper<- simframe_975_fold$PredFold
    Simulation_Frame_Random$plower<- simframe_25_fold$PredFold
    names(Simulation_Frame_Random)[2]<- "p50"
    if(returnparms == FALSE){
      return(Simulation_Frame_Random)
    } else {return(simframe[,1:4])}
  } else if (PredType == "Pop"){
    dir<- ep.compass(endpoint, return.dir = T)
    #### Population median
    if(dir != 0){
      if(dir == 1){
        simframe_zero<-data.frame(y0=(exp(m_y0$m_y0)), x0=(exp(m_x0$m_x0)), Emax=(exp(m_Emax$m_Emax)), n=(exp(m_n$m_n)))
      } else if (dir == -1){
        simframe_zero<-data.frame(y0=(exp(m_y0$m_y0)), x0=(exp(m_x0$m_x0)), Emax=(inv.logit(m_Emax$m_Emax)), n=(exp(m_n$m_n)))
      }
    } else {
      simframe_zero<-data.frame(y0=(exp(m_y0$m_y0)), x0=(exp(m_x0$m_x0)), n=(exp(m_n$m_n)))
    }
    simframe_zero_conc<- data.frame() 
    for (i in 1:nrow(concframe)) {
      simframe_zero$concentration<-concframe$x[i]
      if(GetFreeConc == T){
        simframe_zero$freeconcentration<-concframe$x[i]*InVivoFree$FracFreeMedia
      }
      simframe_zero_conc<- rbind(simframe_zero_conc, simframe_zero)
    }
    #### Fold change - Population median
    if(dir != 0){
      if(y0Norm == F){simframe_zero_conc$PredFold<- DR_func(x=simframe_zero_conc$concentration, y0=simframe_zero_conc$y0, x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n) 
      } else {simframe_zero_conc$PredFold<- DR_func(x=simframe_zero_conc$concentration, y0=1, x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n)}
      if(Scaling == TRUE){simframe_zero_conc$PredFold<- (simframe_zero_conc$PredFold)*scale_factor}
    } else {
      if(y0Norm == F){simframe_zero_conc$PredFold<- DR_func(x=simframe_zero_conc$concentration, y0=simframe_zero_conc$y0, x0=simframe_zero_conc$x0, n=simframe_zero_conc$n) 
      } else {simframe_zero_conc$PredFold<- DR_func(x=simframe_zero_conc$concentration, y0=1, x0=simframe_zero_conc$x0, n=simframe_zero_conc$n)}
      if(Scaling == TRUE){simframe_zero_conc$PredFold<- (simframe_zero_conc$PredFold)*scale_factor}
    }
    if(returnconcframe == T){
      return(simframe_zero_conc)
    }
    ##### Quantiles
    simframe_zero_median_fold<- (aggregate(PredFold~concentration, data=simframe_zero_conc, quantile, prob=0.5))
    simframe_zero_975_fold<- (aggregate(PredFold~concentration, data=simframe_zero_conc, quantile, prob= prob.upper))
    simframe_zero_25_fold<- (aggregate(PredFold~concentration, data=simframe_zero_conc, quantile, prob= prob.lower))
    Simulation_Frame_Pop<- simframe_zero_median_fold
    if(GetFreeConc == T){
      Simulation_Frame_Pop$freeconcentration <- Simulation_Frame_Pop$concentration*InVivoFree$FracFreeMedia
    }
    Simulation_Frame_Pop$pupper<- simframe_zero_975_fold$PredFold
    Simulation_Frame_Pop$plower<- simframe_zero_25_fold$PredFold
    names(Simulation_Frame_Pop)[2]<- "p50"
    if(returnparms == FALSE){
      return(Simulation_Frame_Pop)
    } else {return(simframe_zero[,1:4])}
  }
}

SimFrameGen.Frac <- function(endpoint = endpoint, samples = samples, PredType = "Pop", Scaling = FALSE, y0Norm = FALSE, returnparms = FALSE, 
                             GetFreeConc = TRUE, prob.upper = 0.95, prob.lower = 0.05, returnconcframe = FALSE){
  DR_func <- DR_up_Frac
  dir <- ep.compass(endpoint = endpoint, return.dir = T)
  concframe<-data.frame(x=c(1e-6,10^(seq(-3,2,0.05))))
  ### Pop median prediction
  m_y0<- data.frame(m_y0 = samples$m_y0, sd_y0 = samples$sd_y0)
  ##### x0
  m_x0<- data.frame(m_x0 = samples$m_x0, sd_x0 = samples$sd_x0)
  if(dir != 0){
    ##### Emax
    m_Emax<- data.frame(m_Emax = samples$m_Emax, sd_Emax = samples$sd_Emax)
  }
  ##### n
  m_n<- data.frame(m_n = samples$m_n, sd_n = samples$sd_n)
  if(PredType == "Random"){
    dir<- ep.compass(endpoint, return.dir = T)
    ### Pop median prediction
    m_y0$z_score_sim<- rnorm(8000,0,1)
    simframe<- as.data.frame(exp(m_y0$m_y0+(m_y0$sd_y0*m_y0$z_score_sim)))
    rownames(simframe)<- c(1:8000)
    colnames(simframe)<- "sim_y0"
    ##### x0
    m_x0$z_score_sim<- rnorm(8000,0,1)
    simframe$sim_x0<-exp(m_x0$m_x0+(m_x0$sd_x0*m_x0$z_score_sim))
    if(dir != 0){
      ##### Emax
      m_Emax$z_score_sim<- rnorm(8000,0,1)
      if (dir == -1){
        simframe$sim_Emax<- inv.logit(m_Emax$m_Emax+(m_Emax$sd_Emax*m_Emax$z_score_sim))
      } else if (dir == 1){
        simframe$sim_Emax<- exp(m_Emax$m_Emax+(m_Emax$sd_Emax*m_Emax$z_score_sim))
      }
    }
    ##### n
    m_n$z_score_sim<- rnorm(8000,0,1)
    simframe$sim_n<-exp(m_n$m_n+(m_n$sd_n*m_n$z_score_sim))
    ## set up simframeconc
    simframeconc<- data.frame()
    for (i in 1:nrow(concframe)) {
      simframe$concentration<-concframe$x[i]
      if(GetFreeConc == T){
      simframe$freeconcentration<-concframe$x[i]*InVivoFree$FracFreeMedia
      }
      simframeconc<- rbind(simframeconc, simframe)
    }
    if(dir != 0){
      if(y0Norm == F){simframeconc$PredFrac<- DR_func(x=simframeconc$concentration,x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n)
      }else{simframeconc$PredFrac<- DR_func(x=simframeconc$concentration,x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n)}
      if(Scaling == TRUE){simframeconc$PredFrac<- (simframeconc$PredFrac)*scale_factor}
    } else {
      if(y0Norm == F){simframeconc$PredFrac<- DR_func(x=simframeconc$concentration,x0=simframeconc$sim_x0, n=simframeconc$sim_n)
      }else{simframeconc$PredFrac<- DR_func(x=simframeconc$concentration,x0=simframeconc$sim_x0, n=simframeconc$sim_n)}
      if(Scaling == TRUE){simframeconc$PredFrac<- (simframeconc$PredFrac)*scale_factor}
    }
    if(returnconcframe == T){
      return(simframeconc)
    }
    ##### Quantiles
    simframe_median_Frac<- aggregate(PredFrac~concentration, data=simframeconc, quantile, prob=0.5)
    simframe_975_Frac<- aggregate(PredFrac~concentration, data=simframeconc, quantile, prob= prob.upper)
    simframe_25_Frac<- aggregate(PredFrac~concentration, data=simframeconc, quantile, prob= prob.lower)
    Simulation_Frame_Random<- simframe_median_Frac
    if(GetFreeConc == T){
      Simulation_Frame_Random$freeconcentration <- Simulation_Frame_Random$concentration*InVivoFree$FracFreeMedia
    }
    Simulation_Frame_Random$pupper<- simframe_975_Frac$PredFrac
    Simulation_Frame_Random$plower<- simframe_25_Frac$PredFrac
    names(Simulation_Frame_Random)[2]<- "p50"
    if(returnparms == FALSE){
      return(Simulation_Frame_Random)
    } else {return(simframe[,1:4])}
  } else if (PredType == "Pop"){
    dir<- ep.compass(endpoint, return.dir = T)
    #### Population median
    if(dir != 0){
      if(dir == 1){
        simframe_zero<-data.frame(y0=(exp(m_y0$m_y0)), x0=(exp(m_x0$m_x0)), Emax=(exp(m_Emax$m_Emax)), n=(exp(m_n$m_n)))
      } else if (dir == -1){
        simframe_zero<-data.frame(y0=(exp(m_y0$m_y0)), x0=(exp(m_x0$m_x0)), Emax=(inv.logit(m_Emax$m_Emax)), n=(exp(m_n$m_n)))
      }
    } else {
      simframe_zero<-data.frame(y0=(exp(m_y0$m_y0)), x0=(exp(m_x0$m_x0)), n=(exp(m_n$m_n)))
    }
    simframe_zero_conc<- data.frame() 
    for (i in 1:nrow(concframe)) {
      simframe_zero$concentration<-concframe$x[i]
      if(GetFreeConc == T){
        simframe_zero$freeconcentration<-concframe$x[i]*InVivoFree$FracFreeMedia
      }
      simframe_zero_conc<- rbind(simframe_zero_conc, simframe_zero)
    }
    #### Frac change - Population median
    if(dir != 0){
      if(y0Norm == F){simframe_zero_conc$PredFrac<- DR_func(x=simframe_zero_conc$concentration, x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n) 
      } else {simframe_zero_conc$PredFrac<- DR_func(x=simframe_zero_conc$concentration, x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n)}
      if(Scaling == TRUE){simframe_zero_conc$PredFrac<- (simframe_zero_conc$PredFrac)*scale_factor}
    } else {
      if(y0Norm == F){simframe_zero_conc$PredFrac<- DR_func(x=simframe_zero_conc$concentration, x0=simframe_zero_conc$x0, n=simframe_zero_conc$n) 
      } else {simframe_zero_conc$PredFrac<- DR_func(x=simframe_zero_conc$concentration, x0=simframe_zero_conc$x0, n=simframe_zero_conc$n)}
      if(Scaling == TRUE){simframe_zero_conc$PredFrac<- (simframe_zero_conc$PredFrac)*scale_factor}
    }
    if(returnconcframe == T){
      return(simframe_zero_conc)
    }
    ##### Quantiles
    simframe_zero_median_Frac<- (aggregate(PredFrac~concentration, data=simframe_zero_conc, quantile, prob=0.5))
    simframe_zero_975_Frac<- (aggregate(PredFrac~concentration, data=simframe_zero_conc, quantile, prob= prob.upper))
    simframe_zero_25_Frac<- (aggregate(PredFrac~concentration, data=simframe_zero_conc, quantile, prob= prob.lower))
    Simulation_Frame_Pop<- simframe_zero_median_Frac
    if(GetFreeConc == T){
      Simulation_Frame_Pop$freeconcentration <- Simulation_Frame_Pop$concentration*InVivoFree$FracFreeMedia
    }
    Simulation_Frame_Pop$pupper<- simframe_zero_975_Frac$PredFrac
    Simulation_Frame_Pop$plower<- simframe_zero_25_Frac$PredFrac
    names(Simulation_Frame_Pop)[2]<- "p50"
    if(returnparms == FALSE){
      return(Simulation_Frame_Pop)
    } else {return(simframe_zero[,1:4])}
  }
}


## Calculate EC05 from parms

EC05Calc<- function(df = df){
  if(length(df) > 3){
    warning("\n", "Error in EC05Calc : df must be 3 columns")
    stop()
  }
  colnames(df)<- c("x0", "Emax", "n")
  EC05func<- function(x0 = x0, Emax = Emax, n = n){
  x0*(0.05 * Emax / (Emax - 0.05))^(1/n)
  }
  EC05.array<- array()
  for (i in 1:nrow(df)){
    if (is.na(df$Emax[i]) == T){
      temp.EC05 = 1000
      EC05.array[i]
    }
    if (df$Emax[i] > 0.05){
      x0 = df$x0[i]
      Emax = df$Emax[i]
      n = df$n[i]
      temp.EC05 = EC05func(x0, Emax, n)
      EC05.array[i]<- temp.EC05
    } else {
      temp.EC05 = 1000
      EC05.array[i]<- temp.EC05
    }
  }
  return(EC05.array)
}

EC05Calc<- function(df = df){
  if(length(df) > 3){
    warning("\n", "Error in EC05Calc : df must be 3 columns")
    stop()
  }
  colnames(df)<- c("x0", "Emax", "n")
  EC05func<- function(x0 = x0, Emax = Emax, n = n){
    x0*(0.05 * Emax / (Emax - 0.05))^(1/n)
  }
  EC05.array<- array()
  for (i in 1:nrow(df)){
    if (is.na(df$Emax[i]) == T){
      temp.EC05 = 1000
      EC05.array[i]
    }
    if (df$Emax[i] > 0.05){
      x0 = df$x0[i]
      Emax = df$Emax[i]
      n = df$n[i]
      temp.EC05 = EC05func(x0, Emax, n)
      EC05.array[i]<- temp.EC05
    } else {
      temp.EC05 = 1000
      EC05.array[i]<- temp.EC05
    }
  }
  return(EC05.array)
}

EC95Calc<- function(df = df){
  if(length(df) > 2){
    warning("\n", "Error in EC05Calc : df must be 2 columns")
    stop()
  }
  colnames(df)<- c("x0","n")
  EC95func<- function(x0 = x0, n = n){
    x0*(0.95 * 1 / (1 - 0.95))^(1/n)
  }
  EC95.array<- array()
  for (i in 1:nrow(df)){
      x0 = df$x0[i]
      n = df$n[i]
      temp.EC95 = EC95func(x0, n)
      EC95.array[i]<- temp.EC95
    }
  return(EC95.array)
}

EC05Calc.tcz<- function(df = df){
  if(length(df) > 2){
    warning("\n", "Error in EC05Calc : df must be 3 columns")
    stop()
  }
  colnames(df)<- c("x0", "n")
  EC05.tcz.func<- function(x0 = x0, n = n){
    x0*(0.05 * 1 / (1 - 0.05))^(1/n)
  }
  EC05.array<- array()
  for (i in 1:nrow(df)){
      x0 = df$x0[i]
      n = df$n[i]
      temp.EC05 = EC05.tcz.func(x0, n)
      EC05.array[i]<- temp.EC05
  }
  return(EC05.array)
}

EC10Calc.tcz<- function(df = df){
  if(length(df) > 2){
    warning("\n", "Error in EC05Calc : df must be 2 columns")
    stop()
  }
  colnames(df)<- c("x0", "n")
  EC10.tcz.func<- function(x0 = x0, n = n){
    x0*(0.1 * 1 / (1 - 0.1))^(1/n)
  }
  EC10.array<- array()
  for (i in 1:nrow(df)){
    x0 = df$x0[i]
    n = df$n[i]
    temp.EC10 = EC10.tcz.func(x0, n)
    EC10.array[i]<- temp.EC10
  }
  return(EC10.array)
}

EC50Calc.tcz<- function(df = df){
  if(length(df) > 2){
    warning("\n", "Error in EC05Calc : df must be 2 columns")
    stop()
  }
  colnames(df)<- c("x0", "n")
  EC50.tcz.func<- function(x0 = x0, n = n){
    x0*(0.5 * 1 / (1 - 0.5))^(1/n)
  }
  EC50.array<- array()
  for (i in 1:nrow(df)){
    x0 = df$x0[i]
    n = df$n[i]
    temp.EC50 = EC50.tcz.func(x0, n)
    EC50.array[i]<- temp.EC50
  }
  return(EC50.array)
}

## Generate Percent Change at Cmax in vivo vs in vitro

GenCmaxChange <- function(data.list = data.list){
CmaxVivo<-data.frame()
for (j in 1:length(data.list)) {
  for (k in 1:length(levels(data.list[[j]]$InVivo$Model))) {
    CmaxVivo<- rbind(CmaxVivo,with(data.list[[j]], {
      with(subset(InVivo,Model==levels(InVivo$Model)[k]),{
        if (length(xfree)>0) {
          data.frame(j=j,Chemical.name=Chemical.name[1],
                     Model = Model[1],
                     CmaxFree = xfree[length(xfree)],
                     FracChangeCmaxFree=PredFracChange[length(xfree)])
        }
      }
      )
    }
    )
    )
  }
}
CmaxVitro<- data.frame()
for (i in 1:nrow(CmaxVivo)) {
  j<-CmaxVivo$j[i]
  cmaxfree<-CmaxVivo$CmaxFree[i]
  CmaxVitro<- rbind(CmaxVitro,
                    with(data.list[[j]], {
                      with(SimFramePop.Frac, {
                        data.frame(j=j,Chemical.name=chemical.name,
                                   Model="in vitro Population Median",
                                   FracChangeCmaxFree_p50=(approx(x=freeconcentration,
                                                                  y=p50,
                                                                  xout=cmaxfree)),
                                   FracChangeCmaxFree_p2.5=(approx(x=freeconcentration,
                                                                   y=plower,
                                                                   xout=cmaxfree)),
                                   FracChangeCmaxFree_p97.5=(approx(x=freeconcentration,
                                                                    y=pupper,
                                                                    xout=cmaxfree))
                        ) 
                      }
                      )
                    }
                    )
  )
}
names(CmaxVivo)<-paste("InVivo.",names(CmaxVivo),sep="")
names(CmaxVitro)<-paste("InVitro.",names(CmaxVitro),sep="")
Cmax.df <- cbind(rbind(CmaxVivo,CmaxVivo,CmaxVivo,CmaxVivo,CmaxVivo,CmaxVivo),CmaxVitro)
return(Cmax.df)
}

GenECs <- function(data.list = data.list){
  ECVivo<-data.frame()
  for (j in 1:length(data.list)) {
    for (k in 1:length(levels(data.list[[j]]$InVivo$Model))) {
      ECVivo<- rbind(ECVivo,with(data.list[[j]], {
        with(subset(InVivo,Model==levels(InVivo$Model)[k]),{
          if (length(xfree)>0 & sum(PredFracChange) > 0) {
            data.frame(j=j,Chemical.name=Chemical.name[1],
                       Model = Model[1],
                       EC10=(approx(x=PredFracChange, y=xfree, xout=0.1)),
                       EC05=(approx(x=PredFracChange, y=xfree, xout=0.05)), 
                       EC01=(approx(x=PredFracChange, y=xfree, xout=0.01)))
          }
        }
        )
      }
      )
      )
    }
  }
  
  ECVitro<- data.frame()
  for (i in 1:nrow(ECVivo)) {
    j<-ECVivo$j[i]
    ECVitro<- rbind(ECVitro,
                    with(data.list[[j]], {
                      with(SimFramePop.Frac, {
                        data.frame(j=j,Chemical.name=chemical.name,
                                   Model="in vitro Population Median",
                                   EC10_p50=(approx(x=p50, y=freeconcentration, xout=0.1)), 
                                   EC10_p97.5=(approx(x=pupper,y=freeconcentration, xout=0.1)), 
                                   EC10_p2.5=(approx(x=plower, y=freeconcentration, xout=0.1)),
                                   EC05_p50=(approx(x=p50, y=freeconcentration, xout=0.05)), 
                                   EC05_p97.5=(approx(x=pupper,y=freeconcentration, xout=0.05)), 
                                   EC05_p2.5=(approx(x=plower, y=freeconcentration, xout=0.05)), 
                                   EC01_p50=(approx(x=p50, y=freeconcentration, xout=0.01)), 
                                   EC01_p97.5=(approx(x=pupper,y=freeconcentration, xout=0.01)), 
                                   EC01_p2.5=(approx(x=plower, y=freeconcentration, xout=0.01)) )
                      }
                      )
                    }
                    )
    )
  }
  names(ECVivo)<-paste("InVivo.",names(ECVivo),sep="")
  names(ECVitro)<-paste("InVitro.",names(ECVitro),sep="")
  EC.df <- cbind(rbind(ECVivo,ECVivo,ECVivo,ECVivo,ECVivo,ECVivo),ECVitro)
  return(EC.df)
}  

### Stats on IVIV correlation

IVIV.cor <- function(Cmax.df = Cmax.df, EC.df = EC.df){
Cmax.df$Output<-"ECmax"
allEC.df<-Cmax.df[,c("Output","InVivo.Chemical.name",
                     "InVivo.Model",
                     "InVitro.Model",
                     "InVivo.FracChangeCmaxFree",
                     "InVitro.FracChangeCmaxFree_p50.y",
                     "InVitro.FracChangeCmaxFree_p2.5.y",
                     "InVitro.FracChangeCmaxFree_p97.5.y"
)]
names(allEC.df)[c(2,5,6,7,8)]<-c("Chemical.name","InVivo",
                                 "InVitro.p50",
                                 "InVitro.p2.5",
                                 "InVitro.p97.5")

EC10.df <- data.frame(Output=rep("EC10",nrow(EC.df)))
EC10.df <- cbind(EC10.df,EC.df[,c("InVivo.Chemical.name",
                                  "InVivo.Model",
                                  "InVitro.Model",
                                  "InVivo.EC10.y",
                                  "InVitro.EC10_p50.y",
                                  "InVitro.EC10_p2.5.y",
                                  "InVitro.EC10_p97.5.y")])
names(EC10.df)[c(2,5,6,7,8)]<-c("Chemical.name","InVivo",
                                "InVitro.p50",
                                "InVitro.p2.5",
                                "InVitro.p97.5")
EC05.df <- data.frame(Output=rep("EC05",nrow(EC.df)))
EC05.df <- cbind(EC05.df,EC.df[,c("InVivo.Chemical.name",
                                  "InVivo.Model",
                                  "InVitro.Model",
                                  "InVivo.EC05.y",
                                  "InVitro.EC05_p50.y",
                                  "InVitro.EC05_p2.5.y",
                                  "InVitro.EC05_p97.5.y")])
names(EC05.df)[c(2,5,6,7,8)]<-c("Chemical.name","InVivo",
                                "InVitro.p50",
                                "InVitro.p2.5",
                                "InVitro.p97.5")
EC01.df <- data.frame(Output=rep("EC01",nrow(EC.df)))
EC01.df <- cbind(EC01.df,EC.df[,c("InVivo.Chemical.name",
                                  "InVivo.Model",
                                  "InVitro.Model",
                                  "InVivo.EC01.y",
                                  "InVitro.EC01_p50.y",
                                  "InVitro.EC01_p2.5.y",
                                  "InVitro.EC01_p97.5.y")])
names(EC01.df)[c(2,5,6,7,8)]<-c("Chemical.name","InVivo",
                                "InVitro.p50",
                                "InVitro.p2.5",
                                "InVitro.p97.5")

allEC.df<-rbind(EC01.df,EC05.df,EC10.df,allEC.df)

allEC.df$log10err.p50 <- log10(allEC.df$InVitro.p50)-log10(allEC.df$InVivo)
allEC.df$log10err.p97.5 <- log10(allEC.df$InVitro.p97.5)-log10(allEC.df$InVivo)
allEC.df$log10err.p2.5 <- log10(allEC.df$InVitro.p2.5)-log10(allEC.df$InVivo)

allEC.plotting<-subset(allEC.df,InVitro.Model=="in vitro Population Median")
tmpindx<-allEC.plotting$Chemical.name=="Dofetilide" & 
  is.na(allEC.plotting$log10err.p97.5)
allEC.plotting$log10err.p97.5[tmpindx]<-Inf

allEC.summaries<-allEC.plotting[is.finite(allEC.plotting$log10err.p50) &
                                  is.finite(allEC.plotting$log10err.p97.5),]

EC01.sum <- summary(
  EC01.pop.err.lm<-lm(log10err.p50 ~ 1,
                      weights=1/((log10err.p97.5-log10err.p2.5)/(2*qnorm(0.975)))^2,
                      data=subset(allEC.summaries,Output=="EC01" & 
                                    InVitro.Model=="in vitro Population Median"))
)$coefficients

 EC05.sum <- summary(EC05.pop.err.lm<-lm(log10err.p50 ~ 1,
                            weights=1/((log10err.p97.5-log10err.p2.5)/(2*qnorm(0.975)))^2,
                            data=subset(allEC.summaries,Output=="EC05" & 
                                          InVitro.Model=="in vitro Population Median"))
)$coefficients

EC10.sum <- summary(EC10.pop.err.lm<-lm(log10err.p50 ~ 1,
                            weights=1/((log10err.p97.5-log10err.p2.5)/(2*qnorm(0.975)))^2,
                            data=subset(allEC.summaries,Output=="EC10" & 
                                          InVitro.Model=="in vitro Population Median"))
)$coefficients

ECmax.sum <- summary(ECmax.pop.err.lm<-lm(log10err.p50 ~ 1,
                             weights=1/((log10err.p97.5-log10err.p2.5)/(2*qnorm(0.975)))^2,
                             data=subset(allEC.summaries,Output=="ECmax" & 
                                           InVitro.Model=="in vitro Population Median"))
)$coefficients

allEC.df$pcterr <- 100*(allEC.df$InVitro.p50-allEC.df$InVivo)/allEC.df$InVivo

allEC.PopmedianECMax<-subset(allEC.df,Output=="ECmax" & 
                               InVitro.Model=="in vitro Population Median")
tabECmax<-cbind(as.character(allEC.PopmedianECMax$Output),
                as.character(allEC.PopmedianECMax$Chemical.name),
                as.character(allEC.PopmedianECMax$InVivo.Model),
                paste(signif(100*allEC.PopmedianECMax$InVivo,3),"%",sep=""),
                paste(signif(100*allEC.PopmedianECMax$InVitro.p50,3),"% (",
                      signif(100*allEC.PopmedianECMax$InVitro.p2.5,3),"%,",
                      signif(100*allEC.PopmedianECMax$InVitro.p97.5,3),"%)",
                      sep=""))
allEC.PopmedianEC<-subset(allEC.df,Output!="ECmax" & 
                            InVitro.Model=="in vitro Population Median")
tabEC<-cbind(as.character(allEC.PopmedianEC$Output),
             as.character(allEC.PopmedianEC$Chemical.name),
             as.character(allEC.PopmedianEC$InVivo.Model),
             paste(signif(allEC.PopmedianEC$InVivo,3)),
             paste(signif(allEC.PopmedianEC$InVitro.p50,3)," (",
                   signif(allEC.PopmedianEC$InVitro.p97.5,3),",",
                   signif(allEC.PopmedianEC$InVitro.p2.5,3),")",
                   sep=""))
cor.list <- list(EC01.sum, EC05.sum, EC10.sum, ECmax.sum, allEC.df, tabECmax, tabEC)

return(cor.list)

}


### Functions for Prob GEQ10ms

add_qtc_func <-function(simframe_in,qtc0.m=421.5,qtc0.sd=0 
                        # For random individual, use (458.5-421.5)/qnorm(0.95) 
) {
  n <- nrow(simframe_in)
  qtc0<-qtc0.m + qtc0.sd*rnorm(n)
  simframe_out<-simframe_in
  simframe_out$qtc <- simframe_out$PredFrac * qtc0
  simframe_out
}

qtc_geq_func <- function(qtc,qtcstar=10) {
  p_geq <- sum(qtc>=qtcstar)/length(qtc)
  p_geq
}

GenProbGEQ10ms <- function(data.samples.list = data.samples.list, celllines = celllines){
  names(celllines)<-1:length(celllines)
  probgeq10.df <- data.frame()
  for (j in 1:length(data.samples.list)) {
    # Population Median  
    simframe_zero_conc<-data.samples.list[[j]]$simframe_zero_conc
    simframe_out <- add_qtc_func(simframe_zero_conc)
    probgeq10 <- aggregate(qtc~freeconcentration, data=simframe_out, qtc_geq_func)
    probgeq10$safemin<- -0.05
    probgeq10$safemax<- 0.05
    probgeq10$safemin[probgeq10$qtc > 0.05] <- NA
    probgeq10$safemax[probgeq10$qtc > 0.05] <- NA
    probgeq10$InVitroconcentration <- probgeq10$freeconcentration / data.samples.list[[j]]$InVivoFree$FracFreeMedia
    probgeq10$InVivoconcentration <- probgeq10$freeconcentration / data.samples.list[[j]]$InVivoFree$FracFreePlasma
    probgeq10$Chemical.num <- data.samples.list[[j]]$chemical.num 
    probgeq10$Chemical.name <- data.samples.list[[j]]$chemical.name
    probgeq10$InVitroModel <- "Population median"
    probgeq10$CmaxFreeInVivo.x<-max(data.samples.list[[j]]$InVivo$xfree)
    probgeq10$CmaxInVivo.x<-max(data.samples.list[[j]]$InVivo$xfree)/data.samples.list[[j]]$InVivoFree$FracFreePlasma
    probgeq10$CmaxInVivo.y<-1
    probgeq10.df <- rbind(probgeq10.df,probgeq10)
  }
  names(probgeq10.df)[2]<-"PrQTc10"
  probgeq10.df$Chemical.name <- factor(probgeq10.df$Chemical.name)
  probgeq10.df$InVitroModel <- factor(probgeq10.df$InVitroModel,
                                      levels=c("Population median"))
  return(probgeq10.df)
}
