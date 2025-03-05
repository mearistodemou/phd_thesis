#DSEM diagnostics script

#Preliminary steps:
{#0.1 Set working directory
  rm(list = ls()) #clean global environment
  #setwd("~/COTAPP/Data")
  #0.2 Install packages if necessary
  list.of.packages <- c("MplusAutomation", "summarytools", "tidySEM", "rhdf5", "BiocManager", "psych", "gridExtra", "coda", "dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("rhdf5")
  
  #0.3. Load packages
  library(MplusAutomation)
  library(tidySEM)
  library(rhdf5)
  library(psych)
  library(gridExtra)
  library(coda)
  library(dplyr)
}

######################################################    Full DSEM diagnostics   ######################################################
#For diagnostic script to work you must first run the functions from the file: mplus_plots.R; the mplus_plots script was written by
##Thuy Nguyen, Muthen, & Muthen
######################################################    1.1 Diagnostics for block 1   ################################################

#Read model output from Mplus
dsem_full<-MplusAutomation::readModels("C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Mplus\\DSEM\\MPH\\10K\\dsem_mph_10k.out", recursive = FALSE) #plug the name of your output file here
dsem_full2x<-MplusAutomation::readModels("C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Mplus\\DSEM\\MPH\\20K\\dsem_mph_20k.out", recursive = FALSE) #plug the name of your output file here
#Note: this will only work if you used (PLOT: TYPE IS PLOT 2,3)

#STEP 1: Inspect traceplots
pdf("C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Diagnostics\\Output\\MPH\\MPH_traceplot.pdf", width = 6, height = 10)
mplus.traceplot(dsem_full, rows = 4, cols = 2, parameters_only = TRUE)
#Click enter
dev.off() 

#STEP 2: Inspect autocorrelation plots (you have to click saveparameters)
b1_chains <- read.table("C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Mplus\\DSEM\\MPH\\10K\\DSEM_chains_MPH.dat", comment.char="")
b1_chains <-
  b1_chains %>% 
  rename(
    X                        = V1,
    Y                    = V2,
    PHI_MEAN                      = V3,
    LOGV_MEAN                    = V4,
    TREND_MEAN                   = V5,
    LOGRT_MEAN                   = V6,
    PHI_VAR                 = V7,
    LOGVxPHI = V8,
    LOGV_VAR     = V9,
    TRENDwPHI = V10,
    TRENDwLOGV = V11,
    TREND_VAR             = V12,
    LOGRTwPHI = V13,
    LOGRTwLOGV = V14,
    LOGRTwTREND = V15,
    LOGRT_VAR            = V16)

pdf("C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Diagnostics\\Output\\MPH\\MPH_arplot.pdf", width = 6, height = 10)  
autocorr.plot(b1_chains[3:16], 30, auto.layout = TRUE)
#Click enter
dev.off()


#STEP 3: Inspect histogram plots
pdf("C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Diagnostics\\Output\\MPH\\MPH_hist.pdf", width = 16, height = 12) 
par(mfrow = c(5, 3))
par(mar=c(2, 4, 2, 1))
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Mplus\\DSEM\\MPH\\10K\\dsem_mph_10k.gh5',1)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Mplus\\DSEM\\MPH\\10K\\dsem_mph_10k.gh5',2)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Mplus\\DSEM\\MPH\\10K\\dsem_mph_10k.gh5',3)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Mplus\\DSEM\\MPH\\10K\\dsem_mph_10k.gh5',4)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Mplus\\DSEM\\MPH\\10K\\dsem_mph_10k.gh5',5)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Mplus\\DSEM\\MPH\\10K\\dsem_mph_10k.gh5',6)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Mplus\\DSEM\\MPH\\10K\\dsem_mph_10k.gh5',7)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Mplus\\DSEM\\MPH\\10K\\dsem_mph_10k.gh5',8)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Mplus\\DSEM\\MPH\\10K\\dsem_mph_10k.gh5',9)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Mplus\\DSEM\\MPH\\10K\\dsem_mph_10k.gh5',10)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Mplus\\DSEM\\MPH\\10K\\dsem_mph_10k.gh5',11)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Mplus\\DSEM\\MPH\\10K\\dsem_mph_10k.gh5',12)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Mplus\\DSEM\\MPH\\10K\\dsem_mph_10k.gh5',13)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\CatVar\\Mplus\\DSEM\\MPH\\10K\\dsem_mph_10k.gh5',14)
dev.off()

#STEP 4: Calculate relative bias ((initial converged model - double converged model)/initial convergence)*100
#Relative bias function

#DSEM model comparison code
mplus_output1 = dsem_full
mplus_output2 = dsem_full2x

p1 = mplus_output1$parameters$unstandardized
p2 = mplus_output2$parameters$unstandardized

dsem_modelcompare = function(mplus_output1, mplus_output2, value='median'){
  if(!all(names(mplus_output1$savedata) == names(mplus_output2$savedata))){
    exit('Models do not contain the same parameter names, unable to compare...')
  }
  
  p1 = mplus_output1$parameters$unstandardized
  p2 = mplus_output2$parameters$unstandardized
  
  for (i in 1:nrow(p2)){
    j = as.numeric(which(p2$paramHeader[i] == p1$paramHeader & 
                           p2$param[i] == p1$param))
    p2$rel_bias[i] = 100*(p1[i,'est'] - p2[j,'est']) / p1[j,'est']
    p2$std_bias[i] = (p1[j,'est'] - p2[i,'est']) / p1[j,'posterior_sd']
  }
  return(list(p1,p2))
}

p1$link <- row.names(p1)
p2$link <- row.names(p1)

p3<-dplyr::left_join(p1,p2,by="link")
p3_sub <- c("paramHeader.x", "param.x", "est.x", "est.y", "rel_bias")
p3 <- p3[,p3_sub]

#Rename covariances
b1_chains %>% 
  rename(
    LOGVwPHI = V8,
    TRENDwPHI = V10,
    TRENDwLOGV = V11,
    LOGRTwPHI = V13,
    LOGRTwLOGV = V14,
    LOGRTwTREND = V15,
    LOGRT_VAR            = V16)

#LaTeX table of relative bias
knitr::kable(
  p3, #subset to desired columns
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Block 1 relative bias", #change caption
  col.names = c("", "Parameters", "Initial Estimates", "Estimates After 2x Iterations", "Relative Bias"), #rename your columns
  align = "c", #center columns
  row.names = FALSE, #omit row number
) #%>%
#row_spec(row = 0, align = "c")#%>%
#add_footnote(c("Note: Relative bias = (Initial model - Double iteration model)/Initial model*100"))

######################################################    1.2 Diagnostics for block 2   ################################################

#Read model output from Mplus
dsem_full<-MplusAutomation::readModels("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK2\\D2_MEA\\B2_MC\\b2_dsem-full.out", recursive = FALSE) #plug the name of your output file here
dsem_full2x<-MplusAutomation::readModels("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK2\\D2_MEA\\40K_iterations\\b2_dsem-full40k.out", recursive = FALSE) #plug the name of your output file here
#Note: this will only work if you used (PLOT: TYPE IS PLOT 2,3)

#STEP 1: Inspect traceplots
pdf("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\Diagnostics\\B2\\20K\\b2_traceplot.pdf", width = 6, height = 10)
mplus.traceplot(dsem_full, rows = 4, cols = 2, parameters_only = TRUE)
#Click enter
dev.off() 

#STEP 2: Inspect autocorrelation plots (you have to click saveparameters)
b1_chains <- read.table("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK2\\D2_MEA\\B2_MC\\DSEM_chains_b2.dat", comment.char="")
b1_chains <-
  b1_chains %>% 
  rename(
    X                        = V1,
    Y                    = V2,
    PHI_MEAN                      = V3,
    LOGV_MEAN                    = V4,
    TREND_MEAN                   = V5,
    LOGRT_MEAN                   = V6,
    PHI_VAR                 = V7,
    LOGV_VAR     = V8,
    TREND_VAR             = V9,
    LOGRT_VAR            = V10)

pdf("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\Diagnostics\\B2\\20K\\b2_arplot.pdf", width = 6, height = 10)  
autocorr.plot(b1_chains[3:10], 30, auto.layout = TRUE)
#Click enter
dev.off()


#STEP 3: Inspect histogram plots
pdf("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\Diagnostics\\B2\\20K\\b2_histplot.pdf", width = 16, height = 12) 
par(mfrow = c(3, 3))
par(mar=c(2, 4, 2, 1))
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK2\\D2_MEA\\B2_MC\\b2_dsem-full.gh5',1)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK2\\D2_MEA\\B2_MC\\b2_dsem-full.gh5',2)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK2\\D2_MEA\\B2_MC\\b2_dsem-full.gh5',3)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK2\\D2_MEA\\B2_MC\\b2_dsem-full.gh5',4)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK2\\D2_MEA\\B2_MC\\b2_dsem-full.gh5',5)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK2\\D2_MEA\\B2_MC\\b2_dsem-full.gh5',6)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK2\\D2_MEA\\B2_MC\\b2_dsem-full.gh5',7)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK2\\D2_MEA\\B2_MC\\b2_dsem-full.gh5',8)
dev.off()

#STEP 4: Calculate relative bias ((initial converged model - double converged model)/initial convergence)*100

#DSEM model comparison code
mplus_output1 = dsem_full
mplus_output2 = dsem_full2x

p1 = mplus_output1$parameters$unstandardized
p2 = mplus_output2$parameters$unstandardized

dsem_modelcompare = function(mplus_output1, mplus_output2, value='median'){
  if(!all(names(mplus_output1$savedata) == names(mplus_output2$savedata))){
    exit('Models do not contain the same parameter names, unable to compare...')
  }
  
  p1 = mplus_output1$parameters$unstandardized
  p2 = mplus_output2$parameters$unstandardized
  
  for (i in 1:nrow(p2)){
    j = as.numeric(which(p2$paramHeader[i] == p1$paramHeader & 
                           p2$param[i] == p1$param))
    p2$rel_bias[i] = 100*(p1[i,'est'] - p2[j,'est']) / p1[j,'est']
    p2$std_bias[i] = (p1[j,'est'] - p2[i,'est']) / p1[j,'posterior_sd']
  }
  return(list(p1,p2))
}

p1$link <- row.names(p1)
p2$link <- row.names(p1)

p3<-dplyr::left_join(p1,p2,by="link")
p3_sub <- c("paramHeader.x", "param.x", "est.x", "est.y", "rel_bias")
p3 <- p3[,p3_sub]

#Rename covariances


#LaTeX table of relative bias
knitr::kable(
  p3, #subset to desired columns
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Block 2 relative bias", #change caption
  col.names = c("", "Parameters", "Initial Estimates", "Estimates After 2x Iterations", "Relative Bias"), #rename your columns
  align = "c", #center columns
  row.names = FALSE, #omit row number
) #%>%
#row_spec(row = 0, align = "c")#%>%
#add_footnote(c("Note: Relative bias = (Initial model - Double iteration model)/Initial model*100"))

######################################################    1.2 Diagnostics for block 3   ################################################

#Read model output from Mplus
dsem_full<-MplusAutomation::readModels("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK3\\D2_MEA\\B3_MC\\b3_dsem-full.out", recursive = FALSE) #plug the name of your output file here
dsem_full2x<-MplusAutomation::readModels("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK3\\D2_MEA\\40K_iterations\\dsem-full_40k_b3.out", recursive = FALSE) #plug the name of your output file here
#Note: this will only work if you used (PLOT: TYPE IS PLOT 2,3)

#STEP 1: Inspect traceplots
pdf("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\Diagnostics\\B3\\20K\\b3_traceplot.pdf", width = 6, height = 10)
mplus.traceplot(dsem_full, rows = 4, cols = 2, parameters_only = TRUE)
#Click enter
dev.off() 

#STEP 2: Inspect autocorrelation plots (you have to click saveparameters)
b1_chains <- read.table("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK3\\D2_MEA\\B3_MC\\DSEM_chains_b3.dat", comment.char="")
b1_chains <-
  b1_chains %>% 
  rename(
    X                        = V1,
    Y                    = V2,
    PHI_MEAN                      = V3,
    LOGV_MEAN                    = V4,
    TREND_MEAN                   = V5,
    LOGRT_MEAN                   = V6,
    PHI_VAR                 = V7,
    LOGV_VAR     = V8,
    TREND_VAR             = V9,
    LOGRT_VAR            = V10)

pdf("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\Diagnostics\\B3\\20K\\b3_arplot.pdf", width = 6, height = 10)  
autocorr.plot(b1_chains[3:10], 30, auto.layout = TRUE)
#Click enter
dev.off()


#STEP 3: Inspect histogram plots
pdf("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\Diagnostics\\B3\\20K\\b3_histplot.pdf", width = 16, height = 12) 
par(mfrow = c(3, 3))
par(mar=c(2, 4, 2, 1))
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK3\\D2_MEA\\B3_MC\\b3_dsem-full.gh5',1)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK3\\D2_MEA\\B3_MC\\b3_dsem-full.gh5',2)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK3\\D2_MEA\\B3_MC\\b3_dsem-full.gh5',3)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK3\\D2_MEA\\B3_MC\\b3_dsem-full.gh5',4)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK3\\D2_MEA\\B3_MC\\b3_dsem-full.gh5',5)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK3\\D2_MEA\\B3_MC\\b3_dsem-full.gh5',6)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK3\\D2_MEA\\B3_MC\\b3_dsem-full.gh5',7)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK3\\D2_MEA\\B3_MC\\b3_dsem-full.gh5',8)
dev.off()

#STEP 4: Calculate relative bias ((initial converged model - double converged model)/initial convergence)*100

#DSEM model comparison code
mplus_output1 = dsem_full
mplus_output2 = dsem_full2x

p1 = mplus_output1$parameters$unstandardized
p2 = mplus_output2$parameters$unstandardized

dsem_modelcompare = function(mplus_output1, mplus_output2, value='median'){
  if(!all(names(mplus_output1$savedata) == names(mplus_output2$savedata))){
    exit('Models do not contain the same parameter names, unable to compare...')
  }
  
  p1 = mplus_output1$parameters$unstandardized
  p2 = mplus_output2$parameters$unstandardized
  
  for (i in 1:nrow(p2)){
    j = as.numeric(which(p2$paramHeader[i] == p1$paramHeader & 
                           p2$param[i] == p1$param))
    p2$rel_bias[i] = 100*(p1[i,'est'] - p2[j,'est']) / p1[j,'est']
    p2$std_bias[i] = (p1[j,'est'] - p2[i,'est']) / p1[j,'posterior_sd']
  }
  return(list(p1,p2))
}

p1$link <- row.names(p1)
p2$link <- row.names(p1)

p3<-dplyr::left_join(p1,p2,by="link")
p3_sub <- c("paramHeader.x", "param.x", "est.x", "est.y", "rel_bias")
p3 <- p3[,p3_sub]

#Rename covariances
b1_chains <-
  b1_chains %>% 
  rename(
    LOGV_MEAN                    = V4,
    TREND_MEAN                   = V5,
    LOGRT_MEAN                   = V6,
    PHI_VAR                 = V7,
    LOGV_VAR     = V8,
    TREND_VAR             = V9,
    LOGRT_VAR            = V10)

#LaTeX table of relative bias
knitr::kable(
  p3, #subset to desired columns
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Block 3 relative bias", #change caption
  col.names = c("", "Parameters", "Initial Estimates", "Estimates After 2x Iterations", "Relative Bias"), #rename your columns
  align = "c", #center columns
  row.names = FALSE, #omit row number
) #%>%
#row_spec(row = 0, align = "c")#%>%
#add_footnote(c("Note: Relative bias = (Initial model - Double iteration model)/Initial model*100"))

######################################################    1.5 Diagnostics for block 4   ################################################

#Read model output from Mplus
dsem_full<-MplusAutomation::readModels("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK4\\D2_MEA\\B4_MC\\b4_dsem-full.out", recursive = FALSE) #plug the name of your output file here
dsem_full2x<-MplusAutomation::readModels("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK4\\D2_MEA\\40K_iterations\\dsem-full_40k_b4.out", recursive = FALSE) #plug the name of your output file here
#Note: this will only work if you used (PLOT: TYPE IS PLOT 2,3)

#STEP 1: Inspect traceplots
pdf("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\Diagnostics\\B4\\20K\\b4_traceplot.pdf", width = 6, height = 10)
mplus.traceplot(dsem_full, rows = 4, cols = 2, parameters_only = TRUE)
#Click enter
dev.off() 

#STEP 2: Inspect autocorrelation plots (you have to click saveparameters)
b1_chains <- read.table("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK4\\D2_MEA\\B4_MC\\DSEM_chains_b4.dat", comment.char="")
b1_chains <-
  b1_chains %>% 
  rename(
    X                        = V1,
    Y                    = V2,
    PHI_MEAN                      = V3,
    LOGV_MEAN                    = V4,
    TREND_MEAN                   = V5,
    LOGRT_MEAN                   = V6,
    PHI_VAR                 = V7,
    LOGV_VAR     = V8,
    TREND_VAR             = V9,
    LOGRT_VAR            = V10)

pdf("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\Diagnostics\\B4\\20K\\b4_arplot.pdf", width = 6, height = 10)  
autocorr.plot(b1_chains[3:10], 30, auto.layout = TRUE)
#Click enter
dev.off()


#STEP 3: Inspect histogram plots
pdf("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\Diagnostics\\B4\\20K\\b4_histplot.pdf", width = 16, height = 12) 
par(mfrow = c(3, 3))
par(mar=c(2, 4, 2, 1))
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK4\\D2_MEA\\B4_MC\\b4_dsem-full.gh5',1)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK4\\D2_MEA\\B4_MC\\b4_dsem-full.gh5',2)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK4\\D2_MEA\\B4_MC\\b4_dsem-full.gh5',3)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK4\\D2_MEA\\B4_MC\\b4_dsem-full.gh5',4)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK4\\D2_MEA\\B4_MC\\b4_dsem-full.gh5',5)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK4\\D2_MEA\\B4_MC\\b4_dsem-full.gh5',6)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK4\\D2_MEA\\B4_MC\\b4_dsem-full.gh5',7)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK4\\D2_MEA\\B4_MC\\b4_dsem-full.gh5',8)
dev.off()

#STEP 4: Calculate relative bias ((initial converged model - double converged model)/initial convergence)*100

#DSEM model comparison code
mplus_output1 = dsem_full
mplus_output2 = dsem_full2x

p1 = mplus_output1$parameters$unstandardized
p2 = mplus_output2$parameters$unstandardized

dsem_modelcompare = function(mplus_output1, mplus_output2, value='median'){
  if(!all(names(mplus_output1$savedata) == names(mplus_output2$savedata))){
    exit('Models do not contain the same parameter names, unable to compare...')
  }
  
  p1 = mplus_output1$parameters$unstandardized
  p2 = mplus_output2$parameters$unstandardized
  
  for (i in 1:nrow(p2)){
    j = as.numeric(which(p2$paramHeader[i] == p1$paramHeader & 
                           p2$param[i] == p1$param))
    p2$rel_bias[i] = 100*(p1[i,'est'] - p2[j,'est']) / p1[j,'est']
    p2$std_bias[i] = (p1[j,'est'] - p2[i,'est']) / p1[j,'posterior_sd']
  }
  return(list(p1,p2))
}

p1$link <- row.names(p1)
p2$link <- row.names(p1)

p3<-dplyr::left_join(p1,p2,by="link")
p3_sub <- c("paramHeader.x", "param.x", "est.x", "est.y", "rel_bias")
p3 <- p3[,p3_sub]

#LaTeX table of relative bias
knitr::kable(
  p3, #subset to desired columns
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Block 4 relative bias", #change caption
  col.names = c("", "Parameters", "Initial Estimates", "Estimates After 2x Iterations", "Relative Bias"), #rename your columns
  align = "c", #center columns
  row.names = FALSE, #omit row number
) #%>%
#row_spec(row = 0, align = "c")#%>%
#add_footnote(c("Note: Relative bias = (Initial model - Double iteration model)/Initial model*100"))

######################################################    1.5 Diagnostics for block 7A   ################################################

#Read model output from Mplus
dsem_full<-MplusAutomation::readModels("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7A\\D2_MEA\\B7A_MC\\b7a_dsem-full.out", recursive = FALSE) #plug the name of your output file here
dsem_full2x<-MplusAutomation::readModels("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7A\\D2_MEA\\40K_iterations\\dsem-full_40k_b7a.out", recursive = FALSE) #plug the name of your output file here
#Note: this will only work if you used (PLOT: TYPE IS PLOT 2,3)

#STEP 1: Inspect traceplots
pdf("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\Diagnostics\\B7A\\20K\\b7A_traceplot.pdf", width = 6, height = 10)
mplus.traceplot(dsem_full, rows = 4, cols = 2, parameters_only = TRUE)
#Click enter
dev.off() 

#STEP 2: Inspect autocorrelation plots (you have to click saveparameters)
b1_chains <- read.table("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7A\\D2_MEA\\B7A_MC\\DSEM_chains_b71.dat", comment.char="")
b1_chains <-
  b1_chains %>% 
  rename(
    X                        = V1,
    Y                    = V2,
    PHI_MEAN                      = V3,
    LOGV_MEAN                    = V4,
    TREND_MEAN                   = V5,
    LOGRT_MEAN                   = V6,
    PHI_VAR                 = V7,
    LOGV_VAR     = V8,
    TREND_VAR             = V9,
    LOGRT_VAR            = V10)

pdf("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\Diagnostics\\B7A\\20K\\b7A_arplot.pdf", width = 6, height = 10)  
autocorr.plot(b1_chains[3:10], 30, auto.layout = TRUE)
#Click enter
dev.off()


#STEP 3: Inspect histogram plots
pdf("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\Diagnostics\\B7A\\20K\\b7A_histplot.pdf", width = 16, height = 12) 
par(mfrow = c(3, 3))
par(mar=c(2, 4, 2, 1))
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7A\\D2_MEA\\B7A_MC\\b7a_dsem-full.gh5',1)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7A\\D2_MEA\\B7A_MC\\b7a_dsem-full.gh5',2)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7A\\D2_MEA\\B7A_MC\\b7a_dsem-full.gh5',3)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7A\\D2_MEA\\B7A_MC\\b7a_dsem-full.gh5',4)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7A\\D2_MEA\\B7A_MC\\b7a_dsem-full.gh5',5)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7A\\D2_MEA\\B7A_MC\\b7a_dsem-full.gh5',6)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7A\\D2_MEA\\B7A_MC\\b7a_dsem-full.gh5',7)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7A\\D2_MEA\\B7A_MC\\b7a_dsem-full.gh5',8)
dev.off()

#STEP 4: Calculate relative bias ((initial converged model - double converged model)/initial convergence)*100

#DSEM model comparison code
mplus_output1 = dsem_full
mplus_output2 = dsem_full2x

p1 = mplus_output1$parameters$unstandardized
p2 = mplus_output2$parameters$unstandardized

dsem_modelcompare = function(mplus_output1, mplus_output2, value='median'){
  if(!all(names(mplus_output1$savedata) == names(mplus_output2$savedata))){
    exit('Models do not contain the same parameter names, unable to compare...')
  }
  
  p1 = mplus_output1$parameters$unstandardized
  p2 = mplus_output2$parameters$unstandardized
  
  for (i in 1:nrow(p2)){
    j = as.numeric(which(p2$paramHeader[i] == p1$paramHeader & 
                           p2$param[i] == p1$param))
    p2$rel_bias[i] = 100*(p1[i,'est'] - p2[j,'est']) / p1[j,'est']
    p2$std_bias[i] = (p1[j,'est'] - p2[i,'est']) / p1[j,'posterior_sd']
  }
  return(list(p1,p2))
}

p1$link <- row.names(p1)
p2$link <- row.names(p1)

p3<-dplyr::left_join(p1,p2,by="link")
p3_sub <- c("paramHeader.x", "param.x", "est.x", "est.y", "rel_bias")
p3 <- p3[,p3_sub]

#LaTeX table of relative bias
knitr::kable(
  p3, #subset to desired columns
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Block 7A relative bias", #change caption
  col.names = c("", "Parameters", "Initial Estimates", "Estimates After 2x Iterations", "Relative Bias"), #rename your columns
  align = "c", #center columns
  row.names = FALSE, #omit row number
) #%>%
#row_spec(row = 0, align = "c")#%>%
#add_footnote(c("Note: Relative bias = (Initial model - Double iteration model)/Initial model*100"))

######################################################    1.5 Diagnostics for block 7B   ################################################

#Read model output from Mplus
dsem_full<-MplusAutomation::readModels("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7B\\D2_MEA\\40K_iterations\\dsem-full_40k_b7b.out", recursive = FALSE) #plug the name of your output file here
dsem_full2x<-MplusAutomation::readModels("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7B\\D2_MEA\\40K_iterations\\dsem-full_80k_b7b.out", recursive = FALSE) #plug the name of your output file here
#Note: this will only work if you used (PLOT: TYPE IS PLOT 2,3)

#STEP 1: Inspect traceplots
pdf("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\Diagnostics\\B7B\\40K\\b7B_traceplot.pdf", width = 6, height = 10)
mplus.traceplot(dsem_full, rows = 4, cols = 2, parameters_only = TRUE)
#Click enter
dev.off() 

#STEP 2: Inspect autocorrelation plots (you have to click saveparameters)
b1_chains <- read.table("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7B\\D2_MEA\\40K_iterations\\DSEM_chains_b7240.dat", comment.char="")
b1_chains <-
  b1_chains %>% 
  rename(
    X                        = V1,
    Y                    = V2,
    PHI_MEAN                      = V3,
    LOGV_MEAN                    = V4,
    TREND_MEAN                   = V5,
    LOGRT_MEAN                   = V6,
    PHI_VAR                 = V7,
    LOGV_VAR     = V8,
    TREND_VAR             = V9,
    LOGRT_VAR            = V10)

pdf("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\Diagnostics\\B7B\\40K\\b7B_arplot.pdf", width = 6, height = 10)  
autocorr.plot(b1_chains[3:10], 30, auto.layout = TRUE)
#Click enter
dev.off()


#STEP 3: Inspect histogram plots
pdf("C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\Diagnostics\\B7B\\40K\\b7B_histplot.pdf", width = 16, height = 12) 
par(mfrow = c(3, 3))
par(mar=c(2, 4, 2, 1))
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7B\\D2_MEA\\40K_iterations\\dsem-full_40k_b7b.gh5',1)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7B\\D2_MEA\\40K_iterations\\dsem-full_40k_b7b.gh5',2)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7B\\D2_MEA\\40K_iterations\\dsem-full_40k_b7b.gh5',3)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7B\\D2_MEA\\40K_iterations\\dsem-full_40k_b7b.gh5',4)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7B\\D2_MEA\\40K_iterations\\dsem-full_40k_b7b.gh5',5)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7B\\D2_MEA\\40K_iterations\\dsem-full_40k_b7b.gh5',6)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7B\\D2_MEA\\40K_iterations\\dsem-full_40k_b7b.gh5',7)
mplus.plot.bayesian.distribution('C:\\Users\\micha\\OneDrive\\Documents\\COTAPP\\BLOCK7B\\D2_MEA\\40K_iterations\\dsem-full_40k_b7b.gh5',8)
dev.off()

#STEP 4: Calculate relative bias ((initial converged model - double converged model)/initial convergence)*100

#DSEM model comparison code
mplus_output1 = dsem_full
mplus_output2 = dsem_full2x

p1 = mplus_output1$parameters$unstandardized
p2 = mplus_output2$parameters$unstandardized

dsem_modelcompare = function(mplus_output1, mplus_output2, value='median'){
  if(!all(names(mplus_output1$savedata) == names(mplus_output2$savedata))){
    exit('Models do not contain the same parameter names, unable to compare...')
  }
  
  p1 = mplus_output1$parameters$unstandardized
  p2 = mplus_output2$parameters$unstandardized
  
  for (i in 1:nrow(p2)){
    j = as.numeric(which(p2$paramHeader[i] == p1$paramHeader & 
                           p2$param[i] == p1$param))
    p2$rel_bias[i] = 100*(p1[i,'est'] - p2[j,'est']) / p1[j,'est']
    p2$std_bias[i] = (p1[j,'est'] - p2[i,'est']) / p1[j,'posterior_sd']
  }
  return(list(p1,p2))
}

p1$link <- row.names(p1)
p2$link <- row.names(p1)

p3<-dplyr::left_join(p1,p2,by="link")
p3_sub <- c("paramHeader.x", "param.x", "est.x", "est.y", "rel_bias")
p3 <- p3[,p3_sub]

#LaTeX table of relative bias
knitr::kable(
  p3, #subset to desired columns
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Block 7B relative bias", #change caption
  col.names = c("", "Parameters", "Initial Estimates", "Estimates After 2x Iterations", "Relative Bias"), #rename your columns
  align = "c", #center columns
  row.names = FALSE, #omit row number
) #%>%
#row_spec(row = 0, align = "c")#%>%
#add_footnote(c("Note: Relative bias = (Initial model - Double iteration model)/Initial model*100"))

########################################################################################################################################
