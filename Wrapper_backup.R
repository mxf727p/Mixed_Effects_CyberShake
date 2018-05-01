#########################
### Residuals Analysis For Predicted Ground Motion Equations
### 03/29/2018   Started by Xiaofeng Meng, only ASK14 and a simple model 
###              by Albert Kottke (i.e, KOT) is implemented now.

# Clear the console
cat("\014")  
# Clear the workspace
rm(list=ls())

require(ggplot2)
require(cowplot)
require(nlme)
require(lme4)
require(minpack.lm)
require(scales)
library(RColorBrewer)



########################
### Set up Dataset, Model, Working directories, Plotting flags, Debug flag etc.
Dataname <- "NGA_W2_SC"     # Options: "NGA_W2", "CyberShake", "NGA_W2_SC", "CyberShake_SC_regression"
GM_model <- "VIL"        # Options: "ASK14", "KOT", "VIL"
path_main <- setwd("/Users/xmeng/Documents/MyWork/CyberShake")
regression_flag <- 1
if (regression_flag) {
  regression_model <- 4   # 1: nlme; 2: nls; 3: nlsLM; 4: nlmer
}
residual_flag <- 0
if (residual_flag == 0) {
  read_residual <- 0
}
L2L_flag <- 0
S2S_flag <- 0
P2P_flag <- 1
output_flag <- 0
plot_residual <- 1
plot_response <- 0
debug_flag <- 0
### End set up
########################


#########################
### Create working directories
# Path to dataset 
path_data <- file.path(path_main,"Data") 
# Path to model
path_model <- file.path(path_main,"Model") 
dir.create(file.path(path_model,GM_model),showWarnings = FALSE)
dir.create(file.path(path_model,GM_model,Dataname),showWarnings = FALSE)
# Output files directory
dir.create(file.path(path_model,GM_model,Dataname,"Outputs"),showWarnings = FALSE)
path_outputs <- paste0(path_model,"/",GM_model,"/",Dataname,"/Outputs")
# Figures directory
dir.create(file.path(path_model,GM_model,Dataname,"Figures"),showWarnings = FALSE)
path_figures <- paste0(path_model,"/",GM_model,"/",Dataname,"/Figures")
# Debug directory
dir.create(file.path(path_model,GM_model,Dataname,"Debug"),showWarnings = FALSE)
path_debug <- paste0(path_model,"/",GM_model,"/",Dataname,"/Debug")
### End create working directories
#########################


#########################
### Read dataset
source(file.path(path_data,"Read_data.R"))
Datafile <- paste0(path_data,"/",Dataname,".csv")
Dataset <- Read_data(Dataname, Datafile)
# Global variables
# [eqid]: eq. id; [mag]: eq. magnitude; [rrup]: closest distance to rupture plane (km);
# [dip]: dip angle; [rake]: rake angle; [vs30]: Vs30 (m/s); rx: Horizontal distance (km);
# [w]: rupture width (km); [ztor]: depth to top of fault rupture (km);
# [rjb]: Joyner Boore Dist (km); [z1_true]: depth to Vs=1km/s; #[crjb]: CRjb;
# [sof]: style of faulting; [gmx]: site claasification; [hwflag]: hanging wall flag
# [ry0]: ???; [region]: region flag
### End read dataset
#########################


#########################
### Regression analysis
if (regression_flag) {
  #regression_list <- data.frame(c(0.010,0.020,0.030,0.050,0.075,0.100,0.150,0.200,0.250,0.300,
  #                                0.400,0.500,0.750,1.000,1.500,2.000,3.000,4.000,5.000,6.000,
  #                                7.500,10.000))
  regression_list <- data.frame(c(3.000))
  if (GM_model == "KOT") {
    lnSa <- function(mag, dist, c1, c2, c4, c6, c7, c10, c11) {
    return(c1 + c2 * mag + (c6 + c7 * mag) * log(dist + exp(c4)) +
               c10 * (mag - 6) ^ 2 + c11 * dist)
    }

    colnames(regression_list) <- c("Period")
    selected <- cbind(Dataset$eqid, Dataset$mag, Dataset$rrup, Dataset$stid, Dataset$pathid)
    colnames(selected) <- c("eqid","mag","rrup","stid","pathid")
    coeff <- data.frame(matrix(NA, nrow = 1, ncol = 7))
    coeff <- coeff[-1,]
    colnames(coeff) <- c("c1","c2","c4","c6","c7","c10","c11")
    for (i in 1:nrow(regression_list)) {
      Period <- regression_list[i,1]
      obs_regression <- log(Dataset[colnames(Dataset)==format(round(Period,3),nsmall=3)])
      colnames(obs_regression) <- c("obs")
      selected <- cbind(selected, obs_regression)
      if (regression_model == 1) {
        fitModel <- nlme(obs ~ lnSa(mag, rrup, c1, c2, c4, c6, c7, c10, c11),
                         data = selected[selected$obs!="NaN",],
                         fixed = c1 + c2 + c4 + c6 + c7 + c10 + c11 ~ 1,
                         random = c1 ~ 1|eqid,
                         start=c(c1=0., c2=0.5, c4=2.8, c6=-3.0, c7=0.2, c10=-0.1,c11=-0.02))
        coeff <- rbind(coeff,t(fitModel$coefficients$fixed))
      } else if (regression_model == 2) {
        fitModel <- nls(obs ~ lnSa(mag, rrup, c1, c2, c4, c6, c7, c10, c11),
                        data = selected[selected$obs!="NaN",],
                        start=c(c1=0., c2=0.5, c4=2.8, c6=-3.0, c7=0.2, c10=-0.1,c11=-0.02))
        coeff <- rbind(coeff,t(coef(fitModel)))
      } else if (regression_model == 3) {
        fitModel <- nlsLM(obs ~ lnSa(mag, rrup, c1, c2, c4, c6, c7, c10, c11),
                          data = selected[selected$obs!="NaN",])
        coeff <- rbind(coeff,t(coef(fitModel)))
      }
      selected <- subset(selected, select = -c(obs))
    }
    coeff <- cbind(regression_list,coeff)
    write.csv(coeff,file=file.path(path_model,GM_model,Dataname,"Coeff.csv"))
  } else if (GM_model == "VIL") {
    lnSa <- function(mag, dist, vs30, b1, b2, b3, b4, b5, b6) {
            return(b1 + b2 * mag + b3 * mag ^ 2 + (b4 + b5 * mag) * 
               log(sqrt(4 ^ 2 + dist ^ 2)) + b6 * log(vs30/500))
    }
    colnames(regression_list) <- c("Period")
    coeff <- data.frame(matrix(NA, nrow = 1, ncol = 6))
    coeff <- coeff[-1,]
    colnames(coeff) <- c("b1","b2","b3","b4","b5","b6")
    for (i in 1:nrow(regression_list)) {
      Period <- regression_list[i,1]
      obs_regression <- log(Dataset[colnames(Dataset)==format(round(Period,3),nsmall=3)])
      colnames(obs_regression) <- c("obs")
      Dataset <- cbind(Dataset, obs_regression)
      if (regression_model == 1) {
        fitModel <- nlme(obs ~ lnSa(mag, rrup, vs30, b1, b2, b3, b4, b5, b6),
                         data = Dataset[Dataset$obs!="NaN",],
                         fixed = b1 + b2 + b3 + b4 + b5 + b6 ~ 1,
                         random = b1 ~ 1|eqid,
                         start=c(b1=-6, b2=3, b3=0.0, b4=-2.6, b5=0.3, b6=-1.29))
        coeff <- rbind(coeff,t(fitModel$coefficients$fixed))
      } else if (regression_model == 2) {
        fitModel <- nls(obs ~ lnSa(mag, rrup, vs30, b1, b2, b3, b4, b5, b6),
                       data = Dataset[Dataset$obs!="NaN",],
                       start=list(b1=-6, b2=3, b3=0.0, b4=-2.6, b5=0.3, b6=-1.29))
        coeff <- rbind(coeff,t(coef(fitModel)))
      } else if (regression_model == 3) {
        fitModel <- nlsLM(obs ~ lnSa(mag, rrup, vs30, b1, b2, b3, b4, b5, b6),
                         data = Dataset[Dataset$obs!="NaN",])
        coeff <- rbind(coeff,t(coef(fitModel))) 
      } else if (regression_model == 4) {
        lnSa.deriv = deriv(~ b1 + b2 * mag + b3 * mag ^ 2 + (b4 + b5 * mag) * log(sqrt(4^2 + rrup^2)) + b6 * log(vs30/500), 
                           namevec=c('b1','b2','b3','b4','b5','b6'), 
                           function.arg = c('mag','rrup','vs30','b1','b2','b3','b4','b5','b6'))
        fitModel <- nlmer(obs ~ lnSa.deriv(mag, rrup, vs30, b1, b2, b3, b4, b5, b6) ~ (b1|eqid) + (b1|stid) + (b1|pathid),
                        data=Dataset[Dataset$obs!="NaN",], start=c(b1=-6, b2=3, b3=0.0, b4=-2.6, b5=0.3, b6=-1.29))
        coeff <- rbind(coeff,t(coef(summary(fitModel))[,"Estimate"]))
      }
      Dataset <- subset(Dataset, select = -c(obs))
    }
    coeff <- cbind(regression_list,coeff)
    write.csv(coeff,file=file.path(path_model,GM_model,Dataname,"Coeff.csv"))
  }
  # End select model
  #remove("selected","coeff","Period","obs_regression")
}
### End regression analysis
#########################


#########################
### Read equations and coefficients for the model
source(file.path(path_model,GM_model,"Equations.R"))
coeff <- Read_coeff(file.path(path_model,GM_model,Dataname,"Coeff.csv"))
### End read equations and coefficients
#########################


#########################
### Predict ground motions and compute residuals
PGM <- data.frame(matrix(NA, nrow = nrow(Dataset), ncol = nrow(coeff)))
Obs <- data.frame(matrix(NA, nrow = nrow(Dataset), ncol = 1))
for (i in 1:nrow(coeff)) {
  # Read coefficients for GM_model
  for (j in 1:ncol(coeff)) {
    assign(colnames(coeff)[j], coeff[i,j])
  }
  # End read coefficients
  if (GM_model == "ASK14") {
    lnSa1180 <- apply(Dataset[,c('mag','w','dip','ztor','sof','hwflag','rrup','rjb','rx','ry0','vs30','z1_true','region')],1,
                function(Dataset)  
                lnSa(Dataset['mag'],Dataset['w'],Dataset['dip'],Dataset['ztor'],
                Dataset['sof'],Dataset['hwflag'],Dataset['rrup'],Dataset['rjb'],
                Dataset['rx'],Dataset['ry0'],1180,0,-1,Dataset['region']))
    Sa1180 <- exp(lnSa1180)
    Dataset <- cbind(Dataset,Sa1180) 
    PGM[,i] <- apply(Dataset[,c('mag','w','dip','ztor','sof','hwflag','rrup','rjb','rx','ry0','vs30','Sa1180','z1_true','region')],1,
               function(Dataset)  
               lnSa(Dataset['mag'],Dataset['w'],Dataset['dip'],Dataset['ztor'],
               Dataset['sof'],Dataset['hwflag'],Dataset['rrup'],Dataset['rjb'],
               Dataset['rx'],Dataset['ry0'],Dataset['vs30'],Dataset['Sa1180'],
               Dataset['z1_true'],Dataset['region']))
  } else if (GM_model == "KOT") {
    Sa1180 <- matrix(0, nrow = nrow(Dataset), ncol = 1)
    Dataset <- cbind(Dataset,Sa1180) 
    PGM[,i] <- apply(Dataset[,c('mag','rrup')],1,function(Dataset)  lnSa(Dataset['mag'],Dataset['rrup']))
  } else if (GM_model == "VIL") {
    Sa1180 <- matrix(0, nrow = nrow(Dataset), ncol = 1)
    Dataset <- cbind(Dataset,Sa1180) 
    PGM[,i] <- apply(Dataset[,c('mag','rrup','vs30')],1,function(Dataset)  lnSa(Dataset['mag'],Dataset['rrup'],Dataset['vs30']))
  }
  colnames(PGM)[i] <- coeff$Period[i]
  
  # Read observed ground motions from dataset
  # Obs[,i] <- Dataset[colnames(Dataset)==format(round(coeff$Period[i],3),nsmall=3)]
  
  # Compute residuals
  Period <- coeff$Period[i]
  Obs <- log(Dataset[colnames(Dataset)==format(round(coeff$Period[i],3),nsmall=3)]) 
  colnames(Obs) <- c("obs")
  Resi <- Obs - PGM[,i]
  colnames(Resi) <- c("Resi")
  Dataset <- cbind(Dataset,Resi,Obs)
  Dataset <- Dataset[Dataset$Resi!="NaN",]
  
  
  if (L2L_flag) {
    model <- lme(fixed = Resi ~ 1, data = Dataset2, method="REML", random = ~ 1 | EQID)
    TotalR = residuals(model, level = 0)
    epsR = residuals(model, level = 1)
    etaR = data.frame(model$coefficients$random$EQID)
    #etaR = coef(model)
    colnames(etaR) <- c("L2L")
    p1 <- ggplot(etaR, aes(x=L2L)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
                 stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(etaR$L2L))) * nrow(etaR) * 0.1, color = "red", size = 1) +
                 annotate("text", x=-0.8, y=5, label = paste0("phi[L2L] ==",formatC(sqrt(var(etaR$L2L)),digits=3,format="f",flag="0")), color="red", parse=TRUE) +
                 annotate("text", x=0, y=20, label = "lme(fixed = Residual ~ 1, random = ~ 1 | EQID)", color="blue")
    
    
    if (GM_model == "VIL") {
        lnSa <- function(mag, dist, vs30, b1, b2, b3, b4, b5, b6) {
        return(b1 + b2 * mag + b3 * mag ^ 2 + (b4 + b5 * mag) * 
               log(sqrt(4.5 ^ 2 + dist ^ 2)) + b6 * log(vs30/500))
        }
    }
    model2 <- nlme(obs ~ lnSa(mag, rrup, vs30, b1, b2, b3, b4, b5, b6),
                     data = Dataset2,
                     fixed = b1 + b2 + b3 + b4 + b5 + b6 ~ 1,
                     random = b1 ~ 1 | EQID,
                     start=c(b1=-6, b2=3, b3=0.0, b4=-2.6, b5=0.3, b6=-1.29))
    etaR2 <- data.frame(model2$coefficients$random)
    colnames(etaR2) <- c("L2L")
    p2 <- ggplot(etaR2, aes(x=L2L)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
                 stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(etaR2$L2L))) * nrow(etaR2) * 0.1, color = "red", size = 1) +
                 annotate("text", x=-0.8, y=5, label = paste0("phi[L2L] ==",formatC(sqrt(var(etaR2$L2L)),digits=3,format="f",flag="0")), color="red", parse=TRUE) +
                 annotate("text", x=0, y=20, label = "lnSa ~ (mag,dist,vs30,b1,b2,b3,b4,b5,b6)", color="blue") +
                 annotate("text", x=0, y=18, label = "nlme(obs ~ lnSa,random = b1 ~ 1|EQID)", color="blue")
    
    
    model3 <- lmer(Resi ~ 1 + (1|EQID), data=Dataset2, REML=TRUE)
    TotalR3 <- residuals(model3, level = 0)
    epsR3 <- residuals(model3, level = 1)
    #etaR3 <- coef(model3)$EQID
    etaR3 <- ranef(model3)$EQID
    colnames(etaR3) <- c("L2L")
    p3 <- ggplot(etaR3, aes(x=L2L)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
                stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(etaR3$L2L))) * nrow(etaR3) * 0.1, color = "red", size = 1) +
                annotate("text", x=-0.8, y=5, label = paste0("phi[L2L] ==",formatC(sqrt(var(etaR3$L2L)),digits=3,format="f",flag="0")), color="red", parse=TRUE) +
                annotate("text", x=0, y=20, label = "lmer(Residual ~ 1 + (1|EQID))", color="blue")
    
    
    lnSa.deriv = deriv(~ b1 + b2 * mag + b3 * mag ^ 2 + (b4 + b5 * mag) * log(sqrt(4.5^2 + rrup^2)) + b6 * log(vs30/500), 
                       namevec=c('b1','b2','b3','b4','b5','b6'), 
                       function.arg = c('mag','rrup','vs30','b1','b2','b3','b4','b5','b6'))
    model4 <- nlmer(obs ~ lnSa.deriv(mag, rrup, vs30, b1, b2, b3, b4, b5, b6) ~ (b1|EQID),
                    data=Dataset2, start=c(b1=-6, b2=3, b3=0.0, b4=-2.6, b5=0.3, b6=-1.29))
    etaR4 <- ranef(model4)$EQID
    colnames(etaR4) <- c("L2L")
    p4 <- ggplot(etaR4, aes(x=L2L)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
      stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(etaR4$L2L))) * nrow(etaR4) * 0.1, color = "red", size = 1) +
      annotate("text", x=-0.8, y=5, label = paste0("phi[L2L] ==",formatC(sqrt(var(etaR4$L2L)),digits=3,format="f",flag="0")), color="red", parse=TRUE) +
      annotate("text", x=0, y=20, label = "nlmer(obs ~ deriv(lnSa) ~ 1 | EQID)", color="blue")
    
    etaR_plot <- cbind(etaR, etaR2, etaR3, etaR4)
    colnames(etaR_plot) <- c("L2L","L2L_2","L2L_3","L2L_4")
    p5 <- ggplot(etaR_plot, aes(x = L2L, y = L2L_2)) + geom_point(shape=1,size=1) +
      scale_x_continuous(name="L2L from Panel A", lim=c(-1,1)) +
      scale_y_continuous(lim=c(-1,1),name="L2L from other panels")
    p5 <- p5 + geom_point(data=etaR_plot, aes(x = L2L, y = L2L_3), shape=1, size=1, color="blue") +
      geom_point(data=etaR_plot, aes(x = L2L, y = L2L_4), shape=1, size=1, color="red") +
      geom_segment(x=-1,y=-1,xend=1,yend=1,color = "grey",linetype=2) +
      annotate("text", x=-0.75, y=1, label = "Panel B", color="black") +
      annotate("text", x=-0.75, y=0.8, label = "Panel C", color="blue") +
      annotate("text", x=-0.75, y=0.6, label = "Panel D", color="red") 
    
    
    g <- plot_grid(p1,p2,p3,p4,p5, labels = "AUTO", ncol = 2)
    ggsave(file=paste(path_figures,"/Residual_L2L_",Dataname,".T",
                      formatC(Period,digits=3,width=6,format="f",flag="0"),".png",sep=""),
           g,height=12, width=8, units='in', dpi=300)
  }
  
  if (S2S_flag) {
   
      model <- lme(fixed = Resi ~ 1, data = Dataset2, method="ML", random = ~ 1 | EQID)
      TotalR = residuals(model, level = 0)
      epsR = residuals(model, level = 1)
      etaR = model$coefficients$random$EQID
      Dataset3 <- data.frame(Dataset2$STID,epsR) 
      colnames(Dataset3) <- c("STID","epsR")
      model <- lme(fixed = epsR ~ 1, data = Dataset3, method="ML", random = ~ 1 | STID)
      etaR <- data.frame(model$coefficients$random$STID)
      colnames(etaR) <- c("S2S")
      p1 <- ggplot(etaR, aes(x=S2S)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
                   stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(etaR$S2S))) * nrow(etaR) * 0.1, color = "red", size = 1) +
                   annotate("text", x=-0.8, y=5, label = paste0("phi[S2S] ==",formatC(sqrt(var(etaR$S2S)),digits=3,format="f",flag="0")), color="red", parse=TRUE) +
                   annotate("text", x=0, y=20, label = "lme(fixed = Residual ~ 1, random = ~ 1 | EQID)", color="blue") +
                   annotate("text", x=0, y=18, label = "paste(\"lme(fixed = \",delta,W[es],\" ~ 1, random = ~ 1 | STID)\")", color="blue", parse=TRUE)
      
      
      Dataset3 <- data.frame(Dataset2$STID,epsR)
      colnames(Dataset3) <- c("STID","epsR")
      etaR2 <- aggregate(Dataset3[, 2], list(Dataset3$STID), mean)
      colnames(etaR2) <- c("stid","S2S")
      p2 <- ggplot(etaR2, aes(x=S2S)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
                   stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(etaR2$S2S))) * nrow(etaR2) * 0.1, color = "red", size = 1) +
                   annotate("text", x=-0.8, y=5, label = paste0("phi[S2S] ==",formatC(sqrt(var(etaR2$S2S)),digits=3,format="f",flag="0")), color="red", parse=TRUE) +
                   annotate("text", x=0, y=20, label = "lme(fixed = Residual ~ 1, random = ~ 1 | EQID)", color="blue") +
                   annotate("text", x=0, y=18, label = "paste(S2S[s],\" = mean(\",delta,W[es],\")\")",color="blue", parse=TRUE)

      
      model3 <- lmer(Resi ~ 1 + (1|EQID), data=Dataset2)
      TotalR3 <- residuals(model3, level = 0)
      #colnames(TotalR) <- c("eqid","TotalR")
      epsR3 <- residuals(model3, level = 1)
      #colnames(epsR) <- c("eqid","epsR")
      etaR3 <- ranef(model3)$EQID
      #colnames(etaR2) <- c("etaR")
      Dataset3 <- data.frame(Dataset2$STID,epsR3)
      colnames(Dataset3) <- c("STID","epsR")
      model3 <- lmer(epsR3 ~ 1 + (1|STID), data=Dataset3)
      etaR3 <- coef(model3)$STID
      colnames(etaR3) <- c("S2S")
      p3 <- ggplot(etaR3, aes(x=S2S)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
                   stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(etaR3$S2S))) * nrow(etaR3) * 0.1, color = "red", size = 1) +
                   annotate("text", x=-0.8, y=5, label = paste0("phi[S2S] ==",formatC(sqrt(var(etaR3$S2S)),digits=3,format="f",flag="0")), color="red", parse=TRUE) +
                   annotate("text", x=0, y=20, label = "lmer(Residual ~ 1 + (1|EQID))", color="blue") +
                   annotate("text", x=0, y=18, label = "paste(\"lmer(\",delta,W[es],\" ~ 1 + (1|STID))\")", color="blue", parse=TRUE)

      
      Dataset3 <- data.frame(Dataset2$STID,epsR3)
      colnames(Dataset3) <- c("STID","epsR")
      etaR4 <- aggregate(Dataset3[, 2], list(Dataset3$STID), mean)
      colnames(etaR4) <- c("stid","S2S")
      p4 <- ggplot(etaR4, aes(x=S2S)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
        stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(etaR4$S2S))) * nrow(etaR4) * 0.1, color = "red", size = 1) +
        annotate("text", x=-0.8, y=5, label = paste0("phi[S2S] ==",formatC(sqrt(var(etaR4$S2S)),digits=3,format="f",flag="0")), color="red", parse=TRUE) +
        annotate("text", x=0, y=20, label = "lmer(Residual ~ 1 + (1|EQID))", color="blue") +
        annotate("text", x=0, y=18, label = "paste(S2S[s],\" = mean(\",delta,W[es],\")\")",color="blue", parse=TRUE)
      
      
      model5 <- lmer(Resi ~ 1 + (1|EQID) + (1|STID), data=Dataset2)
      etaR5 <- ranef(model5)$STID
      colnames(etaR5) <- c("S2S")
      p5 <- ggplot(etaR5, aes(x=S2S)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
                   stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(etaR5$S2S))) * nrow(etaR5) * 0.1, color = "red", size = 1) +
                   annotate("text", x=-0.8, y=5, label = paste0("phi[S2S] ==",formatC(sqrt(var(etaR5$S2S)),digits=3,format="f",flag="0")), color="red", parse=TRUE) +
                   annotate("text", x=0, y=20, label = "lmer(Residual ~ 1 + (1|EQID) + (1|STID))", color="blue")

      
      etaR_plot <- cbind(etaR2, etaR, etaR3, etaR5)
      colnames(etaR_plot) <- c("stid","S2S2","S2S1","S2S3","S2S5")
      p6 <- ggplot(etaR_plot, aes(x = S2S5, y = S2S1)) + geom_point(shape=1,size=1) +
                  scale_x_continuous(name="S2S from Panel E", lim=c(-1,1)) +
                  scale_y_continuous(lim=c(-1,1),name="S2S from other panels")
      p6 <- p6 + geom_point(data=etaR_plot, aes(x = S2S5, y = S2S2), shape=1, size=1, color="blue") + 
                 geom_point(data=etaR_plot, aes(x = S2S5, y = S2S3), shape=1, size=1, color="red") +
                 geom_segment(x=-1,y=-1,xend=1,yend=1,color = "grey",linetype=2) +
                 annotate("text", x=-0.75, y=1, label = "Panel A", color="black") +
                 annotate("text", x=-0.75, y=0.8, label = "Panel B", color="blue") +
                 annotate("text", x=-0.75, y=0.6, label = "Panel C", color="red")
      
      
      g <- plot_grid(p1,p2,p3,p4,p5,p6, labels = "AUTO", ncol = 2)
      ggsave(file=paste(path_figures,"/Residual_S2S_",Dataname,".T",
                        formatC(Period,digits=3,width=6,format="f",flag="0"),".png",sep=""),
             g,height=8, width=8, units='in', dpi=300)
      
      
      model <- lme(fixed = Resi ~ 1, data = Dataset2, method="ML", random = ~ 1 | EQID)
      TotalR = residuals(model, level = 0)
      epsR = residuals(model, level = 1)
      etaR = model$coefficients$random$EQID
      etaR <- cbind(data.frame(rownames(etaR)),etaR)
      colnames(etaR) <- c("EQID","L2L")
      lon <- data.frame(matrix(NA, nrow = nrow(etaR), ncol = 1))
      lat <- data.frame(matrix(NA, nrow = nrow(etaR), ncol = 1))
      
      for (i in 1:nrow(etaR)) {
        index <- etaR$EQID[i]
        lon[i,1] <- unique(Dataset[Dataset$eqid==index,]$lon)
        lat[i,1] <- unique(Dataset[Dataset$eqid==index,]$lat)
      }
      map_data <- cbind(etaR,lon,lat)
      colnames(map_data) <- c("EQID","L2L","lon","lat")
      
      ave_data <- data.frame(matrix(NA, nrow = 425, ncol = 4))
      colnames(ave_data) <- c("xx","yy","etaR","L2L")
      xc_list <- seq(from=-119.2, to=-116.8, by=0.1)
      yc_list <- seq(from=33.2, to=34.8, by=0.1)
      for (i in seq_along(xc_list)) {
        for (j in seq_along(yc_list)) {
          k = (i-1) * 17 + j
          xc <- xc_list[i]
          yc <- yc_list[j]
          ave_data[k,1] <- xc
          ave_data[k,2] <- yc
          ave_data[k,3] <- mean(map_data[map_data$lon>=(xc-0.05) & map_data$lon<(xc+0.05) & map_data$lat>=(yc-0.05) & map_data$lat<(yc+0.05),]$L2L)
          if (!is.nan(ave_data[k,3])) {
              if (ave_data[k,3]<(-0.5)) {
                ave_data[k,4] <- (-0.5) 
              } else if (ave_data[k,3]<(-0.25) & ave_data[k,3]>=(-0.5)) {
                ave_data[k,4] <- (-0.25)
              } else if (ave_data[k,3]<(-0.1) & ave_data[k,3]>=(-0.25)) {
                ave_data[k,4] <- (-0.1)
              } else if (ave_data[k,3]<(0.1) & ave_data[k,3]>=(-0.1)) {
                ave_data[k,4] <- 0.1
              } else if (ave_data[k,3]<(0.25) & ave_data[k,3]>=(0.1)) {
                ave_data[k,4] <- 0.25
              } else if (ave_data[k,3]<(0.5) & ave_data[k,3]>=(0.25)) {
                ave_data[k,4] <- 0.5
              } else if (ave_data[k,3]<(0.75) & ave_data[k,3]>=(0.5)) {
                ave_data[k,4] <- 0.75
              } else if (ave_data[k,3]>=(0.75)) {
                ave_data[k,4] <- 1
              }
          }
        }
      }
    
      
      #ca <- get_map(location = c(-119.2, 33.2, -116.8, 34.8), maptype = c("toner"),color = "bw", source = "stamen", crop=FALSE)
      #gg <- ggmap(ca)
      states <- map_data("state")
      ca_df <- subset(states, region == "california")
      gg <- ggplot(data = ca_df, mapping = aes(x = long, y = lat)) + coord_fixed(xlim = c(-119.2, -116.8),  ylim = c(33.2, 34.8), ratio=1.3) + 
        geom_polygon(color = "black", fill = "gray") 
      gg <- gg + geom_point(data=ave_data,aes(x=xx,y=yy,color=factor(L2L),size=2)) + 
                 scale_colour_manual(na.value="white", values = c("red","orange","yellow","green","cyan","blue","purple","darkblue")) +
                 annotate("text", x=-118, y=34.8, label ="paste(L2L[r],\" = mean(\",delta,B[er],\")\")",color="blue", parse=TRUE)
        #scale_colour_gradientn(limits=c(-1,1),colours = brewer.pal(n = 8, name = "RdYlBu"),na.value = "white",
                               #values = rescale(c(-1,-0.5,-0.25,-0.1,0.1,0.25,0.5,1)))
                  #scale_x_continuous(breaks = seq(from=-119.2, to=-116.8,by=0.6), labels = seq(from=-119.2, to=-116.8,by=0.6)) +
                  #scale_y_continuous(breaks = seq(from=33.2, to=34.8,by=0.2), labels = seq(from=33.2, to=34.8,by=0.2))
      ggsave(file=paste(path_figures,"/L2L_map_",Dataname,".T",
                        formatC(Period,digits=3,width=6,format="f",flag="0"),".png",sep=""),
             gg,height=8, width=8, units='in', dpi=300)
  } # End S2S flag
  
  if (P2P_flag) {
    
    ############
    model <- lme(fixed = Resi ~ 1, data = Dataset, random = ~ 1 | eqid, control = lmeControl(opt = "optim"))
    TotalR = residuals(model, level = 0)
    model_eps = residuals(model, level = 1)
    model_L2L <- data.frame(model$coefficients$random$eqid)
    colnames(model_L2L) <- c("L2L")
    Dataset3 <- data.frame(Dataset$stid,model_eps) 
    colnames(Dataset3) <- c("stid","eps")
    
    model <- lme(fixed = eps ~ 1, data = Dataset3, method="REML", random = ~ 1 | stid, control = lmeControl(opt = "optim"))
    model_eps = residuals(model, level = 1)
    model_S2S <- data.frame(model$coefficients$random$stid)
    colnames(model_S2S) <- c("S2S")
    Dataset3 <- data.frame(Dataset$pathid,model_eps) 
    colnames(Dataset3) <- c("pathid","eps")
    
    model <- lme(fixed = eps ~ 1, data = Dataset3, method="REML", random = ~ 1 | pathid, control = lmeControl(opt = "optim"))
    model_eps = residuals(model, level = 1)
    model_P2P <- data.frame(model$coefficients$random$pathid)
    colnames(model_P2P) <- c("P2P")
    
    p1_1 <- ggplot(model_L2L, aes(x=L2L)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
            stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model_L2L$L2L))) * nrow(model_L2L) * 0.1, color = "red", size = 1) + 
            annotate("text", x=-0.6, y=10, label = paste0("phi[L2L] ==",formatC(sqrt(var(model_L2L$L2L)),digits=2,format="f",flag="0")), color="red", parse=TRUE) +
            annotate("text", x=0, y=20, label = "lme(Residual ~ 1 + 1 | EQID)", color="blue", size=2.5) +
            annotate("text", x=0, y=18, label = "paste(\"lme(\",delta,W[es],\" ~ 1 + 1 | STID)\")", color="blue", parse=TRUE, size=2.5) +
            annotate("text", x=0, y=16, label = "paste(\"lme(\",delta,WS[es],\" ~ 1 + 1 | PATHID)\")", color="blue", parse=TRUE, size=2.5)
    
    p1_2 <- ggplot(model_S2S, aes(x=S2S)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
            stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model_S2S$S2S))) * nrow(model_S2S) * 0.1, color = "red", size = 1) + 
            annotate("text", x=-0.6, y=10, label = paste0("phi[S2S] ==",formatC(sqrt(var(model_S2S$S2S)),digits=2,format="f",flag="0")), color="red", parse=TRUE)
    
    p1_3 <- ggplot(model_P2P, aes(x=P2P)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
            stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model_P2P$P2P))) * nrow(model_P2P) * 0.1, color = "red", size = 1) + 
            annotate("text", x=-0.6, y=200, label = paste0("phi[P2P] ==",formatC(sqrt(var(model_P2P$P2P)),digits=2,format="f",flag="0")), color="red", parse=TRUE)
    

    ############
    
    ############
    model2 <- lmer(Resi ~ 1 + (1|eqid), data=Dataset)
    model2_eps <- residuals(model2, level = 1)
    model2_L2L <- ranef(model2)$eqid
    colnames(model2_L2L) <- c("L2L")
    Dataset3 <- data.frame(Dataset$stid,model2_eps)
    colnames(Dataset3) <- c("stid","eps")
    
    model2 <- lmer(eps ~ 1 + (1|stid), data=Dataset3)
    model2_eps <- residuals(model2, level = 1)
    model2_S2S <- ranef(model2)$stid
    colnames(model2_S2S) <- c("S2S")
    Dataset3 <- data.frame(Dataset$pathid,model2_eps)
    colnames(Dataset3) <- c("pathid","eps")
    
    model2 <- lmer(eps ~ 1 + (1|pathid), data=Dataset3)
    model2_eps <- residuals(model2, level = 1)
    model2_P2P <- ranef(model2)$pathid
    colnames(model2_P2P) <- c("P2P")
    
    p2_1 <- ggplot(model2_L2L, aes(x=L2L)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
      stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model2_L2L$L2L))) * nrow(model2_L2L) * 0.1, color = "red", size = 1) + 
      annotate("text", x=-0.6, y=10, label = paste0("phi[L2L] ==",formatC(sqrt(var(model2_L2L$L2L)),digits=2,format="f",flag="0")), color="red", parse=TRUE) +
      annotate("text", x=0, y=20, label = "lmer(Residual ~ 1 + 1 | EQID)", color="blue", size=2.5) +
      annotate("text", x=0, y=18, label = "paste(\"lmer(\",delta,W[es],\" ~ 1 + 1 | STID)\")", color="blue", parse=TRUE, size=2.5) +
      annotate("text", x=0, y=16, label = "paste(\"lmer(\",delta,WS[es],\" ~ 1 + 1 | PATHID)\")", color="blue", parse=TRUE, size=2.5)
    
    p2_2 <- ggplot(model2_S2S, aes(x=S2S)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
      stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model2_S2S$S2S))) * nrow(model2_S2S) * 0.1, color = "red", size = 1) + 
      annotate("text", x=-0.6, y=10, label = paste0("phi[S2S] ==",formatC(sqrt(var(model2_S2S$S2S)),digits=2,format="f",flag="0")), color="red", parse=TRUE)
    
    p2_3 <- ggplot(model2_P2P, aes(x=P2P)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
      stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model2_P2P$P2P))) * nrow(model2_P2P) * 0.1, color = "red", size = 1) + 
      annotate("text", x=-0.6, y=200, label = paste0("phi[P2P] ==",formatC(sqrt(var(model2_P2P$P2P)),digits=2,format="f",flag="0")), color="red", parse=TRUE)
    ############
    
    ############
    model3 <- lmer(Resi ~ 1 + (1|eqid) + (1|stid) + (1|pathid), data=Dataset)
    model3_L2L <- ranef(model3)$eqid
    colnames(model3_L2L) <- c("L2L")
    model3_S2S <- ranef(model3)$stid
    colnames(model3_S2S) <- c("S2S")
    model3_P2P <- ranef(model3)$pathid
    colnames(model3_P2P) <- c("P2P")
    
    p3_1 <- ggplot(model3_L2L, aes(x=L2L)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
      stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model3_L2L$L2L))) * nrow(model3_L2L) * 0.1, color = "red", size = 1) + 
      annotate("text", x=-0.6, y=10, label = paste0("phi[L2L] ==",formatC(sqrt(var(model3_L2L$L2L)),digits=2,format="f",flag="0")), color="red", parse=TRUE) +
      annotate("text", x=0, y=20, label = "lmer(Residual ~ 1 + 1 | EQID", color="blue", size=2.5) +
      annotate("text", x=0, y=18, label = "+ 1 | STID + 1 | PATHID)", color="blue", size=2.5) 
    
    p3_2 <- ggplot(model3_S2S, aes(x=S2S)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
      stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model3_S2S$S2S))) * nrow(model3_S2S) * 0.1, color = "red", size = 1) + 
      annotate("text", x=-0.6, y=10, label = paste0("phi[S2S] ==",formatC(sqrt(var(model3_S2S$S2S)),digits=2,format="f",flag="0")), color="red", parse=TRUE)
    
    p3_3 <- ggplot(model3_P2P, aes(x=P2P)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
      stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model3_P2P$P2P))) * nrow(model3_P2P) * 0.1, color = "red", size = 1) + 
      annotate("text", x=-0.6, y=200, label = paste0("phi[P2P] ==",formatC(sqrt(var(model3_P2P$P2P)),digits=2,format="f",flag="0")), color="red", parse=TRUE)
    

    ############
    
    ############
    if (GM_model == "VIL") {
        lnSa.deriv = deriv(~ b1 + b2 * mag + b3 * mag ^ 2 + (b4 + b5 * mag) * log(sqrt(4^2 + rrup^2)) + b6 * log(vs30/500), 
                           namevec=c('b1','b2','b3','b4','b5','b6'), 
                           function.arg = c('mag','rrup','vs30','b1','b2','b3','b4','b5','b6'))
        model4 <- nlmer(obs ~ lnSa.deriv(mag, rrup, vs30, b1, b2, b3, b4, b5, b6) ~ (b1|eqid) + (b1|stid) + (b1|pathid),
                        data=Dataset, start=c(b1=-6, b2=3, b3=0.0, b4=-2.6, b5=0.3, b6=-1.29))
        model4_L2L <- ranef(model4)$eqid
        colnames(model4_L2L) <- c("L2L")
        model4_S2S <- ranef(model4)$stid
        colnames(model4_S2S) <- c("S2S")
        model4_P2P <- ranef(model4)$pathid
        colnames(model4_P2P) <- c("P2P")
	      epsR <- residuals(model4, level = 1)
   
        p4_1 <- ggplot(model4_L2L, aes(x=L2L)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
          stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model4_L2L$L2L))) * nrow(model4_L2L) * 0.1, color = "red", size = 1) + 
          annotate("text", x=-0.6, y=10, label = paste0("phi[L2L] ==",formatC(sqrt(var(model4_L2L$L2L)),digits=2,format="f",flag="0")), color="red", parse=TRUE) +
          annotate("text", x=0, y=20, label = "nlmer(obs ~ deriv(lnSa) ~ 1 | EQID", color="blue", size=2.5) +
          annotate("text", x=0, y=18, label = "+ 1 | STID + 1 | PATHID)", color="blue", size=2.5) 
        
        p4_2 <- ggplot(model4_S2S, aes(x=S2S)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
          stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model4_S2S$S2S))) * nrow(model4_S2S) * 0.1, color = "red", size = 1) + 
          annotate("text", x=-0.6, y=10, label = paste0("phi[S2S] ==",formatC(sqrt(var(model4_S2S$S2S)),digits=2,format="f",flag="0")), color="red", parse=TRUE)
        
        p4_3 <- ggplot(model4_P2P, aes(x=P2P)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=c(-1,1)) +
          stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model4_P2P$P2P))) * nrow(model4_P2P) * 0.1, color = "red", size = 1) + 
          annotate("text", x=-0.6, y=200, label = paste0("phi[P2P] ==",formatC(sqrt(var(model4_P2P$P2P)),digits=2,format="f",flag="0")), color="red", parse=TRUE)
        
        g <- plot_grid(p1_1,p2_1,p3_1,p4_1,p1_2,p2_2,p3_2,p4_2,p1_3,p2_3,p3_3,p4_3, labels = "AUTO", nrow = 3, align = 'v')
        ggsave(file=paste(path_figures,"/Residual_P2P_",Dataname,".T", formatC(Period,digits=3,width=6,format="f",flag="0"),".pdf",sep=""), 
               g,height=9, width=12, units='in', dpi=300)
    } else if (GM_model=="ASK14") {
        g <- plot_grid(p1_1,p2_1,p3_1,p1_2,p2_2,p3_2,p1_3,p2_3,p3_3, labels = "AUTO", nrow = 3, align = 'v')
        ggsave(file=paste(path_figures,"/Residual_P2P_",Dataname,".T", formatC(Period,digits=3,width=6,format="f",flag="0"),".pdf",sep=""), 
               g,height=6, width=6, units='in', dpi=300)
    } # End GM model selection
    ############
    
    ############
    if (GM_model == "VIL") {
        L2L_plot <- cbind(model4_L2L, model2_L2L,model3_L2L)
        colnames(L2L_plot) <- c("L2L_4","L2L_2","L2L_3")
        p5_1 <- ggplot(L2L_plot, aes(x = L2L_4, y = L2L_2)) + geom_point(shape=1,size=1) +
                scale_x_continuous(name="L2L from model 4", lim=c(-1,1)) +
                scale_y_continuous(lim=c(-1,1),name="L2L from model 2, 3")
        p5_1 <- p5_1 + geom_point(data=L2L_plot, aes(x = L2L_4, y = L2L_3), shape=1, size=1, color="red") +
                       geom_segment(x=-1,y=-1,xend=1,yend=1,color = "grey",linetype=2) +
                       annotate("text", x=-0.75, y=1, label = "model 2", color="black") +
                       annotate("text", x=-0.75, y=0.6, label = "model 3", color="red")
        
        S2S_plot <- cbind(model4_S2S, model2_S2S,model3_S2S)
        colnames(S2S_plot) <- c("S2S_4","S2S_2","S2S_3")
        p5_2 <- ggplot(S2S_plot, aes(x = S2S_4, y = S2S_2)) + geom_point(shape=1,size=1) +
          scale_x_continuous(name="S2S from model 4", lim=c(-1,1)) +
          scale_y_continuous(lim=c(-1,1),name="S2S from model 2, 3")
        p5_2 <- p5_2 + geom_point(data=S2S_plot, aes(x = S2S_4, y = S2S_3), shape=1, size=1, color="red") +
          geom_segment(x=-1,y=-1,xend=1,yend=1,color = "grey",linetype=2) +
          annotate("text", x=-0.75, y=1, label = "model 2", color="black") +
          annotate("text", x=-0.75, y=0.6, label = "model 3", color="red")
        
        P2P_plot <- cbind(model4_P2P, model2_P2P,model3_P2P)
        colnames(P2P_plot) <- c("P2P_4","P2P_2","P2P_3")
        p5_3 <- ggplot(P2P_plot, aes(x = P2P_4, y = P2P_2)) + geom_point(shape=1,size=1) +
          scale_x_continuous(name="P2P from model 4", lim=c(-1,1)) +
          scale_y_continuous(lim=c(-1,1),name="P2P from model 2, 3")
        p5_3 <- p5_3 + geom_point(data=P2P_plot, aes(x = P2P_4, y = P2P_3), shape=1, size=1, color="red") +
          geom_segment(x=-1,y=-1,xend=1,yend=1,color = "grey",linetype=2) +
          annotate("text", x=-0.75, y=1, label = "model 2", color="black") +
          annotate("text", x=-0.75, y=0.6, label = "model 3", color="red")
        
        g <- plot_grid(p5_1,p5_2,p5_3, labels = "AUTO", nrow = 2, align = 'v')
        ggsave(file=paste(path_figures,"/Compare_P2P_",Dataname,".T", formatC(Period,digits=3,width=6,format="f",flag="0"),".pdf",sep=""), 
               g,height=6, width=6, units='in', dpi=300)
    } # End GM model selection
    ############
    
    ############
    if (GM_model == "VIL") {
        lon <- data.frame(matrix(NA, nrow = nrow(model_L2L), ncol = 1))
        lat <- data.frame(matrix(NA, nrow = nrow(model_L2L), ncol = 1))
        for (i in 1:nrow(model_L2L)) {
          index <- rownames(model_L2L)[i]
          lon[i,1] <- unique(Dataset[Dataset$eqid==index,]$lon)
          lat[i,1] <- unique(Dataset[Dataset$eqid==index,]$lat)
        }
        xc_list <- seq(from=-119.2, to=-116.8, by=0.1)
        yc_list <- seq(from=33.2, to=34.8, by=0.1)
         
        # map_data <- cbind(model_L2L,lon,lat)
        # colnames(map_data) <- c("L2L","lon","lat")
        # ave_data <- data.frame(matrix(NA, nrow = 425, ncol = 4))
        # colnames(ave_data) <- c("xx","yy","L2L","Ave")
        # for (i in seq_along(xc_list)) {
        #   for (j in seq_along(yc_list)) {
        #     k = (i-1) * 17 + j
        #     xc <- xc_list[i]
        #     yc <- yc_list[j]
        #     ave_data[k,1] <- xc
        #     ave_data[k,2] <- yc
        #     ave_data[k,3] <- mean(map_data[map_data$lon>=(xc-0.05) & map_data$lon<(xc+0.05) & map_data$lat>=(yc-0.05) & map_data$lat<(yc+0.05),]$L2L)
        #     if (!is.nan(ave_data[k,3])) {
        #       if (ave_data[k,3]<(-0.5)) {
        #         ave_data[k,4] <- (-0.5) 
        #       } else if (ave_data[k,3]<(-0.25) & ave_data[k,3]>=(-0.5)) {
        #         ave_data[k,4] <- (-0.25)
        #       } else if (ave_data[k,3]<(-0.1) & ave_data[k,3]>=(-0.25)) {
        #         ave_data[k,4] <- (-0.1)
        #       } else if (ave_data[k,3]<(0.1) & ave_data[k,3]>=(-0.1)) {
        #         ave_data[k,4] <- 0.1
        #       } else if (ave_data[k,3]<(0.25) & ave_data[k,3]>=(0.1)) {
        #         ave_data[k,4] <- 0.25
        #       } else if (ave_data[k,3]<(0.5) & ave_data[k,3]>=(0.25)) {
        #         ave_data[k,4] <- 0.5
        #       } else if (ave_data[k,3]<(0.75) & ave_data[k,3]>=(0.5)) {
        #         ave_data[k,4] <- 0.75
        #       } else if (ave_data[k,3]>=(0.75)) {
        #         ave_data[k,4] <- 1
        #       }
        #     }
        #   }
        # }
        # states <- map_data("state")
        # ca_df <- subset(states, region == "california")
        # p1 <- ggplot(data = ca_df, mapping = aes(x = long, y = lat)) + coord_fixed(xlim = c(-119.2, -116.8),  ylim = c(33.2, 34.8), ratio=1.3) + 
        #   geom_polygon(color = "black", fill = "gray") 
        # p1 <- p1 + geom_point(data=ave_data,aes(x=xx,y=yy,color=factor(Ave)),size=2) + 
        #   scale_colour_manual(na.value="white", values = c("red","orange","yellow","green","cyan","blue","purple","darkblue"),guide=FALSE) +
        #   annotate("text", x=-118, y=34.8, label ="Model 1",color="blue")
        # 
        # map_data <- cbind(model3_L2L,lon,lat)
        # colnames(map_data) <- c("L2L","lon","lat")
        # ave_data <- data.frame(matrix(NA, nrow = 425, ncol = 4))
        # colnames(ave_data) <- c("xx","yy","L2L","Ave")
        # for (i in seq_along(xc_list)) {
        #   for (j in seq_along(yc_list)) {
        #     k = (i-1) * 17 + j
        #     xc <- xc_list[i]
        #     yc <- yc_list[j]
        #     ave_data[k,1] <- xc
        #     ave_data[k,2] <- yc
        #     ave_data[k,3] <- mean(map_data[map_data$lon>=(xc-0.05) & map_data$lon<(xc+0.05) & map_data$lat>=(yc-0.05) & map_data$lat<(yc+0.05),]$L2L)
        #     if (!is.nan(ave_data[k,3])) {
        #       if (ave_data[k,3]<(-0.5)) {
        #         ave_data[k,4] <- (-0.5) 
        #       } else if (ave_data[k,3]<(-0.25) & ave_data[k,3]>=(-0.5)) {
        #         ave_data[k,4] <- (-0.25)
        #       } else if (ave_data[k,3]<(-0.1) & ave_data[k,3]>=(-0.25)) {
        #         ave_data[k,4] <- (-0.1)
        #       } else if (ave_data[k,3]<(0.1) & ave_data[k,3]>=(-0.1)) {
        #         ave_data[k,4] <- 0.1
        #       } else if (ave_data[k,3]<(0.25) & ave_data[k,3]>=(0.1)) {
        #         ave_data[k,4] <- 0.25
        #       } else if (ave_data[k,3]<(0.5) & ave_data[k,3]>=(0.25)) {
        #         ave_data[k,4] <- 0.5
        #       } else if (ave_data[k,3]<(0.75) & ave_data[k,3]>=(0.5)) {
        #         ave_data[k,4] <- 0.75
        #       } else if (ave_data[k,3]>=(0.75)) {
        #         ave_data[k,4] <- 1
        #       }
        #     }
        #   }
        # }
        # states <- map_data("state")
        # ca_df <- subset(states, region == "california")
        # p2 <- ggplot(data = ca_df, mapping = aes(x = long, y = lat)) + coord_fixed(xlim = c(-119.2, -116.8),  ylim = c(33.2, 34.8), ratio=1.3) + 
        #   geom_polygon(color = "black", fill = "gray") 
        # p2 <- p2 + geom_point(data=ave_data,aes(x=xx,y=yy,color=factor(Ave)),size=2) + 
        #   scale_colour_manual(na.value="white", values = c("red","orange","yellow","green","cyan","blue","purple","darkblue"),guide=FALSE) +
        #   annotate("text", x=-118, y=34.8, label = "Model 3",color="blue")
        
        map_data <- cbind(model4_L2L,lon,lat)
        colnames(map_data) <- c("L2L","lon","lat")
        ave_data <- data.frame(matrix(NA, nrow = 425, ncol = 4))
        colnames(ave_data) <- c("xx","yy","L2L","Ave")
        for (i in seq_along(xc_list)) {
          for (j in seq_along(yc_list)) {
            k = (i-1) * 17 + j
            xc <- xc_list[i]
            yc <- yc_list[j]
            ave_data[k,1] <- xc
            ave_data[k,2] <- yc
            ave_data[k,3] <- mean(map_data[map_data$lon>=(xc-0.05) & map_data$lon<(xc+0.05) & map_data$lat>=(yc-0.05) & map_data$lat<(yc+0.05),]$L2L)
            if (!is.nan(ave_data[k,3])) {
              if (ave_data[k,3]<(-0.5)) {
                ave_data[k,4] <- (-0.5) 
              } else if (ave_data[k,3]<(-0.25) & ave_data[k,3]>=(-0.5)) {
                ave_data[k,4] <- (-0.25)
              } else if (ave_data[k,3]<(-0.1) & ave_data[k,3]>=(-0.25)) {
                ave_data[k,4] <- (-0.1)
              } else if (ave_data[k,3]<(0.1) & ave_data[k,3]>=(-0.1)) {
                ave_data[k,4] <- 0.1
              } else if (ave_data[k,3]<(0.25) & ave_data[k,3]>=(0.1)) {
                ave_data[k,4] <- 0.25
              } else if (ave_data[k,3]<(0.5) & ave_data[k,3]>=(0.25)) {
                ave_data[k,4] <- 0.5
              } else if (ave_data[k,3]<(0.75) & ave_data[k,3]>=(0.5)) {
                ave_data[k,4] <- 0.75
              } else if (ave_data[k,3]>=(0.75)) {
                ave_data[k,4] <- 1
              }
            }
          }
        }
        states <- map_data("state")
        ca_df <- subset(states, region == "california")
        p3 <- ggplot(data = ca_df, mapping = aes(x = long, y = lat)) + coord_fixed(xlim = c(-119.2, -116.8),  ylim = c(33.2, 34.8), ratio=1.3) + 
          geom_polygon(color = "black", fill = "gray") 
        p3 <- p3 + geom_point(data=ave_data,aes(x=xx,y=yy,color=factor(Ave)),size=2) + 
          scale_colour_manual(na.value="white", values = c("red","orange","yellow","green","cyan","blue","purple","darkblue"),guide=FALSE) +
          annotate("text", x=-118, y=34.8, label ="Model 4",color="blue")
        
        #g <- plot_grid(p1,p2,p3, labels = "AUTO", nrow = 2, align = 'v')
        ggsave(file=paste(path_figures,"/L2L_map_",Dataname,".T", formatC(Period,digits=3,width=6,format="f",flag="0"),".pdf",sep=""),
               p3,height=6, width=6, units='in', dpi=300)
    } # End GM model selection
    ############
    
    ############
    if (GM_model == "VIL") {
        ptid <- data.frame(as.integer(rownames(model4_P2P)))
        colnames(ptid) <- c("ptid")
        P2P_map <- cbind(model4_P2P, ptid)
        P2P_map <- P2P_map[P2P_map$ptid > 298000 & substr(P2P_map$ptid,1,3)==298,]
        rgid <- data.frame(as.numeric(substr(P2P_map$ptid,4,6)))
        colnames(rgid) <- c("rgid")
        map_data <- cbind(P2P_map, rgid)
        ave_data <- data.frame(matrix(NA, nrow = 825, ncol = 3))
        colnames(ave_data) <- c("xx","yy","P2P")
        xc_list <- seq(from=-119.2, to=-116.0, by=0.1)
        yc_list <- seq(from=32.6, to=35.0, by=0.1)
        for (i in seq_along(xc_list)) {
          for (j in seq_along(yc_list)) {
            k = (i-1) * 25 + j
            xc <- xc_list[i]
            yc <- yc_list[j]
            ave_data[k,1] <- xc
            ave_data[k,2] <- yc
            if (length(map_data[map_data$rgid==k,]$P2P) == 1) {
              ave_data[k,3] <- map_data[map_data$rgid==k,]$P2P
              if (ave_data[k,3]<(-0.5)) {
                ave_data[k,3] <- (-0.5) 
              } else if (ave_data[k,3]<(-0.25) & ave_data[k,3]>=(-0.5)) {
                ave_data[k,3] <- (-0.25)
              } else if (ave_data[k,3]<(-0.1) & ave_data[k,3]>=(-0.25)) {
                ave_data[k,3] <- (-0.1)
              } else if (ave_data[k,3]<(0.1) & ave_data[k,3]>=(-0.1)) {
                ave_data[k,3] <- 0.1
              } else if (ave_data[k,3]<(0.25) & ave_data[k,3]>=(0.1)) {
                ave_data[k,3] <- 0.25
              } else if (ave_data[k,3]<(0.5) & ave_data[k,3]>=(0.25)) {
                ave_data[k,3] <- 0.5
              } else if (ave_data[k,3]<(0.75) & ave_data[k,3]>=(0.5)) {
                ave_data[k,3] <- 0.75
              } else if (ave_data[k,3]>=(0.75)) {
                ave_data[k,3] <- 1
              }
            }
          }
        }
        states <- map_data("state")
        ca_df <- subset(states, region == "california")
        p1 <- ggplot(data = ca_df, mapping = aes(x = long, y = lat)) + coord_fixed(xlim = c(-119.2, -116.8),  ylim = c(33.2, 34.8), ratio=1.3) + 
          geom_polygon(color = "black", fill = "gray") 
        p1 <- p1 + geom_point(data=ave_data,aes(x=xx,y=yy,color=factor(P2P)),size=2) + 
          scale_colour_manual(na.value="white", values = c("red","orange","yellow","green","cyan","blue","purple","darkblue"),guide=FALSE) +
          annotate("text", x=-118, y=34.8, label ="Model 1",color="blue")
        
        g <- plot_grid(p1, labels = "AUTO", nrow = 2, align = 'v')
        ggsave(file=paste(path_figures,"/P2P_map_",Dataname,".T", formatC(Period,digits=3,width=6,format="f",flag="0"),".pdf",sep=""),
               g,height=6, width=6, units='in', dpi=300)
    } # End GM model selection
    ############
  } # End P2P 
  
  
  
    
  if (residual_flag) {
        # Mixed effects analysis
        
        model <- lme(fixed = Resi ~ 1, data = Dataset2, method="ML", random = ~ 1 | EQID)
        TotalR = residuals(model, level = 0)
        epsR = residuals(model, level = 1)
        etaR = model$coefficients$random$EQID
        
    if (output_flag) {    
      write.csv(TotalR,file=paste(path_outputs,"/Total_resids.",Dataname, ".T",
                                  formatC(Period,digits=3,width=6,format="f",flag="0"),".csv", sep=""), 
                row.names=TRUE)
      write.csv(epsR,file=paste(path_outputs,"/Intra_resids.",Dataname, ".T",
                                formatC(Period,digits=3,width=6,format="f",flag="0"),".csv", sep=""), 
                row.names=TRUE)
      write.csv(etaR,file=paste(path_outputs,"/Inter_resids.",Dataname, ".T",
                                formatC(Period,digits=3,width=6,format="f",flag="0"),".csv", sep=""), 
                row.names=TRUE)
    }
  }
  
  
  ### Plot figures if needed
  if (plot_residual) {
    source(file.path(path_main,"Plot_figures.R"))
    Intra_residual_VIL(Dataname, Dataset, epsR, path_figures, Period)
    Inter_residual_VIL(Dataname, Dataset, model4_L2L, model4_S2S, model4_P2P, path_figures, Period)
    Compare_residual_Villani(Dataname, Dataset, model4_L2L, epsR, path_figures, Period)
  }
  # End plot figures
  Dataset <- subset(Dataset, select = -c(Sa1180))
} 
#End loop over all periods


### Plot more figures if needed
if (plot_response) {
  source(file.path(path_main,"Plot_figures.R"))
  #Plot response spectrum
  mag_list <- c(5,6,7,8)
  w <- 10
  dip <- 90
  ztor <- 0
  sof <- 0
  hwflag <- 0
  ry0 <- 0
  sa1180 <- 0
  z1 <- 0
  region <- 0
  
  vs30 <- 760
  vs30_2 <- 270
  rjb <- 30
  rrup <- rjb
  rx <- rjb
  
  plot_data <- data.frame(coeff[coeff$Period>0,]$Period)
  plot_data2 <- data.frame(coeff[coeff$Period>0,]$Period)
  resp <- data.frame(matrix(1, nrow = nrow(coeff[coeff$Period>0,]), ncol = length(mag_list)))
  resp2 <- data.frame(matrix(1, nrow = nrow(coeff[coeff$Period>0,]), ncol = length(mag_list)))
  plot_data <- cbind(plot_data,resp)
  plot_data2 <- cbind(plot_data2,resp)
  colnames(plot_data)[1] <- c("Period")
  colnames(plot_data2)[1] <- c("Period")
  for (k in 1:4) {
    mag <- mag_list[k]
    colnames(plot_data)[k+1] <- paste("M",mag,sep="")
    colnames(plot_data2)[k+1] <- paste("M",mag,sep="")
    for (i in 1:nrow(coeff[coeff$Period>0,])) {
      for (j in 1:ncol(coeff)) {
        assign(colnames(coeff)[j], coeff[i,j])
      }
      if (GM_model == "ASK14") {
        plot_data[i,k+1] <- exp(lnSa(mag,w,dip,ztor,sof,hwflag,rrup,rjb,rx,ry0,vs30,sa1180,z1,region))
        plot_data2[i,k+1] <- exp(lnSa(mag,w,dip,ztor,sof,hwflag,rrup,rjb,rx,ry0,vs30_2,sa1180,z1,region))
      } else if (GM_model=="KOT") {
        plot_data[i,k+1] <- exp(lnSa(mag,rrup))
        plot_data2[i,k+1] <- exp(lnSa(mag,rrup))
      } else if (GM_model=="VIL") {
        plot_data[i,k+1] <- exp(lnSa(mag,rrup,vs30))
        plot_data2[i,k+1] <- exp(lnSa(mag,rrup,vs30_2))
      }
    }
  }
  Response_spectrum(Dataname,path_figures,plot_data,plot_data2,vs30,vs30_2,rjb)
  
  # Plot distance and magnitude scaling at desired period
  desired_period <- 2
  dist <- c(0:400)
  mag_list <- c(5,6,7,8) 
  ztor <- 0
  vs30 <- 760
  rjb <- 30
  rx <- 30
  sof <- 0
  z1 <- 0
  region <- 0
  hwflag <- 0
  ry0 <- 0
  sa1180 <- 0
  w <- 10
  dip <- 90
  
  plot_data <- data.frame(dist)
  resp <- data.frame(matrix(NA, nrow = length(dist), ncol = length(mag_list)))
  plot_data <- cbind(plot_data,resp)
  colnames(plot_data)[1] <- c("Distance")
  
  for (j in 1:ncol(coeff)) {
    assign(colnames(coeff)[j], coeff[coeff$Period==desired_period,][j])
  }
  for (k in 1:4) {
    mag <- mag_list[k]
    colnames(plot_data)[k+1] <- paste("M",mag,sep="")
    for (i in 1:length(dist)) {
      if (GM_model=="ASK14") {
        plot_data[i,k+1] <- exp(lnSa(mag,w,dip,ztor,sof,hwflag,dist[i],rjb,rx,ry0,vs30,sa1180,z1,region))
      } else if (GM_model=="KOT") {
        plot_data[i,k+1] <- exp(lnSa(mag,dist[i]))
      } else if (GM_model=="VIL") {
        plot_data[i,k+1] <- exp(lnSa(mag,dist[i],vs30))
        plot_data2[i,k+1] <- exp(lnSa(mag,dist[i],vs30))
      }
    }
  }
  
  
  rjb_list <- c(1,30,200)
  mag_list <- seq(from=4, to=9,by=0.0125)
  
  plot_data2 <- data.frame(mag_list)
  resp2 <- data.frame(matrix(NA, nrow = length(mag_list), ncol = length(rjb_list)))
  plot_data2 <- cbind(plot_data2,resp2)
  colnames(plot_data2)[1] <- c("Magnitude")
  
  for (k in 1:3) {
    rjb <- rjb_list[k]
    colnames(plot_data2)[k+1] <- paste("T",rjb,sep="")
    for (i in 1:length(mag_list)) {
      if (GM_model=="ASK14") {
        plot_data2[i,k+1] <- exp(lnSa(mag_list[i],w,dip,ztor,sof,hwflag,rjb,rjb,rjb,ry0,vs30,sa1180,z1,region)) 
      } else if (GM_model=="KOT") {
        plot_data2[i,k+1] <- exp(lnSa(mag_list[i],rjb))
      } else if (GM_model=="VIL") {
        plot_data2[i,k+1] <- exp(lnSa(mag_list[i],rjb,vs30))
      }
    }
  }
  Mag_Dist_Scaling(Dataname,path_figures,plot_data,plot_data2,vs30,desired_period)
  
  ASK14_data <- read.csv(paste0(path_model,"/",GM_model,"/",Dataname,"/Distance_Scaling_vs760_t2.000s.csv"))
  ASK14_data2 <- read.csv(paste0(path_model,"/",GM_model,"/",Dataname,"/Magnitude_Scaling_vs760_t2.000s.csv"))
  Compare_Mag_Dist_Scaling(Dataname,path_figures,plot_data,plot_data2,vs30,desired_period,ASK14_data,ASK14_data2)
} # End plot more figures
    



