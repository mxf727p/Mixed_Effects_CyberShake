#########################
### Residuals Analysis For Predicted Ground Motion Equations
### 03/29/2018   Started by Xiaofeng Meng, only ASK14 and a simple model 
###              by Albert Kottke (i.e, KOT) is implemented now.

# Clear the console
cat("\014")  
# Clear the workspace
rm(list=ls())


library(ggplot2)
library(cowplot)
library(nlme)
library(lme4)
library(minpack.lm)
library(scales)
library(tictoc)


tic("total")

########################
### Set up Dataset, Model, Working directories, Plotting flags, Debug flag etc.
Dataname <- "CyberShake_SC_regression_old"     # Options: "NGA_W2", "CyberShake", "NGA_W2_SC", "CyberShake_SC_regression", "CyberShake_SC_full"
lon_e <- -115
lon_w <- -122
lat_s <- 31
lat_n <- 37
GM_model <- "VIL"        # Options: "ASK14", "KOT", "VIL"
ME_model <- 1
path_main <- setwd("/auto/scec-00/xiaofenm/CyberShake")
regression_flag <- 1
if (regression_flag) {
  regression_model <- 4   # 1: nlme; 2: nls; 3: nlsLM; 4: nlmer
}
output_flag <- 0
plot_residual <- 0
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


tic("Read data")
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
toc()


#########################
### Read equations for the model
source(file.path(path_model,GM_model,"Equations.R"))
### End read equations
#########################


#########################
### Plot map of events
  states <- map_data("state")
  ca_df <- subset(states, region == "california")
  map <- ggplot(data = ca_df, mapping = aes(x = long, y = lat)) + coord_fixed(xlim = c(lon_w, lon_e),  ylim = c(lat_s, lat_n), ratio=1.3) + geom_polygon(color = "black", fill = "gray")
  map <- map + geom_point(data=Dataset,aes(x=lon,y=lat,size=mag), shape=21)
  ggsave(file=paste(path_figures,"/Map_",Dataname,".png",sep=""),map,height=6, width=6, units='in', dpi=600)	
### End plot map
#########################


tic("Regression")
#########################
### Regression analysis
if (regression_flag) {
  #regression_list <- data.frame(c(0.010,0.020,0.030,0.050,0.075,0.100,0.150,0.200,0.250,0.300,
  #                                0.400,0.500,0.750,1.000,1.500,2.000,3.000,4.000,5.000,6.000,
  #                                7.500,10.000))
  regression_list <- data.frame(c(3.000))
  if (GM_model == "KOT") {
    colnames(regression_list) <- c("Period")
    coeff <- data.frame(matrix(NA, nrow = 1, ncol = 7))
    coeff <- coeff[-1,]
    colnames(coeff) <- c("c1","c2","c4","c6","c7","c10","c11")
    for (i in 1:nrow(regression_list)) {
      Period <- regression_list[i,1]
      obs_regression <- log(Dataset[colnames(Dataset)==format(round(Period,3),nsmall=3)])
      colnames(obs_regression) <- c("obs")
      Dataset <- cbind(Dataset, obs_regression)
      if (regression_model == 1) {
        fitModel <- nlme(obs ~ lnSa2(mag, rrup, c1, c2, c4, c6, c7, c10, c11),
                         data = Dataset[Dataset$obs!="NaN",],
                         fixed = c1 + c2 + c4 + c6 + c7 + c10 + c11 ~ 1,
                         random = c1 ~ 1 | eqid,
                         start=c(c1=-6., c2=0.9, c4=0.35, c6=-3.0, c7=0.3, c10=-0.25,c11=-0.004))
        coeff <- rbind(coeff,t(fitModel$coefficients$fixed))
      } else if (regression_model == 2) {
        fitModel <- nls(obs ~ lnSa2(mag, rrup, c1, c2, c4, c6, c7, c10, c11),
                        data = Dataset[Dataset$obs!="NaN",],
                        start=list(c1=-6., c2=0.9, c4=0.35, c6=-3.0, c7=0.3, c10=-0.25,c11=-0.004))
        coeff <- rbind(coeff,t(coef(fitModel)))
      } else if (regression_model == 3) {
        fitModel <- nlsLM(obs ~ lnSa2(mag, rrup, c1, c2, c4, c6, c7, c10, c11),
                          data = Dataset[Dataset$obs!="NaN",])
        coeff <- rbind(coeff,t(coef(fitModel)))
      } else if (regression_model == 4) {
        lnSa.deriv = deriv(~ c1 + c2 * mag + (c6 + c7 * mag) * log(rrup + exp(c4)) + c10 * (mag - 6) ^ 2 + c11 * rrup, 
                           namevec=c('c1','c2','c4','c6','c7','c10','c11'), 
                           function.arg = c('mag','rrup','c1','c2','c4','c6','c7','c10','c11'))
        fitModel <- nlmer(obs ~ lnSa.deriv(mag,rrup,c1,c2,c4,c6,c7,c10,c11) ~ (c1|eqid) + (c1|stid) + (c1|pathid),
                          data=Dataset[Dataset$obs!="NaN",], start=c(c1=-6, c2=0.9, c4=0.35, c6=-3.0, c7=0.3, c10=-0.25,c11=-0.004))
        coeff <- rbind(coeff,t(coef(summary(fitModel))[,"Estimate"]))
      }
      Dataset <- subset(Dataset, select = -c(obs))
    }
    coeff <- cbind(regression_list,coeff)
    write.csv(coeff,file=file.path(path_model,GM_model,Dataname,"Coeff.csv"))
  } else if (GM_model == "VIL") {
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
        fitModel <- nlme(obs ~ lnSa2(mag, rrup, vs30, b1, b2, b3, b4, b5, b6),
                         data = Dataset[Dataset$obs!="NaN",],
                         fixed = b1 + b2 + b3 + b4 + b5 + b6 ~ 1,
                         random = b1 ~ 1|eqid,
                         start=c(b1=-6, b2=3, b3=0.0, b4=-2.6, b5=0.3, b6=-1.29))
        coeff <- rbind(coeff,t(fitModel$coefficients$fixed))
      } else if (regression_model == 2) {
        fitModel <- nls(obs ~ lnSa2(mag, rrup, vs30, b1, b2, b3, b4, b5, b6),
                       data = Dataset[Dataset$obs!="NaN",],
                       start=list(b1=-6, b2=3, b3=0.0, b4=-2.6, b5=0.3, b6=-1.29))
        coeff <- rbind(coeff,t(coef(fitModel)))
      } else if (regression_model == 3) {
        fitModel <- nlsLM(obs ~ lnSa2(mag, rrup, vs30, b1, b2, b3, b4, b5, b6),
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
  }  # End select model
}  ### End regression analysis
#########################
toc()



#########################
### Predict ground motions and compute residuals
coeff <- Read_coeff(file.path(path_model,GM_model,Dataname,"Coeff.csv"))
PGM <- data.frame(matrix(NA, nrow = nrow(Dataset), ncol = nrow(coeff)))
Obs <- data.frame(matrix(NA, nrow = nrow(Dataset), ncol = 1))
tic("Compute residual")
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
  Obs <- log(Dataset[colnames(Dataset)==format(round(coeff$Period[i],3),nsmall=3)]) 
  
  # Compute residuals
  Period <- coeff$Period[i]
  colnames(Obs) <- c("obs")
  Resi <- Obs - PGM[,i]
  colnames(Resi) <- c("Resi")
  Dataset <- cbind(Dataset,Resi,Obs)
  Dataset <- Dataset[Dataset$Resi!="NaN",]
  Dataset <- subset(Dataset, select = -c(Sa1180))
} #End loop over all periods for residuals
toc()

  
if (ME_model == 1) {
      tic("model")
      model <- lme(fixed = Resi ~ 1, data = Dataset, random = ~ 1 | eqid, control = lmeControl(opt = "optim"))
      model_eps = residuals(model, level = 1)
      model_L2L <- data.frame(model$coefficients$random$eqid)
      colnames(model_L2L) <- c("L2L")
      Dataset3 <- data.frame(Dataset$stid,model_eps)
      colnames(Dataset3) <- c("stid","eps")
  
      model <- lme(fixed = eps ~ 1, data = Dataset3, method="REML", random = ~ 1 | stid, control = lmeControl(opt = "optim"))
      model_eps2 = residuals(model, level = 1)
      model_S2S <- data.frame(model$coefficients$random$stid)
      colnames(model_S2S) <- c("S2S")
      Dataset3 <- data.frame(Dataset$pathid,model_eps2)
      colnames(Dataset3) <- c("pathid","eps")
  
      model <- lme(fixed = eps ~ 1, data = Dataset3, method="REML", random = ~ 1 | pathid, control = lmeControl(opt = "optim"))
      epsR = residuals(model, level = 1)
      model_P2P <- data.frame(model$coefficients$random$pathid)
      colnames(model_P2P) <- c("P2P")
      toc()
} else if (ME_model == 4) {
      tic("model4")
      model4 <- lmer(Resi ~ 1 + (1|eqid) + (1|stid) + (1|pathid), data=Dataset) 
      model_L2L <- ranef(model4)$eqid
      colnames(model_L2L) <- c("L2L")
      model_S2S <- ranef(model4)$stid
      colnames(model_S2S) <- c("S2S")
      model_P2P <- ranef(model4)$pathid
      colnames(model_P2P) <- c("P2P")
      epsR <- residuals(model4, level = 1)
      toc()
} else if (ME_model == 5) {
      tic("model5")
      lnSa.deriv = deriv(~ b1 + b2 * mag + b3 * mag ^ 2 + (b4 + b5 * mag) * log(sqrt(4^2 + rrup^2)) + b6 * log(vs30/500), 
                         namevec=c('b1','b2','b3','b4','b5','b6'), 
                         function.arg = c('mag','rrup','vs30','b1','b2','b3','b4','b5','b6'))
      model4 <- nlmer(obs ~ lnSa.deriv(mag, rrup, vs30, b1, b2, b3, b4, b5, b6) ~ (b1|eqid) + (b1|stid) + (b1|pathid),
                      data=Dataset, start=c(b1=-6, b2=3, b3=0.0, b4=-2.6, b5=0.3, b6=-1.29))
      model_L2L <- ranef(model4)$eqid
      colnames(model_L2L) <- c("L2L")
      model_S2S <- ranef(model4)$stid
      colnames(model_S2S) <- c("S2S")
      model_P2P <- ranef(model4)$pathid
      colnames(model_P2P) <- c("P2P")
      epsR <- residuals(model4, level = 1)
      toc()
}



### Plot figures if needed
if (plot_residual) {
  source(file.path(path_main,"Plot_figures.R"))
  tic("Residual P2P")
  x_axis_range <- c(-1.5, 1.5)
  ymax <- max(hist(model_L2L$L2L, breaks = seq(from=-15,to=15,by=0.1), plot=FALSE)$counts)
  p1_1 <- ggplot(model_L2L, aes(x=L2L)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=x_axis_range) +
          stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model_L2L$L2L))) * nrow(model_L2L) * 0.1, color = "red", size = 1) +
          annotate("text", x=-1, y=ymax*0.1, label = paste0("phi[L2L] ==",formatC(sqrt(var(model_L2L$L2L)),digits=2,format="f",flag="0")), color="red", parse=TRUE) +
          annotate("text", x=0, y=ymax*1.0, label = "lme(Residual ~ 1 + 1 | EQID)", color="blue", size=3.5) +
          annotate("text", x=0, y=ymax*0.9, label = "paste(\"lme(\",delta,W[es],\" ~ 1 + 1 | STID)\")", color="blue", parse=TRUE, size=3.5) +
          annotate("text", x=0, y=ymax*0.8, label = "paste(\"lme(\",delta,WS[es],\" ~ 1 + 1 | PATHID)\")", color="blue", parse=TRUE, size=3.5)

  ymax <- max(hist(model_S2S$S2S, breaks = seq(from=-15,to=15,by=0.1), plot=FALSE)$counts)
  p1_2 <- ggplot(model_S2S, aes(x=S2S)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=x_axis_range) +
          stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model_S2S$S2S))) * nrow(model_S2S) * 0.1, color = "red", size = 1) +
          annotate("text", x=-1, y=ymax*0.1, label = paste0("phi[S2S] ==",formatC(sqrt(var(model_S2S$S2S)),digits=2,format="f",flag="0")), color="red", parse=TRUE)

  ymax <- max(hist(model_P2P$P2P, breaks = seq(from=-15,to=15,by=0.1), plot=FALSE)$counts)
  p1_3 <- ggplot(model_P2P, aes(x=P2P)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=x_axis_range) +
          stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model_P2P$P2P))) * nrow(model_P2P) * 0.1, color = "red", size = 1) +
          annotate("text", x=-1, y=ymax*0.1, label = paste0("phi[P2P] ==",formatC(sqrt(var(model_P2P$P2P)),digits=2,format="f",flag="0")), color="red", parse=TRUE)


  ymax <- max(hist(model4_L2L$L2L, breaks = seq(from=-15,to=15,by=0.1), plot=FALSE)$counts)
  p4_1 <- ggplot(model4_L2L, aes(x=L2L)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=x_axis_range) +
    stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model4_L2L$L2L))) * nrow(model4_L2L) * 0.1, color = "red", size = 1) +
    annotate("text", x=-1, y=ymax*0.1, label = paste0("phi[L2L] ==",formatC(sqrt(var(model4_L2L$L2L)),digits=2,format="f",flag="0")), color="red", parse=TRUE) +
    annotate("text", x=0, y=ymax*1, label = "nlmer(obs ~ deriv(lnSa) ~ 1 | EQID", color="blue", size=3.5) +
    annotate("text", x=0, y=ymax*0.9, label = "+ 1 | STID + 1 | PATHID)", color="blue", size=3.5)
  ymax <- max(hist(model4_S2S$S2S, breaks = seq(from=-15,to=15,by=0.1), plot=FALSE)$counts)
  p4_2 <- ggplot(model4_S2S, aes(x=S2S)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=x_axis_range) +
    stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model4_S2S$S2S))) * nrow(model4_S2S) * 0.1, color = "red", size = 1) +
    annotate("text", x=-1, y=ymax*0.1, label = paste0("phi[S2S] ==",formatC(sqrt(var(model4_S2S$S2S)),digits=2,format="f",flag="0")), color="red", parse=TRUE)
  ymax <- max(hist(model4_P2P$P2P, breaks = seq(from=-15,to=15,by=0.1), plot=FALSE)$counts)
  p4_3 <- ggplot(model4_P2P, aes(x=P2P)) + geom_histogram(color="black", fill="white", binwidth = 0.1) + scale_x_continuous(lim=x_axis_range) +
    stat_function(fun = function(x) dnorm(x, mean = 0, sd = sqrt(var(model4_P2P$P2P))) * nrow(model4_P2P) * 0.1, color = "red", size = 1) +
    annotate("text", x=-1, y=ymax*0.1, label = paste0("phi[P2P] ==",formatC(sqrt(var(model4_P2P$P2P)),digits=2,format="f",flag="0")), color="red", parse=TRUE)
  g <- plot_grid(p1_1,p4_1,p1_2,p4_2,p1_3,p4_3, labels = "AUTO", nrow = 3, align = 'v')
  ggsave(file=paste(path_figures,"/Residual_P2P_",Dataname,".T", formatC(Period,digits=3,width=6,format="f",flag="0"),".png",sep=""),
             g,height=9, width=6, units='in', dpi=600)    
  toc()


  ############
  tic("L2L map")
  xc_list <- seq(from=lon_w, to=lon_e, by=0.1)
  yc_list <- seq(from=lat_s, to=lat_n, by=0.1)
      eqid <- rownames(model4_L2L)
      model4_L2L <- cbind(eqid, model4_L2L)
      lon <- Dataset[match(model4_L2L$eqid,Dataset$eqid),]$lon
      lat <- Dataset[match(model4_L2L$eqid,Dataset$eqid),]$lat
      map_data <- cbind(model4_L2L,lon,lat)
      ave_data <- data.frame(matrix(NA, nrow = length(xc_list)*length(yc_list), ncol = 4))
      colnames(ave_data) <- c("xx","yy","L2L","Ave")
      for (i in seq_along(xc_list)) {
        for (j in seq_along(yc_list)) {
          k = (i-1) * length(yc_list) + j
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
      p3 <- ggplot(data = ca_df, mapping = aes(x = long, y = lat)) + coord_fixed(xlim = c(lon_w, lon_e),  ylim = c(lat_s, lat_n), ratio=1.3) +
        geom_polygon(color = "black", fill = "gray")
      p3 <- p3 + geom_point(data=ave_data,aes(x=xx,y=yy),size=1, shape=21)
      p3 <- p3 + geom_point(data=ave_data[ave_data[,3]!="NaN",],aes(x=xx,y=yy,color=factor(Ave)),size=2) +
        scale_colour_manual(na.value="white", values = c("red","orange","yellow","green","cyan","blue","purple","darkblue"),guide=FALSE)

      ggsave(file=paste(path_figures,"/L2L_map_",Dataname,".T", formatC(Period,digits=3,width=6,format="f",flag="0"),".png",sep=""),
             p3,height=6, width=6, units='in', dpi=600)
      toc()
  ############

  ############
  	tic("P2P map") 
  	xc_list <- seq(from=lon_w, to=lon_e, by=0.1)
  	yc_list <- seq(from=lat_s, to=lat_n, by=0.1)
      ptid <- data.frame(as.integer(rownames(model4_P2P)))
      colnames(ptid) <- c("ptid")
      P2P_map <- cbind(model4_P2P, ptid)
      P2P_map <- P2P_map[P2P_map$ptid > 2980000 & substr(P2P_map$ptid,1,3)==298,]
      rgid <- data.frame(as.numeric(substr(P2P_map$ptid,4,7)))
      colnames(rgid) <- c("rgid")
      map_data <- cbind(P2P_map, rgid)
      ave_data <- data.frame(matrix(NA, nrow = length(xc_list)*length(yc_list), ncol = 4))
      colnames(ave_data) <- c("xx","yy","P2P","Ave")
      for (i in seq_along(xc_list)) {
        for (j in seq_along(yc_list)) {
          k = (i-1) * length(yc_list) + j
          xc <- xc_list[i]
          yc <- yc_list[j]
          ave_data[k,1] <- xc
          ave_data[k,2] <- yc
          if (length(map_data[map_data$rgid==k,]$P2P) == 1) {
            ave_data[k,3] <- map_data[map_data$rgid==k,]$P2P
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
      p1 <- ggplot(data = ca_df, mapping = aes(x = long, y = lat)) + coord_fixed(xlim = c(lon_w, lon_e),  ylim = c(lat_s, lat_n), ratio=1.3) +
        geom_polygon(color = "black", fill = "gray")
      p1 <- p1 + geom_point(data=ave_data,aes(x=xx,y=yy),size=1, shape=21)
      p1 <- p1 + geom_point(data=ave_data[ave_data[,3]!="NaN",],aes(x=xx,y=yy,color=factor(Ave)),size=2) +
        scale_colour_manual(na.value="white", values = c("red","orange","yellow","green","cyan","blue","purple","darkblue"),guide=FALSE)

      ggsave(file=paste(path_figures,"/P2P_map_",Dataname,".T", formatC(Period,digits=3,width=6,format="f",flag="0"),".png",sep=""),
             p1,height=6, width=6, units='in', dpi=600)
      toc()
  ############



  ############
  para_list <- c("rrup")
  label_list <- c("Rupture distance (km)")
  Intra_residual_linear(Dataname, Dataset, epsR, path_figures, Period, para_list, label_list)
  Path_term_linear(Dataname, Dataset, model4_P2P, path_figures, Period, para_list, label_list)
  para_list <- c("mag")
  label_list <- c("Magnitude")
  Intra_residual_linear(Dataname, Dataset, epsR, path_figures, Period, para_list, label_list)
  Source_term_linear(Dataname, Dataset, model4_L2L, path_figures, Period, para_list, label_list)
  para_list <- c("vs30")
  label_list <- c("Vs30 (m/s)")
  Intra_residual_linear(Dataname, Dataset, epsR, path_figures, Period, para_list, label_list)
  Site_term_linear(Dataname, Dataset, model4_S2S, path_figures, Period, para_list, label_list)
  ############
}
# End plot figures


if (output_flag) {
    Resi <- subset(Dataset,select = c(Resi,vs30,rrup,mag))
    Resi <- cbind(Resi,epsR)
    colnames(Resi) <- c("Resi","vs30","rrup","mag","aleatory")
    write.csv(Resi,file=paste(path_outputs,"/Residuals.",Dataname, ".T",
                                formatC(Period,digits=3,width=6,format="f",flag="0"),".csv", sep=""), row.names=TRUE)

    fit <- coef(lm(Resi ~ log10(vs30), data = Resi))
    write.csv(fit,file=paste0(path_outputs,"/Resi_vs30.csv"))
    fit <- coef(lm(aleatory ~ log10(vs30), data = Resi))
    write.csv(fit,file=paste0(path_outputs,"/Aleatory_vs30.csv"))
    nn <- seq(from=0.1,to=1,by=0.1)
    xmed <- rep(0,length(nn))
    ymed <- rep(0,length(nn))
    err1 <- rep(0,length(nn))
    err2 <- rep(0,length(nn))
    ymed2 <- rep(0,length(nn))
    for (i in seq_along(nn)) {
        x1 <- min(log10(Resi$vs30)) + (max(log10(Resi$vs30))-min(log10(Resi$vs30))) * (nn[i]-0.1)
        x2 <- min(log10(Resi$vs30)) + (max(log10(Resi$vs30))-min(log10(Resi$vs30))) * nn[i]
        temp <- Resi[Resi$vs30>10^x1 & Resi$vs30<=10^x2,]
        xmed[i] <- 10^((x1+x2)/2)
        ymed[i] <- median(temp$Resi)
        err1[i] <- sd(temp$Resi)
        ymed2[i] <- median(temp$aleatory)
        err2[i] <- sd(temp$aleatory)
    }
    err_data <- data.frame(cbind(xmed,ymed,err1,ymed2,err2))
    write.csv(err_data,file=paste0(path_outputs,"/Resi_vs30_err.csv"))

    fit <- coef(lm(Resi ~ log10(rrup), data = Resi))
    write.csv(fit,file=paste0(path_outputs,"/Resi_rrup.csv"))
    fit <- coef(lm(aleatory ~ log10(rrup), data = Resi))
    write.csv(fit,file=paste0(path_outputs,"/Aleatory_rrup.csv"))
    nn <- seq(from=0.1,to=1,by=0.1)
    xmed <- rep(0,length(nn))
    ymed <- rep(0,length(nn))
    err1 <- rep(0,length(nn))
    err2 <- rep(0,length(nn))
    ymed2 <- rep(0,length(nn))
    for (i in seq_along(nn)) {
        x1 <- 1 + (max(log10(Resi$rrup))-1) * (nn[i]-0.1)
        x2 <- 1 + (max(log10(Resi$rrup))-1) * nn[i]
        temp <- Resi[Resi$rrup>10^x1 & Resi$rrup<=10^x2,]
        xmed[i] <- 10^((x1+x2)/2)
        ymed[i] <- median(temp$Resi)
        err1[i] <- sd(temp$Resi)
	ymed2[i] <- median(temp$aleatory)
        err2[i] <- sd(temp$aleatory)
    }
    err_data <- data.frame(cbind(xmed,ymed,err1,ymed2,err2))
    write.csv(err_data,file=paste0(path_outputs,"/Resi_rrup_err.csv"))

    fit <- coef(lm(Resi ~ mag, data = Resi))
    write.csv(fit,file=paste0(path_outputs,"/Resi_mag.csv"))
    fit <- coef(lm(aleatory ~ mag, data = Resi))
    write.csv(fit,file=paste0(path_outputs,"/Aleatory_mag.csv"))
    nn <- seq(from=0.1,to=1,by=0.1)
    xmed <- rep(0,length(nn))
    ymed <- rep(0,length(nn))
    err1 <- rep(0,length(nn))
    err2 <- rep(0,length(nn))
    ymed2 <- rep(0,length(nn))
    for (i in seq_along(nn)) {
      x1 <- min(Resi$mag) + (max(Resi$mag)-min(Resi$mag)) * (nn[i]-0.1)
      x2 <- min(Resi$mag) + (max(Resi$mag)-min(Resi$mag)) * nn[i]
      temp <- Resi[Resi$mag>x1 & Resi$mag<=x2,]
      xmed[i] <- (x1+x2)/2
      ymed[i] <- median(temp$Resi)
      err1[i] <- sd(temp$Resi)
      ymed2[i] <- median(temp$aleatory)
      err2[i] <- sd(temp$aleatory)
    }
    err_data <- data.frame(cbind(xmed,ymed,err1,ymed2,err2))
    write.csv(err_data,file=paste0(path_outputs,"/Resi_mag_err.csv"))

    lon <- aggregate(Dataset[["hypo_lon"]],list(Dataset$eqid),mean)
    colnames(lon) <- c("eqid","hypo_lon")
    lat <- aggregate(Dataset[["hypo_lat"]],list(Dataset$eqid),mean)
    colnames(lat) <- c("eqid","hypo_lat")
    mag <- aggregate(Dataset[["mag"]],list(Dataset$eqid),mean)
    colnames(mag) <- c("eqid","mag")
    model_L2L <- cbind(model_L2L, lon, lat,mag)
    model_L2L <- subset(model_L2L,select=c(eqid,L2L,hypo_lon,hypo_lat,mag))
    write.csv(model_L2L,file=paste(path_outputs,"/Source_term.",Dataname, ".T",
                                    formatC(Period,digits=3,width=6,format="f",flag="0"),".csv", sep=""), row.names=TRUE)

    st_lon <- aggregate(Dataset[["st_lon"]],list(Dataset$stid),mean)
    colnames(st_lon) <- c("stid","st_lon")
    st_lat <- aggregate(Dataset[["st_lat"]],list(Dataset$stid),mean)
    colnames(st_lat) <- c("stid","st_lat")
    vs30 <- aggregate(Dataset[["vs30"]],list(Dataset$stid),mean)
    colnames(vs30) <- c("stid","vs30")
    model_S2S <- cbind(model_S2S, st_lon, st_lat, vs30)
    model_S2S <- subset(model_S2S,select=c(stid,S2S,st_lon,st_lat,vs30))
    write.csv(model_S2S,file=paste(path_outputs,"/Site_term.",Dataname, ".T",
                                    formatC(Period,digits=3,width=6,format="f",flag="0"),".csv", sep=""), row.names=TRUE)

    lon <- aggregate(Dataset[["lon"]],list(Dataset$pathid),mean)
    colnames(lon) <- c("pathid","lon")
    lat <- aggregate(Dataset[["lat"]],list(Dataset$pathid),mean)
    colnames(lat) <- c("pathid","lat")
    hypo_lon <- aggregate(Dataset[["hypo_lon"]],list(Dataset$pathid),mean)
    colnames(hypo_lon) <- c("eqid","hypo_lon")
    hypo_lat <- aggregate(Dataset[["hypo_lat"]],list(Dataset$pathid),mean)
    colnames(hypo_lat) <- c("eqid","hypo_lat")
    st_lon <- aggregate(Dataset[["st_lon"]],list(Dataset$pathid),mean)
    colnames(st_lon) <- c("pathid","st_lon")
    st_lat <- aggregate(Dataset[["st_lat"]],list(Dataset$pathid),mean)
    colnames(st_lat) <- c("pathid","st_lat")
    rrup <- aggregate(Dataset[["rrup"]],list(Dataset$pathid),mean)
    colnames(rrup) <- c("pathid","rrup")
    model_P2P <- cbind(model_P2P, st_lon, st_lat, lon, lat, hypo_lon, hypo_lat,rrup)
    model_P2P <- subset(model_P2P,select=c(pathid,P2P,lon,lat,hypo_lon,hypo_lat,st_lon,st_lat,rrup))
    write.csv(model_P2P,file=paste(path_outputs,"/Path_term.",Dataname, ".T",
                                    formatC(Period,digits=3,width=6,format="f",flag="0"),".csv", sep=""), row.names=TRUE)

    CC_100000 <- Resi[sample(nrow(Resi),100000),]
    write.csv(CC_100000,file=paste(path_outputs,"/Residuals_100000.",Dataname, ".T",
                                formatC(Period,digits=3,width=6,format="f",flag="0"),".csv", sep=""), row.names=TRUE)
 
} ### End output

