###################
### Intra event residual versus Rrup, Vs30, Sa1180 and Z1
Intra_residual <- function(Dataname, Dataset, Resi, epsR, path_figures, Period) {
  plot_data <- cbind(Dataset,Resi)
  plot_data <- plot_data[plot_data$Resi!="NaN",]
  plot_data <- cbind(plot_data,epsR)
  p1 <- ggplot(plot_data, aes(x = rrup, y = epsR)) + geom_point(shape=1) +
    scale_x_log10(name="Rupture distance (km)", lim=c(0.1,400), breaks=c(0.1,1,10,100,400),
                  label=c(0.1,1,10,100,400)) +
    scale_y_continuous(lim=c(-3,3),name="Epsilon") +
    theme(panel.grid.major = element_line(colour = "grey")) +
    geom_smooth(data=plot_data[plot_data$rrup>=5,], se = FALSE) +
    ggtitle(paste("T=",Period,"s"))
  
  p2 <- ggplot(plot_data, aes(x = vs30, y = epsR)) + geom_point(shape=1) +
    scale_x_log10(name="Vs30 (m/sec)", lim=c(100,3000), breaks=c(100,1000,3000),
                  label=c(100,1000,3000)) +
    scale_y_continuous(lim=c(-3,3),name="Epsilon") +
    theme(panel.grid.major = element_line(colour = "grey")) +
    geom_smooth(data=plot_data[plot_data$vs30>=100,], se = FALSE)
  
  p3 <- ggplot(plot_data[plot_data$vs30>180 & plot_data$vs30<=360,], aes(x = Sa1180, y = epsR)) +
    geom_point(shape=1) + scale_x_log10(name="Sa1180 (g)", lim=c(0.0001,2), breaks=c(0.0001,0.001,0.01,1),
                                        label=c(0.0001,0.001,0.01,1)) +
    scale_y_continuous(lim=c(-3,3),name="Epsilon") +
    theme(panel.grid.major = element_line(colour = "grey")) +
    geom_smooth(data=plot_data[plot_data$Sa1180>=0.001,], se = FALSE)
  
  p4 <- ggplot(plot_data, aes(x = z1_true, y = epsR)) + geom_point(shape=1) +
    scale_x_log10(name="Z1 (km)", lim=c(0.01,5),breaks=c(0.01,0.1,1,4),
                  label=c(0.01,0.1,1,4)) +
    scale_y_continuous(lim=c(-3,3),name="Epsilon") +
    theme(panel.grid.major = element_line(colour = "grey")) +
    geom_smooth(data=plot_data[plot_data$z1_true>=0.02,], se = FALSE)
  g <- plot_grid(p1,p2,p3,p4, labels = "AUTO", ncol = 1, align = 'v')
  
  ggsave(file=paste(path_figures,"/Intra_resids.",Dataname, ".T",
                    formatC(Period,digits=3,width=6,format="f",flag="0"),".png",sep=""),g,
         height=16, width=6, units='in', dpi=300)
}
###
###################



###################
### Intra event residual versus Rrup, Mag, Vs30 for GMPE VIL
Intra_residual_VIL <- function(Dataname, plot_data, epsR, path_figures, Period) {
  plot_data <- cbind(plot_data,epsR)
  
  p1 <- ggplot() + geom_point(data=plot_data, aes(x = rrup, y = epsR), shape=21, size=0.01)
  nn <- seq(from=-2,to=5,by=0.5)
  xmed <- rep(0,length(nn))
  ymed <- rep(0,length(nn))
  ystd <- rep(0,length(nn))
  for (i in seq_along(nn)) {
    x1 <- exp(nn[i])
    x2 <- exp(nn[i]+0.5)
    temp <- plot_data[plot_data$rrup>x1 & plot_data$rrup<=x2,]
    xmed[i] <- (x1+x2)/2
    ymed[i] <- median(temp$epsR)
    ystd[i] <- sd(temp$epsR)
  }
  err_data <- data.frame(cbind(xmed,ymed,ystd))
  p1 <- p1 + geom_point(data=err_data, aes(x = xmed, y = ymed), color="blue",shape=22) +
    geom_errorbar(data=err_data, aes(x=xmed, ymin=ymed-ystd, ymax=ymed+ystd), width=0.1, color="blue") +
    scale_x_log10(name="Rupture distance (km)", lim=c(0.1,200), breaks=c(0.1,1,10,100,200),
                  label=c(0.1,1,10,100,200)) +
    scale_y_continuous(lim=c(-3,3),name="Epsilon") +
    theme(panel.grid.major = element_line(colour = "grey")) +
    ggtitle(paste("T=",Period,"s")) 
  
  p2 <- ggplot() + geom_point(data=plot_data, aes(x = vs30, y = epsR), shape=21, size=0.01)
  nn <- seq(from=4,to=7,by=0.25)
  xmed <- rep(0,length(nn))
  ymed <- rep(0,length(nn))
  ystd <- rep(0,length(nn))
  for (i in seq_along(nn)) {
    x1 <- exp(nn[i])
    x2 <- exp(nn[i]+0.25)
    temp <- plot_data[plot_data$vs30>x1 & plot_data$vs30<=x2,]
    xmed[i] <- (x1+x2)/2
    ymed[i] <- median(temp$epsR)
    ystd[i] <- sd(temp$epsR)
  }
  err_data <- data.frame(cbind(xmed,ymed,ystd))
  p2 <- p2 + geom_point(data=err_data, aes(x = xmed, y = ymed), color="blue",shape=22) + 
    geom_errorbar(data=err_data, aes(x=xmed, ymin=ymed-ystd, ymax=ymed+ystd), width=0.1, color="blue") +
    scale_x_log10(name="Vs30 (m/sec)", lim=c(100,1000), breaks=c(100,200,300,400,500,600,700,800,900,1000),
                  label=c(100,200,300,400,500,600,700,800,900,1000)) +
    scale_y_continuous(lim=c(-3,3),name="Epsilon") +
    theme(panel.grid.major = element_line(colour = "grey")) 
  
  p3 <- ggplot() + geom_point(data=plot_data, aes(x = mag, y = epsR), shape=21, size=0.01) 
  nn <- seq(from=6.0,to=8.5,by=0.25)
  xmed <- rep(0,length(nn))
  ymed <- rep(0,length(nn))
  ystd <- rep(0,length(nn))
  for (i in seq_along(nn)) {
    x1 <- nn[i]
    x2 <- nn[i]+0.25
    temp <- plot_data[plot_data$mag>x1 & plot_data$mag<=x2,]
    xmed[i] <- (x1+x2)/2
    ymed[i] <- median(temp$epsR)
    ystd[i] <- sd(temp$epsR)
  }
  err_data <- data.frame(cbind(xmed,ymed,ystd))
  p3 <- p3 + geom_point(data=err_data, aes(x = xmed, y = ymed), color="blue",shape=22) + 
    geom_errorbar(aes(x=xmed, ymin=ymed-ystd, ymax=ymed+ystd), width=0.1, color="blue") +
    scale_x_continuous(name="Magnitude", lim=c(6,9),breaks=c(6,7,8,9),
                       label=c(6,7,8,9)) +
    scale_y_continuous(lim=c(-3,3),name="Epsilon") +
    theme(panel.grid.major = element_line(colour = "grey"))
  
  g <- plot_grid(p1,p2,p3, labels = "AUTO", ncol = 1, align = 'v')
  
  ggsave(file=paste(path_figures,"/Intra_resids.",Dataname, ".T",
                    formatC(Period,digits=3,width=6,format="f",flag="0"),".png",sep=""),g,
         height=9, width=6, units='in', dpi=600)
}
###
###################



###################
### Inter event residual versus Rrup, Vs30, Sa1180 and Z1
Inter_residual <- function(Dataname, Dataset, etaR, path_figures, Period) {
  plot_data <- data.frame(matrix(NA, nrow = nrow(etaR), ncol = 4))
  colnames(plot_data) <- c("mag","ztor","rake","etaR")
  for (k in 1:nrow(etaR)) {
    index <- as.numeric(rownames(etaR)[k])
    plot_data$mag[k] <- Dataset[Dataset$eqid==index,]$mag[1]
    plot_data$ztor[k] <- Dataset[Dataset$eqid==index,]$ztor[1]
    plot_data$rake[k] <- Dataset[Dataset$eqid==index,]$rake[1]
    plot_data$etaR[k] <- etaR[k]
  }
  p1 <- ggplot(plot_data, aes(x = mag, y = etaR)) + geom_point(shape=1) +
    scale_x_continuous(name="Magnitude", lim=c(3,8)) +
    scale_y_continuous(name="Eta",lim=c(-2,2)) +
    theme(panel.grid.major = element_line(colour = "grey")) +
    ggtitle(paste("T=",Period,"s"))
  p2 <- ggplot(plot_data, aes(x = ztor, y = etaR)) + geom_point(shape=1) +
    scale_x_continuous(name="Depth to top of rupture (km)", lim=c(0,30)) +
    scale_y_continuous(name="Eta",lim=c(-2,2)) +
    theme(panel.grid.major = element_line(colour = "grey"))
  p3 <- ggplot(plot_data, aes(x = rake, y = etaR)) + geom_point(shape=1) +
    scale_x_continuous(name="Rake", lim=c(-90,90)) +
    scale_y_continuous(name="Eta",lim=c(-2,2)) +
    theme(panel.grid.major = element_line(colour = "grey"))
  g <- plot_grid(p1,p2,p3, labels = "AUTO", ncol = 1, align = 'v')
  ggsave(file=paste(path_figures,"/Inter_resids.",Dataname, ".T",
                    formatC(Period,digits=3,width=6,format="f",flag="0"),".png",sep=""),g,
         height=10, width=8, units='in', dpi=300)
}
###
###################



###################
### Inter event residual versus Rrup, Vs30, Magnitude
Inter_residual_VIL <- function(Dataname, Dataset, model4_L2L, model4_S2S, model4_P2P, path_figures, Period) {
  
  plot_data <- aggregate(Dataset$mag,list(Dataset$eqid),mean)
  colnames(plot_data) <- c("eqid","mag")
  plot_data <- cbind(plot_data,model4_L2L)
  p1 <- ggplot() + geom_point(data=plot_data, aes(x = mag, y = L2L), shape=21, size=0.01) 
  nn <- seq(from=6.0,to=8.5,by=0.25)
  xmed <- rep(0,length(nn))
  ymed <- rep(0,length(nn))
  ystd <- rep(0,length(nn))
  for (i in seq_along(nn)) {
    x1 <- nn[i]
    x2 <- nn[i]+0.25
    temp <- plot_data[plot_data$mag>x1 & plot_data$mag<=x2,]
    xmed[i] <- (x1+x2)/2
    ymed[i] <- median(temp$L2L)
    ystd[i] <- sd(temp$L2L)
  }
  err_data <- data.frame(cbind(xmed,ymed,ystd))
  p1 <- p1 + geom_point(data=err_data, aes(x = xmed, y = ymed), color="blue",shape=22) + 
    geom_errorbar(data=err_data, aes(x=xmed, ymin=ymed-ystd, ymax=ymed+ystd), width=0.1, color="blue") +
    scale_x_continuous(name="Magnitude", lim=c(6,8)) +
    scale_y_continuous(name="Source term",lim=c(-2,2)) +
    theme(panel.grid.major = element_line(colour = "grey")) +
    ggtitle(paste("T=",Period,"s"))
  
  
  plot_data <- aggregate(Dataset$vs30,list(Dataset$stid),mean)
  colnames(plot_data) <- c("stid","vs30")
  plot_data <- cbind(plot_data,model4_S2S)
  p2 <- ggplot() + geom_point(data=plot_data, aes(x = vs30, y = S2S), shape=21, size=0.01) 
  nn <- seq(from=4,to=7,by=0.25)
  xmed <- rep(0,length(nn))
  ymed <- rep(0,length(nn))
  ystd <- rep(0,length(nn))
  for (i in seq_along(nn)) {
    x1 <- exp(nn[i])
    x2 <- exp(nn[i]+0.25)
    temp <- plot_data[plot_data$vs30>x1 & plot_data$vs30<=x2,]
    xmed[i] <- (x1+x2)/2
    ymed[i] <- median(temp$S2S)
    ystd[i] <- sd(temp$S2S)
  }
  err_data <- data.frame(cbind(xmed,ymed,ystd))
  p2 <- p2 + geom_point(data=err_data, aes(x = xmed, y = ymed), color="blue",shape=22) + 
    geom_errorbar(data=err_data, aes(x=xmed, ymin=ymed-ystd, ymax=ymed+ystd), width=0.1, color="blue") +
    scale_x_log10(name="Vs30 (m/sec)", lim=c(100,1000), breaks=c(100,200,300,400,500,1000),
                  label=c(100,200,300,400,500,1000)) +
    scale_y_continuous(name="Site term",lim=c(-2,2)) +
    theme(panel.grid.major = element_line(colour = "grey"))
  
  
  plot_data <- aggregate(Dataset$rrup,list(Dataset$pathid),mean)
  colnames(plot_data) <- c("pathid","rrup")
  plot_data <- cbind(plot_data,model4_P2P)
  p3 <- ggplot() + geom_point(data=plot_data, aes(x = rrup, y = P2P), shape=21, size=0.01) 
  nn <- seq(from=-2,to=5,by=0.5)
  xmed <- rep(0,length(nn))
  ymed <- rep(0,length(nn))
  ystd <- rep(0,length(nn))
  for (i in seq_along(nn)) {
    x1 <- exp(nn[i])
    x2 <- exp(nn[i]+0.5)
    temp <- plot_data[plot_data$rrup>x1 & plot_data$rrup<=x2,]
    xmed[i] <- (x1+x2)/2
    ymed[i] <- median(temp$P2P)
    ystd[i] <- sd(temp$P2P)
  }
  err_data <- data.frame(cbind(xmed,ymed,ystd))
  p3 <- p3 + geom_point(data=err_data, aes(x = xmed, y = ymed), color="blue",shape=22) + 
    geom_errorbar(data=err_data, aes(x=xmed, ymin=ymed-ystd, ymax=ymed+ystd), width=0.1, color="blue") +
    scale_x_log10(name="Rupture distance (km)", lim=c(0.1,200), breaks=c(0.1,1,10,100,200),
                  label=c(0.1,1,10,100,200)) +
    scale_y_continuous(name="Path term",lim=c(-2,2)) +
    theme(panel.grid.major = element_line(colour = "grey"))
  
  
  g <- plot_grid(p1,p2,p3, labels = "AUTO", ncol = 1, align = 'v')
  ggsave(file=paste(path_figures,"/Inter_resids.",Dataname, ".T",
                    formatC(Period,digits=3,width=6,format="f",flag="0"),".png",sep=""),g,
         height=9, width=6, units='in', dpi=600)
}
###
###################



###################
### Residuals plots comparing to Villani slides
Compare_residual_Villani <- function(Dataname, Dataset, model4_L2L, epsR, path_figures, Period) {
    p1 <- ggplot() + geom_point(data=Dataset, aes(x = mag, y = Resi), shape=21, size=0.01) + 
      scale_x_continuous(name="Magnitude", lim=c(6,9),breaks=c(6,7,8,9), label=c(6,7,8,9)) +
      scale_y_continuous(lim=c(-4,4),name="residual")
    miu <- mean(Dataset$Resi)
    sigma <- sd(Dataset$Resi)
    plot_data2 <- data.frame(x1=c(6,9),y1=rep(miu,2),y2=rep(miu + sigma, 2),y3=rep(miu-sigma, 2))
    p1 <- p1 + geom_line(data=plot_data2, aes(x=x1,y=y1), color="red") + 
      geom_line(data=plot_data2, aes(x=x1,y=y2), linetype="dotted") + 
      geom_line(data=plot_data2, aes(x=x1,y=y3), linetype="dotted") +
      annotate("text", x=6.5, y=3, label = paste0("miu ==",formatC(miu,digits=2,format="f",flag="0")), parse=TRUE) +
      annotate("text", x=7, y=3, label = paste0("sigma ==",formatC(sigma,digits=2,format="f",flag="0")), parse=TRUE)
                                            
    
    plot_data <- aggregate(Dataset$mag,list(Dataset$eqid),mean)
    colnames(plot_data) <- c("eqid","mag")
    plot_data <- cbind(plot_data,model4_L2L)
    p2 <- ggplot() + geom_point(data=plot_data, aes(x = mag, y = L2L), shape=21, size=0.01) +
      scale_x_continuous(name="Magnitude", lim=c(6,9),breaks=c(6,7,8,9), label=c(6,7,8,9)) +
      scale_y_continuous(lim=c(-2,2),name="L2L")
    miu <- mean(plot_data$L2L)
    sigma <- sd(plot_data$L2L)
    plot_data2 <- data.frame(x1=c(6,9),y1=rep(miu,2),y2=rep(miu + sigma, 2),y3=rep(miu-sigma, 2))
    p2 <- p2 + geom_line(data=plot_data2, aes(x=x1,y=y1), color="red") + 
      geom_line(data=plot_data2, aes(x=x1,y=y2), linetype="dotted") + 
      geom_line(data=plot_data2, aes(x=x1,y=y3), linetype="dotted") +
      annotate("text", x=6.5, y=1, label = paste0("miu ==",formatC(miu,digits=2,format="f",flag="0")), parse=TRUE) +
      annotate("text", x=7, y=1, label = paste0("sigma ==",formatC(sigma,digits=2,format="f",flag="0")), parse=TRUE)
    
    plot_data <- data.frame(cbind(Dataset$mag,epsR))
    colnames(plot_data) <- c("mag","Wes")
    p3 <- ggplot() + geom_point(data=plot_data, aes(x = mag, y = Wes), shape=21, size=0.01) +
      scale_x_continuous(name="Magnitude", lim=c(6,9),breaks=c(6,7,8,9), label=c(6,7,8,9)) +
      scale_y_continuous(lim=c(-4,4),name="Wes")
    miu <- mean(plot_data$Wes)
    sigma <- sd(plot_data$Wes)
    plot_data2 <- data.frame(x1=c(6,9),y1=rep(miu,2),y2=rep(miu + sigma, 2),y3=rep(miu-sigma, 2))
    p3 <- p3 + geom_line(data=plot_data2, aes(x=x1,y=y1), color="red") + 
      geom_line(data=plot_data2, aes(x=x1,y=y2), linetype="dotted") + 
      geom_line(data=plot_data2, aes(x=x1,y=y3), linetype="dotted") +
      annotate("text", x=6.5, y=3, label = paste0("miu ==",formatC(miu,digits=2,format="f",flag="0")), parse=TRUE) +
      annotate("text", x=7, y=3, label = paste0("sigma ==",formatC(sigma,digits=2,format="f",flag="0")), parse=TRUE)
    
    g <- plot_grid(p1,p2,p3, labels = "AUTO", ncol = 1, align = 'v')
    ggsave(file=paste(path_figures,"/Compare_Villani_mag.",Dataname, ".T",
                      formatC(Period,digits=3,width=6,format="f",flag="0"),".png",sep=""),g,
                      height=9, width=6, units='in', dpi=600)
    
    
    p1 <- ggplot() + geom_point(data=Dataset, aes(x = rrup, y = Resi), shape=21, size=0.01) + 
      scale_x_continuous(name="Rrup", lim=c(0,100),breaks=c(0,20,40,60,80,100), label=c(0,20,40,60,80,100)) +
      scale_y_continuous(lim=c(-4,4),name="residual")
    miu <- mean(Dataset$Resi)
    sigma <- sd(Dataset$Resi)
    plot_data2 <- data.frame(x1=c(0,100),y1=rep(miu,2),y2=rep(miu + sigma, 2),y3=rep(miu-sigma, 2))
    p1 <- p1 + geom_line(data=plot_data2, aes(x=x1,y=y1), color="red") + 
      geom_line(data=plot_data2, aes(x=x1,y=y2), linetype="dotted") + 
      geom_line(data=plot_data2, aes(x=x1,y=y3), linetype="dotted") +
      annotate("text", x=5.5, y=3, label = paste0("miu ==",formatC(miu,digits=2,format="f",flag="0")), parse=TRUE) +
      annotate("text", x=20, y=3, label = paste0("sigma ==",formatC(sigma,digits=2,format="f",flag="0")), parse=TRUE)
    
    plot_data <- aggregate(Dataset$mag,list(Dataset$eqid),mean)
    colnames(plot_data) <- c("eqid","mag")
    plot_data <- cbind(plot_data,model4_L2L)
    Dataset$L2L <- plot_data[match(Dataset$eqid,plot_data$eqid),3]
    p2 <- ggplot() + geom_point(data=Dataset, aes(x = rrup, y = L2L), shape=21, size=0.01) +
      scale_x_continuous(name="Rrup", lim=c(0,100),breaks=c(0,20,40,60,80,100), label=c(0,20,40,60,80,100)) +
      scale_y_continuous(lim=c(-2,2),name="L2L")
    miu <- mean(Dataset$L2L)
    sigma <- sd(Dataset$L2L)
    plot_data2 <- data.frame(x1=c(0,100),y1=rep(miu,2),y2=rep(miu + sigma, 2),y3=rep(miu-sigma, 2))
    p2 <- p2 + geom_line(data=plot_data2, aes(x=x1,y=y1), color="red") + 
      geom_line(data=plot_data2, aes(x=x1,y=y2), linetype="dotted") + 
      geom_line(data=plot_data2, aes(x=x1,y=y3), linetype="dotted") +
      annotate("text", x=5.5, y=1, label = paste0("miu ==",formatC(miu,digits=2,format="f",flag="0")), parse=TRUE) +
      annotate("text", x=20, y=1, label = paste0("sigma ==",formatC(sigma,digits=2,format="f",flag="0")), parse=TRUE)
    
    plot_data <- data.frame(cbind(Dataset$rrup,epsR))
    colnames(plot_data) <- c("rrup","Wes")
    p3 <- ggplot() + geom_point(data=plot_data, aes(x = rrup, y = Wes), shape=21, size=0.01) +
      scale_x_continuous(name="Rrup", lim=c(0,100),breaks=c(0,20,40,60,80,100), label=c(0,20,40,60,80,100))+
      scale_y_continuous(lim=c(-4,4),name="Wes")
    miu <- mean(plot_data$Wes)
    sigma <- sd(plot_data$Wes)
    plot_data2 <- data.frame(x1=c(0,100),y1=rep(miu,2),y2=rep(miu + sigma, 2),y3=rep(miu-sigma, 2))
    p3 <- p3 + geom_line(data=plot_data2, aes(x=x1,y=y1), color="red") + 
      geom_line(data=plot_data2, aes(x=x1,y=y2), linetype="dotted") + 
      geom_line(data=plot_data2, aes(x=x1,y=y3), linetype="dotted") +
      annotate("text", x=5.5, y=3, label = paste0("miu ==",formatC(miu,digits=2,format="f",flag="0")), parse=TRUE) +
      annotate("text", x=20, y=3, label = paste0("sigma ==",formatC(sigma,digits=2,format="f",flag="0")), parse=TRUE)
    
    g <- plot_grid(p1,p2,p3, labels = "AUTO", ncol = 1, align = 'v')
    ggsave(file=paste(path_figures,"/Compare_Villani_rrup.",Dataname, ".T",
                      formatC(Period,digits=3,width=6,format="f",flag="0"),".png",sep=""),g,
           height=9, width=6, units='in', dpi=600)
}



###################
### Plot response spectrums
Response_spectrum <- function(Dataname,path_figures,plot_data,plot_data2,vs30,vs30_2,rjb) {
    p1 <- ggplot(plot_data, aes(x = Period, y = M5)) + geom_line(color = "blue") +
      scale_x_log10(name="Period (s)", lim=c(0.001,10), breaks=c(0.001,0.01,0.1,1,10),
                    label=c(0.001,0.01,0.1,1,10)) +
      scale_y_log10(lim=c(0.0001,1),name="lnSa", breaks=c(0.0001,0.001,0.01,0.1,1),
                    label=c(0.0001,0.001,0.01,0.1,1)) +
      annotation_logticks(base = 10) +
      theme(panel.grid.major = element_line(colour = "grey")) +
      ggtitle(paste0("Vs30=",vs30,"m/s   ","Rjb=",rjb,"km")) +
      geom_line(data=plot_data, aes(x = Period, y = M6),color = "cyan") +
      geom_line(data=plot_data, aes(x = Period, y = M7),color = "green") +
      geom_line(data=plot_data, aes(x = Period, y = M8),color = "red") +
      annotate("text", x=0.002, y=0.0001, label = "M=5", color="blue") +
      annotate("text", x=0.002, y=0.0002, label = "M=6", color="cyan") +
      annotate("text", x=0.002, y=0.0005, label = "M=7", color="green") +
      annotate("text", x=0.002, y=0.001, label = "M=8", color="red")
    
    p2 <- ggplot(plot_data2, aes(x = Period, y = M5)) + geom_line(color = "blue") +
      scale_x_log10(name="Period (s)", lim=c(0.001,10), breaks=c(0.001,0.01,0.1,1,10),
                    label=c(0.001,0.01,0.1,1,10)) +
      scale_y_log10(lim=c(0.0001,1),name="lnSa", breaks=c(0.0001,0.001,0.01,0.1,1),
                    label=c(0.0001,0.001,0.01,0.1,1)) +
      annotation_logticks(base = 10) +
      theme(panel.grid.major = element_line(colour = "grey")) +
      ggtitle(paste0("Vs30=",vs30_2,"m/s   ","Rjb=",rjb,"km")) +
      geom_line(data=plot_data2, aes(x = Period, y = M6),color = "cyan") +
      geom_line(data=plot_data2, aes(x = Period, y = M7),color = "green") +
      geom_line(data=plot_data2, aes(x = Period, y = M8),color = "red") +
      annotate("text", x=0.002, y=0.0001, label = "M=5", color="blue") +
      annotate("text", x=0.002, y=0.0002, label = "M=6", color="cyan") +
      annotate("text", x=0.002, y=0.0005, label = "M=7", color="green") +
      annotate("text", x=0.002, y=0.001, label = "M=8", color="red")
    
    g <- plot_grid(p1,p2, labels = "AUTO", ncol = 2)
    ggsave(file=paste(path_figures,"/Response_spectrum.",Dataname,".Rjb=",rjb,".png",sep=""),
           g,height=8, width=12, units='in', dpi=300)
}
###
###################




###################
### Plot magnitude and distance scaling
Mag_Dist_Scaling <- function(Dataname,path_figures,plot_data,plot_data2,vs30,desired_period) {
    p1 <- ggplot(plot_data, aes(x = Distance, y = M5)) + geom_line(color = "blue") +
      scale_x_log10(name="Distance (km)", lim=c(1,400), breaks=c(1,10,100,200),
                    label=c(1,10,100,200)) +
      scale_y_log10(lim=c(0.00001,1),name="lnSa", breaks=c(0.001,0.01,0.1,1,10),
                    label=c(0.001,0.01,0.1,1,10)) +
      annotation_logticks(base = 10) +
      theme(panel.grid.major = element_line(colour = "grey")) +
      ggtitle(paste0("Vs30=",vs30,"m/s   ","T=",desired_period,"s")) + 
      geom_line(data=plot_data, aes(x = Distance, y = M6),color = "cyan") + 
      geom_line(data=plot_data, aes(x = Distance, y = M7),color = "green") +
      geom_line(data=plot_data, aes(x = Distance, y = M8),color = "red") + 
      annotate("text", x=2, y=0.001, label = "M=5", color="blue") + 
      annotate("text", x=2, y=0.002, label = "M=6", color="cyan") + 
      annotate("text", x=2, y=0.005, label = "M=7", color="green") + 
      annotate("text", x=2, y=0.01, label = "M=8", color="red")
    
    p2 <- ggplot(plot_data2, aes(x = Magnitude, y = T1)) + geom_line(color = "blue") +
      scale_x_continuous(name="Magnitude", lim=c(4,9), breaks=c(3,4,5,6,7,8),
                    label=c(3,4,5,6,7,8)) +
      scale_y_log10(lim=c(0.00001,1),name="lnSa", breaks=c(0.0001,0.001,0.01,0.1,1,10),
                    label=c(0.0001,0.001,0.01,0.1,1,10)) +
      theme(panel.grid.major = element_line(colour = "grey")) +
      geom_line(data=plot_data2, aes(x = Magnitude, y = T30),color = "green") + 
      geom_line(data=plot_data2, aes(x = Magnitude, y = T200),color = "red") +
      annotate("text", x=3.5, y=1, label = "Rjb=1", color="blue") + 
      annotate("text", x=3.5, y=2, label = "Rjb=30", color="green") + 
      annotate("text", x=3.5, y=5, label = "Rjb=200", color="red")
    g <- plot_grid(p1,p2, labels = "AUTO", ncol = 1, align = 'v')
    ggsave(file=paste(path_figures,"/Distance_magnitude_scaling.",Dataname,"T",desired_period,"s.png",sep=""),
           g,height=12, width=6, units='in', dpi=300)
}
###
###################


###################
### Compare magnitude and distance scaling with ASK14
Compare_Mag_Dist_Scaling <- function(Dataname,path_figures,plot_data,plot_data2,
                                     vs30,desired_period,ASK14_data,ASK14_data2) {
  p1 <- ggplot(plot_data, aes(x = Distance, y = M5)) + geom_line(color = "blue") +
    scale_x_log10(name="Distance (km)", lim=c(1,400), breaks=c(1,10,100,200), label=c(1,10,100,200)) +
    scale_y_log10(lim=c(0.00001,1),name="lnSa", breaks=c(0.001,0.01,0.1,1,10), label=c(0.001,0.01,0.1,1,10)) +
    annotation_logticks(base = 10) +
    theme(panel.grid.major = element_line(colour = "grey")) +
    ggtitle(paste0("Vs30=",vs30," T=",desired_period," Rjb=30"," ZTOR=",ztor," W=",w," Dip=",dip," Z1=",z1)) + 
    geom_line(data=plot_data, aes(x = Distance, y = M6),color = "cyan") + 
    geom_line(data=plot_data, aes(x = Distance, y = M7),color = "green") +
    geom_line(data=plot_data, aes(x = Distance, y = M8),color = "red") + 
    geom_line(data=ASK14_data, aes(x = Distance, y = M5),color = "blue",linetype=2) + 
    geom_line(data=ASK14_data, aes(x = Distance, y = M6),color = "cyan",linetype=2) + 
    geom_line(data=ASK14_data, aes(x = Distance, y = M7),color = "green",linetype=2) +
    geom_line(data=ASK14_data, aes(x = Distance, y = M8),color = "red",linetype=2) + 
    annotate("text", x=2, y=0.00002, label = "M=5", color="blue", size = 6) + 
    annotate("text", x=2, y=0.0002, label = "M=6", color="cyan", size = 6) + 
    annotate("text", x=2, y=0.002, label = "M=7", color="green", size = 6) + 
    annotate("text", x=2, y=0.02, label = "M=8", color="red", size = 6)
  
 
  compare_data <- cbind(plot_data, ASK14_data)
  colnames(compare_data) <- c("Distance","M5","M6","M7","M8","Distance_ASK","M5_ASK","M6_ASK","M7_ASK","M8_ASK")
  p2 <- ggplot(compare_data, aes(x=M5_ASK, y=M5)) + geom_point(shape=1,size=1,color="blue") +
        scale_y_log10(lim=c(0.00001,1),name="Prediction from R code", breaks=c(0.001,0.01,0.1,1,10), label=c(0.001,0.01,0.1,1,10)) +
        scale_x_log10(lim=c(0.00001,1),name="From OpenSHA", breaks=c(0.001,0.01,0.1,1,10), label=c(0.001,0.01,0.1,1,10)) +
        annotation_logticks(base = 10) +
        theme(panel.grid.major = element_line(colour = "grey"))
  p2 <- p2 + geom_point(data=compare_data, aes(x=M6_ASK, y=M6),shape=1,size=1,color="cyan") +
             geom_point(data=compare_data, aes(x=M7_ASK, y=M7),shape=1,size=1,color="green") +
             geom_point(data=compare_data, aes(x=M8_ASK, y=M8),shape=1,size=1,color="red") 
    
  
  p3 <- ggplot(plot_data2, aes(x = Magnitude, y = T1)) + geom_line(color = "blue") +
              scale_x_log10(name="Magnitude", lim=c(4,9), breaks=c(3,4,5,6,7,8),
                  label=c(3,4,5,6,7,8)) +
              scale_y_log10(lim=c(0.00001,1),name="lnSa", breaks=c(0.0001,0.001,0.01,0.1,1,10),
                  label=c(0.0001,0.001,0.01,0.1,1,10)) +
              annotation_logticks(base = 10) +
              theme(panel.grid.major = element_line(colour = "grey")) +
              ggtitle(paste0("Vs30=",vs30," T=",desired_period," ZTOR=",ztor," W=",w," Dip=",dip," Z1=",z1)) + 
              geom_line(data=plot_data2, aes(x = Magnitude, y = T30),color = "green") + 
              geom_line(data=plot_data2, aes(x = Magnitude, y = T200),color = "red") +
              geom_line(data=ASK14_data2, aes(x = Magnitude, y = T1),color = "blue",linetype=2) +
              geom_line(data=ASK14_data2, aes(x = Magnitude, y = T30),color = "green",linetype=2) +
              geom_line(data=ASK14_data2, aes(x = Magnitude, y = T200),color = "red",linetype=2) + 
              annotate("text", x=4.5, y=0.00002, label = "R[jb]==1", color="blue", size = 6, parse=TRUE) + 
              annotate("text", x=4.5, y=0.0002, label = "R[jb]==30", color="green", size = 6, parse=TRUE) + 
              annotate("text", x=4.5, y=0.002, label = "R[jb]==200", color="red", size = 6, parse=TRUE)
  
  
  compare_data <- cbind(plot_data2, ASK14_data2)
  colnames(compare_data) <- c("Magnitude","T1","T30","T200","Magnitude_ASK","T1_ASK","T30_ASK","T200_ASK")
  p4 <- ggplot(compare_data, aes(x=T1_ASK, y=T1)) + geom_point(shape=1,size=1,color="blue") +
        scale_y_log10(lim=c(0.00001,1),name="Prediction from R code", breaks=c(0.001,0.01,0.1,1,10), label=c(0.001,0.01,0.1,1,10)) +
        scale_x_log10(lim=c(0.00001,1),name="From OpenSHA", breaks=c(0.001,0.01,0.1,1,10), label=c(0.001,0.01,0.1,1,10)) +
        annotation_logticks(base = 10) +
        theme(panel.grid.major = element_line(colour = "grey"))
  p4 <- p4 + geom_point(data=compare_data, aes(x=T30_ASK, y=T30),shape=1,size=1,color="green") +
        geom_point(data=compare_data, aes(x=T200_ASK, y=T200),shape=1,size=1,color="red")
  
  
  g <- plot_grid(p1,p2,p3,p4, labels = "AUTO", ncol = 2, nrow = 2)
  ggsave(file=paste(path_figures,"/Compare_Distance_magnitude_scaling.",Dataname,"T",desired_period,"s.png",sep=""),
         g,height=12, width=12, units='in', dpi=300)
}
###
###################



###################
###
Compare_response_Villani <- function(vs30) {
  lnSa <- function(mag, dist) {
    return(b1 + b2 * mag + b3 * mag ^ 2 + (b4 + b5 * mag) * 
             log(sqrt(4.0 ^ 2 + dist ^ 2)) + b6 * log(vs30/500))
  }
  
  b1 <- -15.8
  b2 <- 5.36
  b3 <- -0.33
  b4 <- -2.05
  b5 <- 0.20
  b6 <- 1.11
  
  dist <- seq(0,200,1)
  mag <- rep(6,length(dist))
  Dataset <- cbind(dist, mag)
  colnames(Dataset) <- c("dist", "mag")
  sa <- exp(apply(Dataset[,c("dist", "mag")], 1, function(Dataset) lnSa(Dataset['mag'],Dataset['dist'])))/100
  Dataset <- data.frame(cbind(Dataset,sa))
  p1 <- ggplot() + geom_line(data=Dataset,aes(x=dist,y=sa),color="red") + 
    scale_x_log10(name="Distance (km)", lim=c(1,200), breaks=c(1,10,100,200), label=c(1,10,100,200)) +
    scale_y_log10(lim=c(0.001,10),name="Sa (g)", breaks=c(0.001,0.01,0.1,1,10), label=c(0.001,0.01,0.1,1,10)) +
    annotation_logticks(base = 10) + theme(panel.grid.major = element_line(colour = "grey"))
  
  dist <- seq(0,200,1)
  mag <- rep(7,length(dist))
  Dataset <- cbind(dist, mag)
  colnames(Dataset) <- c("dist", "mag")
  sa <- exp(apply(Dataset[,c("dist", "mag")], 1, function(Dataset) lnSa(Dataset['mag'],Dataset['dist'])))/100
  Dataset <- data.frame(cbind(Dataset,sa))
  p1 <- p1 + geom_line(data=Dataset,aes(x=dist,y=sa),color="green")
  
  dist <- seq(0,200,1)
  mag <- rep(8,length(dist))
  Dataset <- cbind(dist, mag)
  colnames(Dataset) <- c("dist", "mag")
  sa <- exp(apply(Dataset[,c("dist", "mag")], 1, function(Dataset) lnSa(Dataset['mag'],Dataset['dist'])))/100
  Dataset <- data.frame(cbind(Dataset,sa))
  p1 <- p1 + geom_line(data=Dataset,aes(x=dist,y=sa),color="blue")
  
  b1 <- 0.947
  b2 <- 2.148
  b3 <- -0.159
  b4 <- -4.286
  b5 <- 0.435
  b6 <- -1.005
  
  dist <- seq(0,200,1)
  mag <- rep(6,length(dist))
  Dataset <- cbind(dist, mag)
  colnames(Dataset) <- c("dist", "mag")
  sa <- exp(apply(Dataset[,c("dist", "mag")], 1, function(Dataset) lnSa(Dataset['mag'],Dataset['dist'])))/100
  Dataset <- data.frame(cbind(Dataset,sa))
  p1 <- p1 + geom_line(data=Dataset,aes(x=dist,y=sa),color="red",linetype="dotted")
  
  dist <- seq(0,200,1)
  mag <- rep(7,length(dist))
  Dataset <- cbind(dist, mag)
  colnames(Dataset) <- c("dist", "mag")
  sa <- exp(apply(Dataset[,c("dist", "mag")], 1, function(Dataset) lnSa(Dataset['mag'],Dataset['dist'])))/100
  Dataset <- data.frame(cbind(Dataset,sa))
  p1 <- p1 + geom_line(data=Dataset,aes(x=dist,y=sa),color="green",linetype="dotted")
  
  dist <- seq(0,200,1)
  mag <- rep(8,length(dist))
  Dataset <- cbind(dist, mag)
  colnames(Dataset) <- c("dist", "mag")
  sa <- exp(apply(Dataset[,c("dist", "mag")], 1, function(Dataset) lnSa(Dataset['mag'],Dataset['dist'])))/100
  Dataset <- data.frame(cbind(Dataset,sa))
  p1 <- p1 + geom_line(data=Dataset,aes(x=dist,y=sa),color="blue",linetype="dotted") +
    annotate("text", x=4.5, y=0.002, label = "M=6", color="red", size = 6) + 
    annotate("text", x=4.5, y=0.005, label = "M=7", color="green", size = 6) + 
    annotate("text", x=4.5, y=0.02, label = "M=8", color="blue", size = 6)
  
}