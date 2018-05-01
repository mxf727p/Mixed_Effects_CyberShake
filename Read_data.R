Read_data <- function(dataname, filename) {
  if (dataname == "NGA_W2") {
    # Load data file
    dataset <- read.csv(filename)
    dataset <- dataset[dataset$EQID!=-999 & dataset$Earthquake.Magnitude!=-999 & dataset$Dip..deg.!=-999 & 
                         dataset$Rake.Angle..deg.!=-999 & dataset$Depth.to.Top.Of.Fault.Rupture.Model!=-999 &
                         dataset$ClstD..km.!=-999 & dataset$Joyner.Boore.Dist...km.!=-999 & 
                         dataset$Rx!=-999 & dataset$Vs30..m.s..selected.for.analysis!=-999,]
    # Remove certain events
    dataset<-dataset[dataset$EQID!=93,]
    dataset<-dataset[dataset$EQID!=142,]
    dataset<-dataset[dataset$EQID!=109,]
    dataset<-dataset[dataset$EQID!=3,]
    dataset<-dataset[dataset$EQID!=7,]
    dataset<-dataset[dataset$EQID!=22,]
    dataset<-dataset[dataset$EQID!=26,]
    dataset<-dataset[dataset$EQID!=67,]
    dataset<-dataset[dataset$EQID!=1009,]
    dataset<-dataset[dataset$EQID!=1256,]
    dataset<-dataset[dataset$EQID!=1244,]
    dataset<-dataset[dataset$EQID!=242,]
    dataset<-dataset[dataset$EQID!=1136,]
    dataset<-dataset[dataset$EQID!=154,]
    dataset<-dataset[dataset$EQID!=214,]
    dataset<-dataset[dataset$EQID!=222,]
    dataset<-dataset[dataset$EQID!=250,]
    dataset<-dataset[dataset$EQID!=259,]
    dataset<-dataset[dataset$EQID!=203,]
    dataset<-dataset[dataset$EQID!=71,]
    dataset<-dataset[dataset$EQID!=86,]
    dataset<-dataset[dataset$EQID!=95,]
    dataset<-dataset[dataset$EQID!=100,]
    dataset<-dataset[dataset$EQID<282 | dataset$EQID>345,]
    dataset<-dataset[dataset$GMX.s.C1!="C",]
    dataset<-dataset[dataset$GMX.s.C1!="D",]
    dataset<-dataset[dataset$GMX.s.C1!="E",]
    dataset<-dataset[dataset$GMX.s.C1!="F",]
    dataset<-dataset[dataset$GMX.s.C1!="G",]
    dataset<-dataset[dataset$GMX.s.C1!="P",]
    dataset<-dataset[dataset$GMX.s.C1!="Q",]
    dataset<-dataset[dataset$GMX.s.C1!="R",]
    dataset<-dataset[dataset$GMX.s.C1!="S",]
    dataset<-dataset[dataset$GMX.s.C1!="T",]
    #dataset<-dataset[dataset$T0.200S!="NaN",]
    Dcensor=numeric(length(dataset$EQID))
    for (i in 1:length(dataset$EQID)) {
      YEAR<-dataset$YEAR[i] 
      Mag<-dataset$Earthquake.Magnitude[i] 
      if (YEAR<=2000) {
        D5=50
        D6=100
        D7=200
      } else if (YEAR>=2001 && YEAR<=2005) {
        D5=100
        D6=150
        D7=250
      } else {
        D5=200
        D6=250
        D7=350
      }
      
      if (Mag<=5) {
        Dcensor[i]=D5
      } else if (Mag>5 && Mag<=6) {
        Dcensor[i]=D5+(D6-D5)*(Mag-5) 
      } else if (Mag>6 && Mag<=7) {
        Dcensor[i]=D5+(D7-D6)*(Mag-6)
      } else {
        Dcensor[i]=D7
      }
    }
    dataset<-cbind(dataset,Dcensor)
    dataset<-dataset[dataset$ClstD..km.<=dataset$Dcensor,]
    
    aa <- table(dataset$Station.Sequence.Number)
    dataset <- dataset[dataset$Station.Sequence.Number %in% names(aa[aa >= 3]),]
    dataset <- dataset[dataset$Quality_Flag!=1,]
    dataset <- dataset[dataset$Spectra.Quality.Flag!=1,]
    dataset <- dataset[dataset$Late.S.trigger!="Y",]
    
    # Output necessary parameters for ASK model
    eqid<-dataset$EQID
    stid<-dataset$Station.Sequence.Number
    year<-dataset$YEAR 
    mag<-dataset$Earthquake.Magnitude
    rrup<-dataset$ClstD..km.
    dip<-dataset$Dip..deg.
    rake<-dataset$Rake.Angle..deg.
    vs30<-dataset$Vs30..m.s..selected.for.analysis
    rx<-dataset$Rx
    w<-dataset$Fault.Rupture.Width..km.
    ztor<-dataset$Depth.to.Top.Of.Fault.Rupture.Model
    rjb<-dataset$Joyner.Boore.Dist...km.
    z1<-dataset$Northern.CA.Southern.CA...S4.Z1..m.
    crjb<-dataset$CRjb
    sof<-dataset$Mechanism.Based.on.Rake.Angle
    gmx<-dataset$GMX.s.C1
    pga<-dataset$PGA..g.
    pgv<-dataset$PGV..cm.sec.
    pgd<-dataset$PGD..cm.
    lon<-dataset$Hypocenter.Longitude..deg.
    lat<-dataset$Hypocenter.Latitude..deg.
    nn <- length(eqid)
    hwflag <- rep(0,nn)
    ry0 <- rep(-999,nn)
    sarock <- rep(0,nn)
    region <- rep(1,nn)
    z1_true=numeric(nn)
    for (i in 1:nn) {
      if (z1[i]<0) {
        z1_true[i] = exp(-7.67/4*log((vs30[i]^4+610^4)/(1360^4+610^4)))/1000
      } else {
        z1_true[i]=z1[i]
      }
    }
    df <- data.frame(eqid,stid,year,mag,rrup,dip,rake,vs30,rx,w,ztor,rjb,z1,z1_true,crjb,sof,gmx,hwflag,ry0,sarock,region,lon,lat)
    df2 <- dataset[,c(135:245)]
    for (i in 1:ncol(df2)) {
      if (as.numeric(substr(colnames(df2)[i],2,6)) < 10) {
        colnames(df2)[i] <- substr(colnames(df2)[i],2,6)
      } else {
        colnames(df2)[i] <- substr(colnames(df2)[i],2,7)
      }
    }
    df <- cbind(df,df2,pga,pgv)
    colnames(df)[length(df)-1] <- format(round(0,3),nsmall=3)
    colnames(df)[length(df)] <- format(round(-1,3),nsmall=3)
    return(df)
  } else if (dataname == "NGA_W2_SC") {
    dataset <- read.csv(filename)
    
    eqid<-dataset$EQID
    stid<-dataset$Station.Sequence.Number
    rgid<-dataset$rgid
    pathid<-dataset$pathid
    year<-dataset$YEAR 
    mag<-dataset$Earthquake.Magnitude
    rrup<-dataset$ClstD..km.
    dip<-dataset$Dip..deg.
    rake<-dataset$Rake.Angle..deg.
    vs30<-dataset$Vs30..m.s..selected.for.analysis
    rx<-dataset$Rx
    w<-dataset$Fault.Rupture.Width..km.
    ztor<-dataset$Depth.to.Top.Of.Fault.Rupture.Model
    rjb<-dataset$Joyner.Boore.Dist...km.
    z1<-dataset$Northern.CA.Southern.CA...S4.Z1..m.
    crjb<-dataset$CRjb
    sof<-dataset$Mechanism.Based.on.Rake.Angle
    gmx<-dataset$GMX.s.C1
    pga<-dataset$PGA..g.
    pgv<-dataset$PGV..cm.sec.
    pgd<-dataset$PGD..cm.
    lon<-dataset$Hypocenter.Longitude..deg.
    lat<-dataset$Hypocenter.Latitude..deg.
    nn <- length(eqid)
    hwflag <- rep(0,nn)
    ry0 <- rep(-999,nn)
    sarock <- rep(0,nn)
    region <- rep(0,nn)
    z1_true <- rep(-1,nn)
    # z1_true=numeric(nn)
    # for (i in 1:nn) {
    #   if (z1[i]<0) {
    #     z1_true[i] = exp(-7.67/4*log((vs30[i]^4+610^4)/(1360^4+610^4)))/1000
    #   } else {
    #     z1_true[i]=z1[i]
    #   }
    # }
    df <- data.frame(eqid,stid,rgid,pathid,year,mag,rrup,dip,rake,vs30,rx,w,ztor,rjb,z1,z1_true,crjb,sof,gmx,hwflag,ry0,sarock,region,lon,lat)
    
    df2 <- dataset[,c(135:245)]
    for (i in 1:ncol(df2)) {
      if (as.numeric(substr(colnames(df2)[i],2,6)) < 10) {
        colnames(df2)[i] <- substr(colnames(df2)[i],2,6)
      } else {
        colnames(df2)[i] <- substr(colnames(df2)[i],2,7)
      }
    }
    df <- cbind(df,df2,pga,pgv)
    colnames(df)[length(df)-1] <- format(round(0,3),nsmall=3)
    colnames(df)[length(df)] <- format(round(-1,3),nsmall=3)
    return(df)
  } else if (dataname == "CyberShake") {
    # Load data file
    dataset <- read.csv(filename)
    
    eqid <- dataset$EQID
    mag <- dataset$Magnitude
    dip <- dataset$Dip
    rake <- dataset$Rake
    ztor <- dataset$ZTOR..km.
    rrup <- dataset$Rrup..km.
    rjb <- dataset$Rjb..km.
    rx <- dataset$Rx..km.
    ry0 <- dataset$Ry0..km.
    z1_true <- dataset$Site.Z1.0..km.
    vs30 <- dataset$Site.Vs30..m.s.
    
    df <- data.frame(eqid,mag,rrup,dip,rake,vs30,rx,ztor,rjb,z1_true,ry0)
    df2 <- dataset[,c(28:32)]
    colnames(df2)[1] <- format(round(2,3),nsmall=3)
    colnames(df2)[2] <- format(round(3,3),nsmall=3)
    colnames(df2)[3] <- format(round(5,3),nsmall=3)
    colnames(df2)[4] <- format(round(7,3),nsmall=3)
    colnames(df2)[5] <- format(round(10,3),nsmall=3)
    df <- cbind(df,df2)
    return(df)
  } else if (dataname == "CyberShake_SC_regression") {
    df <- read.csv(filename)
    colnames(df)[11] <- substr(colnames(df)[11],2,6)
    return(df)
  } else if (dataname == "CyberShake_SC_full") {
    df <- read.csv(filename)
    colnames(df)[11] <- substr(colnames(df)[11],2,6)
    return(df)
  }# End if
} # End function Read_data
