# Clear the console
cat("\014")  
# Clear the workspace
rm(list=ls())


#Runs produced with Study 15.4 have 
#SGT_Variation_ID=8, Rup_Var_Scenario_ID=6, ERF_ID=36, Low_Frequency_Cutoff=1.0, Velocity_Model_ID=5 
#in the CyberShake database.

#select CS_Site_Name,CS_Site_Type_ID  from CyberShake_Sites S 
#where S.CS_Site_Lat>32.6 AND S.CS_Site_Lat<34.8 AND S.CS_Site_Lon>-119.2 AND 
#S.CS_Site_Lon<-116 AND S.CS_Site_Type_ID<7;


require(RMySQL)

con <- dbConnect(MySQL(), dbname = "CyberShake", user = "cybershk_ro", host = "focal.usc.edu", password = "CyberShake2007")
rs = dbGetQuery(con,statement = paste ("select CS_Short_Name,CS_Site_Lon,CS_Site_Lat",  
                                       "from CyberShake_Sites S", 
                                       "where S.CS_Site_Lat>32.6 AND S.CS_Site_Lat<35.2", 
                                       "AND S.CS_Site_Lon>-119.2 AND S.CS_Site_Lon<-116"))

for (i in 1:nrow(rs)) {
  print(i)
  sta <- rs$CS_Short_Name[i]
  data1 = dbGetQuery(con,statement = paste0 ("select RV.Hypocenter_Lon, RV.Hypocenter_Lat, RV.Hypocenter_Depth,", 
                     " CS.Site_Rupture_Dist, RU.Mag, P.Source_ID, P.Rupture_ID, P.Rup_Var_ID, P.IM_Value,",
                     " S.CS_Short_Name, S.CS_site_ID, S.CS_Site_Lon, S.CS_Site_Lat, R.Vs30, P.Run_ID",
                     " from CyberShake_Sites S, CyberShake_Runs R, PeakAmplitudes P,",
                     " Rupture_Variations RV, CyberShake_Site_Ruptures CS, Ruptures RU",
                     " where S.CS_Short_Name=\"",sta,"\"",
                     " and S.CS_Site_ID=R.Site_ID",
                     " and R.Low_Frequency_Cutoff=1.0 and R.ERF_ID=36 and R.Status=\"Verified\"", 
                     " and R.SGT_Variation_ID=8 and R.Velocity_Model_ID=5 and R.Rup_Var_Scenario_ID=6",
                     " and P.Run_ID=R.Run_ID and P.IM_Type_ID=21",
                     " and RU.ERF_ID=36 and RU.Source_ID=P.Source_ID and RU.Rupture_ID=P.Rupture_ID",
                     " and CS.CS_Site_ID=S.CS_site_ID and CS.ERF_ID=36 and CS.Source_ID=P.Source_ID and CS.Rupture_ID=P.Rupture_ID",
                     " and RV.ERF_ID=36 and RV.Rup_Var_Scenario_ID=6 and RV.Source_ID=P.Source_ID and RV.Rupture_ID=P.Rupture_ID and RV.Rup_Var_ID=P.Rup_Var_ID"))
  
  data1 <- na.omit(data1[data1$Site_Rupture_Dist>0 & data1$Site_Rupture_Dist<200,])
  write.csv(data1,file=paste0("/home/rcf-40/xiaofenm//CyberShake/Old/15_4/",sta,".csv"))
}

on.exit(dbDisconnect(con))

