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
library(ggplot2)
library(ggmap)

con <- dbConnect(MySQL(), dbname = "CyberShake", user = "cybershk_ro", host = "focal.usc.edu", password = "CyberShake2007")

# rs = dbGetQuery(con,statement = paste ("select P.Source_ID, P.Rupture_ID, P.Rup_Var_ID, P.IM_Value",
#                   "from CyberShake_Sites S, CyberShake_Runs R, PeakAmplitudes P",
#                   "where S.CS_Short_Name=\"LADT\"",
#                   "and S.CS_Site_ID=R.Site_ID",
#                   "and R.Max_Frequency=10.0 and R.ERF_ID=36 and R.Status=\"Verified\"", 
#                   "and R.SGT_Variation_ID=8 and R.Velocity_Model_ID=5 and R.Rup_Var_Scenario_ID=6",
#                   "and P.Run_ID=R.Run_ID",
#                   "and P.IM_Type_ID=21",
#                   "order by P.IM_Value desc limit 100"))


rs = dbGetQuery(con,statement = paste ("select CS_Short_Name,CS_Site_Lon,CS_Site_Lat",  
                                       "from CyberShake_Sites S", 
                                       "where S.CS_Site_Lat>32.6 AND S.CS_Site_Lat<35.2", 
                                       "AND S.CS_Site_Lon>-119.2 AND S.CS_Site_Lon<-116", 
                                       "AND S.CS_Site_Type_ID=7"))

for (i in 1:nrow(rs)) {
#for (i in 3:3) {
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
  
  #    data2 <- data.frame(matrix(-999, nrow = nrow(data1), ncol = 5))
  #    colnames(data2) <- c("Rrup","Mag","Hypocenter_Lon","Hypocenter_Lat","Hypocenter_Dep")
  #     for (j in 1:nrow(data1)) {
  #           src_id <- data1$Source_ID[j]
  #           rup_id <- data1$Rupture_ID[j]
  #           var_id <- data1$Rup_Var_ID[j]
  #           sta_id <- data1$CS_site_ID[j]
  #           data3 = dbGetQuery(con,statement=paste0 ("select CS.Site_Rupture_Dist, RU.Mag, RV.Hypocenter_Lon, RV.Hypocenter_Lat, RV.Hypocenter_Depth",
  #                                                   " from CyberShake_Site_Ruptures CS, Rupture_Variations RV, Ruptures RU",
  #                                                   " where RV.ERF_ID=36 and RU.ERF_ID=36 and CS.ERF_ID=36",
  #                                                   " and RV.Rup_Var_Scenario_ID=6 and CS_Site_ID=",sta_id,
  #                                                   " and RV.Source_ID=",src_id," and RU.Source_ID=",src_id," and CS.Source_ID=",src_id,
  #                                                   " and RV.Rupture_ID=",rup_id," and RU.Rupture_ID=",rup_id," and CS.Rupture_ID=",rup_id, 
  #                                                   " and RV.Rup_Var_ID=",var_id))
  #           data2$Rrup[j] <- data3$Site_Rupture_Dist
  #           data2$Mag[j] <- data3$Mag
  #           data2$Hypocenter_Lon[j] <- data3$Hypocenter_Lon
  #           data2$Hypocenter_Lat[j] <- data3$Hypocenter_Lat
  #           data2$Hypocenter_Dep[j] <- data3$Hypocenter_Depth
  #     }
  # data1 <- cbind(data1, data2)
  data1 <- na.omit(data1[data1$Site_Rupture_Dist>0 & data1$Site_Rupture_Dist<200,])
  write.csv(data1,file=paste0("/Users/xmeng/Documents/MyWork/CyberShake/Old/15_4/",sta,".csv"))
}

on.exit(dbDisconnect(con))


ca <- get_map(location = c(-119.2, 32.6, -116, 34.8), maptype = c("terrain"),color = "bw")
gg <- ggmap(ca)
p1 <- gg + geom_point(data=rs,aes(x=CS_Site_Lon,y=CS_Site_Lat),
                      shape=21, color="red",show.legend =FALSE)
