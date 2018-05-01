# Clear the console
cat("\014")  
# Clear the workspace
rm(list=ls())

sta_list <- c("BAK","BCW","BLM","CAR","ECH","FIG","GATR","LMY","LOWE","MPP","OCUT","OSI",
              "SCYB","SMM","STC","SYP","TAFT","TRAM","TSCN","UCSB")

CS <- data.frame(Date=as.Date("01/01/2000", format="%m/%d/%Y"), 
                    File="", User="", stringsAsFactors=FALSE)
CS <- CS[-1,]
for (i in 1:length(sta_list)) {
  sta <- sta_list[i]
  filename <- paste("/Users/xmeng/Documents/MyWork/CyberShake/Old/CyberShake/",sta,".csv",sep="")
  filename2 <- paste("/Users/xmeng/Documents/MyWork/CyberShake/Old/CyberShake/",sta,"_00.csv",sep="")
  CC <- read.csv(filename2)
  CC <- CC[CC$Rupture.ID==0,]
  #CC$EQID <- as.numeric((formatC(CC$Source.ID,width=3,format="d",flag="0")))
  #CC$EQID <- paste(formatC(CC$Source.ID,width=3,format="d",flag="0"),
  #                 formatC(CC$Rupture.ID,width=4,format="d",flag="0"),sep="")
  #CC <- CC[CC$Rupture.Variation.ID==0,]
  #write.csv(CC,filename2)
  
  CS <- rbind(CS,CC)
}
CS <- CS[CS$EQID%%2==0,]
#write.csv(CS,"/Users/xmeng/Documents/MyWork/CyberShake/Old/CyberShake/CyberShake.csv")

library(nlme)
library(plyr)
library(reshape2)
library(ggplot2)


eqId<-CS$EQID
mag<-CS$Magnitude
dist<-CS$Rrup..km.
obs<-log(CS$X2s.SA.RotD50..cm.s.s./100)
selected<-data.frame(eqId,mag,dist,obs)
selected<-selected[order(selected$eqId),]
selected<-selected[selected$eqId%%1==0,]

# Functions defining the components of the model
magAndDistanceComponents <- function(mag, dist, c1, c2, c4, c6, c7, c10, c11) {
  return(c1 + c2 * mag + (c6 + c7 * mag) * log(dist + exp(c4)) +
           c10 * (mag - 6) ^ 2 + c11 * dist)
}

siteComponent <- function(region, vs30, c12g, c12n) {
  vsRef = 3000.
  
  return(ifelse(region == 'glacial', c12g, c12n) * log(vs30 / vsRef))
}

#PGM <- data.frame(matrix(NA, nrow = nrow(selected), ncol = 1))
#for (i in 1:nrow(selected)) {
#  PGM[i,1] <- magAndDistanceComponents(selected$mag[i], selected$dist[i],-0.06,0.28,0.60,-2.06,0.20,-0.29,-0.01)
#}


# control = nlmeControl(pnlsTol = 0.02, msVerbose = TRUE)
 fitModel <- nlme(obs ~ magAndDistanceComponents(mag, dist, c1, c2, c4, c6, c7, c10, c11),
                  data = selected,
                  fixed = c1 + c2 + c4 + c6 + c7 + c10 + c11 ~ 1,
                  random = c1 ~ 1|eqId,
                  start=c(c1=-0.06, c2=0.28, c4=0.6, c6=-2.06, c7=0.2, c10=-0.29,
                          c11=-0.01)
 )

#   fitModel <- function(selected) {
#        nlme(obs ~ magAndDistanceComponents(mag, dist, c1, c2, c4, c6, c7, c10, c11),
#        data = selected,
#        fixed = c1 + c2 + c4 + c6 + c7 + c10 + c11 ~ 1,
#        random = c1 ~ 1|eqId,
#        start=c(c1=0., c2=0.5, c4=2.8, c6=-3.0, c7=0.2, c10=-0.1,
#                c11=-0.02)
#     )
#   }
# 
# 
# models <- dlply(selected, 'obs', failwith(NULL, fitModel), .progress='text')
# 
# # Remove NULL entries
# models <- compact(models)
# coefs <- dlply(models, fixed.effects)