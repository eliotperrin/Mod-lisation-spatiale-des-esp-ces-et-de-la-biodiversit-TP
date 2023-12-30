setwd("/Users/eliotperrin/Desktop/Uni/Biogéosciences/Etude de cas modélisation")
set.seed(209)
# Visualisation des datas
occ<-read.table('/Users/eliotperrin/Desktop/Uni/Biogéosciences/Etude de cas modélisation/Botrychium_lunaria.txt')
head(occ)


# Lecture des fichiers raster
library(terra)
#list files in the Covariate folder
l<- list.files("/Users/eliotperrin/Desktop/Uni/Biogéosciences/Etude de cas modélisation/Covariate",pattern=".tif",full.names=T)
l<-l[-grep("RCP",l)]# remove future layers
env.stack<-rast(l)

env<-extract(env.stack,occ[,1:2],ID=FALSE) #takes some time
occ_env<-cbind(occ,env) # join the occ and env tables
occ_env<-na.omit(occ_env) # remove lines with NAs
occ<-occ_env[,1:3]
head(occ)
env<-occ_env[,-c(1:3)]
head(env)

correlation.matrix<-cor(env, method=c("spearman"))
hc <- hclust(as.dist(1-abs(correlation.matrix)), method="centroid")
plot (hc, sub="", xlab="Variables", hang=-1, ylab="1- abs(correlation)")
abline(0.3,0)
rect.hclust(hc,h=0.3)

print(sum(occ$resp))

# Choix des prédicteurs: bio1, srad_mean_060708_summer, r

var<-c("bio1_19812010","srad_mean_060708_summer","r")
B.pyr<-cbind(occ,env[,var])

# presences to calibrate the model
cal.pres<-sample(which(B.pyr$resp==1),sum(B.pyr$resp==1)*0.8)

# absences to calibrate the model
cal.abs<-sample(which(B.pyr$resp==0),sum(B.pyr$resp==0)*0.8)

#lines of the dataset used for calibration
cal<-c(cal.pres,cal.abs)


formula<-as.formula(paste0("resp ~", paste(var,collapse="+")))
glm1<-glm(formula,data=B.pyr[cal,], family ="binomial")

formula<-as.formula(paste0("resp ~", paste(paste0("poly(",var,",2)"), collapse="+")))
glm2<-glm(formula,data=B.pyr[cal,], family ="binomial")

library(pROC)

# Prédire les valeurs avec le premier modèle
glm.pred1 <- predict(glm1, newdata = B.pyr[-cal, ], type = "response")
roc.object1 <- roc(B.pyr$resp[-cal], glm.pred1)
roc.object1$auc

# Prédire les valeurs avec le deuxième modèle
glm.pred2 <- predict(glm2, newdata = B.pyr[-cal, ], type = "response")
roc.object2 <- roc(B.pyr$resp[-cal], glm.pred2)
roc.object2$auc

# Créer une nouvelle fenêtre graphique
par(mfrow = c(1, 2))  # Une ligne, deux colonnes

# Tracer le premier graphique ROC
plot(roc.object1, main = "ROC Curve - Model 1", col = "blue", lwd = 2)

# Tracer le deuxième graphique ROC sur la deuxième colonne
plot(roc.object2, main = "ROC Curve - Model 2", col = "red", lwd = 2)

source("resplot.R") #load function from file

resplot(glm1,B.pyr)
resplot(glm2,B.pyr)

library(randomForest)

formula<-as.formula(paste0("as.factor(resp) ~", paste(var,collapse="+")))
RF<-randomForest(formula, data= B.pyr[cal,], ntree=1000,importance= TRUE)


importance(RF,type=1)
resplot(RF,B.pyr)

RF.pred01<-predict(RF,newdata=B.pyr[-cal,],type="prob")
RF.pred<-RF.pred01[,2] #probability of presence in 2nd column

roc.object3 <- roc(B.pyr$resp[-cal],RF.pred)
roc.object3$auc

plot(roc.object3, main = "ROC Curve - Random forest", lwd = 2)

pred = predict(env.stack, RF, type="prob", index=2)

par(mfrow = c(1, 1))

plot(pred,main = "RF" )
points(occ[occ$resp==1,])

writeRaster(pred,"prediction.tif",overwrite=T)

var_RCP45<-sub("19812010","20702099_RCP45",var)
var_RCP85<-sub("19812010","20702099_RCP85",var)


env.stack.RCP45<-rast(paste0("Covariate/",var_RCP45,".tif"))
env.stack.RCP85<-rast(paste0("Covariate/",var_RCP85,".tif"))

names(env.stack.RCP45)<-var
names(env.stack.RCP85)<-var

pred_RCP45<-predict(env.stack.RCP45,RF, type="prob", index=2)
plot(pred_RCP45,main = "RF RCP45")
writeRaster(pred_RCP45,"prediction_RCP45.tif", overwrite=T)

pred_RCP85<-predict(env.stack.RCP85,RF, type="prob", index=2)
plot(pred_RCP85,main = "RF RCP85")
writeRaster(pred_RCP85,"prediction_RCP85.tif", overwrite=T)

boxplot(c(pred,pred_RCP45,pred_RCP85),names=c("présent","RCP 4.5","RCP 8.5"))

diff_RCP45 <- pred_RCP45-pred

writeRaster(diff_RCP45,"diff_RCP45.tif",overwrite=T)

prev<-sum(occ$resp)/length(occ$resp)
pred_bin <- pred>prev
pred_RCP45_bin <- pred_RCP45>prev
plot(pred_RCP45_bin, main = "RCP45 Binary predictions")

change.RCP45<- pred_RCP45_bin - 2*pred_bin
plot(change.RCP45, main = "Changes RCP45")

writeRaster(change.RCP45,"chg_RCP45.tif",overwrite=T)

diff_RCP85 <- pred_RCP85-pred

writeRaster(diff_RCP85,"diff_RCP85.tif",overwrite=T)

prev<-sum(occ$resp)/length(occ$resp)
pred_bin <- pred>prev
pred_RCP85_bin <- pred_RCP85>prev
plot(pred_RCP85_bin, main = "RCP85 Binary predictions")

change.RCP85<- pred_RCP85_bin - 2*pred_bin
plot(change.RCP85, main = "Changes RCP85")

writeRaster(change.RCP85,"chg_RCP85.tif",overwrite=T)