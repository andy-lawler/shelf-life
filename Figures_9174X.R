############
## 9174X ###
############
source("shelflife_function.R")
stdata <- read.csv("9174X Stability Data.csv", header=T, na=".")
stdata <- na.omit(stdata)

response <- "Result"
colnames(stdata)[1] <- "Month"
stdata$Batch <- as.factor(stdata$Batch)

UU <- 55
LL <- 35
var1 <- ((UU-LL)/4)^2
var2 <- ((UU-LL)/4/2.5)^2
sigmab <- sd(stdata[,response], na.rm=TRUE)


### random effect model
library(lme4)
rfit <- lmer(Result ~ Month + (1 + Month|Batch), data = stdata) 
summary(rfit)

apply(coef(rfit)[[1]], 2, sd)

### ICH method
DT <- as.data.table(stdata[!is.na(stdata[,response]), ])
colnames(DT)[3] <- "batch"
colnames(DT)[1] <- "month"
colnames(DT)[which(colnames(stdata)==response)] <- "resp"

ich_osm <- ich_shelf(DT, byvar="batch", spec.up=UU, spec.lo=LL, xmax=200, xmin=-100)
ich_osm
tb_osm <- ich_osm$res[,c(1,2,3,5,11,12,13,14,19)]
colnames(tb_osm)[9] <- "shelf_life"
tb_osm

### Baymethod
#fit0 <- bay_shelflife(stdata, response=response, spec_up=UU, spec_lo=LL)  ##default: the same as SAS
#fit1 <- bay_shelflife(stdata, response=response, spec_up=UU, spec_lo=LL, sigmab=sigmab) 
set.seed(1500)
fit2 <- bay_shelflife(stdata, response=response, spec_up=UU, spec_lo=LL, var2=((UU-LL)/4/3.5)^2)
#fit3 <- bay_shelflife(stdata, response=response, spec_up=UU, spec_lo=LL, var2=((UU-LL)/4/4)^2)
#quantile(fit0$sl, 0.05, na.rm=TRUE)
#quantile(fit1$sl, 0.05, na.rm=TRUE)
quantile(fit2$sl, 0.05, na.rm=TRUE)
#quantile(fit3$sl, 0.05, na.rm=TRUE)

regress <- fit2$st.out$BUGSoutput$sims.list[c("mu.int", "mu.slope")]
##################
### figures 4 ####
##################
tiff("Figure4.tiff", width = 8, height = 8, units = 'in', res = 600)
par(mfrow=c(2, 1), mar=c(4, 4, 3, 1))
plot(resp~month, data=DT, pch=16, col="gray50", cex=0.7, ylim=c(32, 55), xlim=c(0, 80), ylab="Assey", type="n",
     main="Distribution of mean responses over time (take 200 sample lines as example)")
# for(ii in 1:nrow(tb_osm))abline(tb_osm[ii, 2], tb_osm[ii, 3], col=ii)
abline(h=c(35, 55), col=2, lwd=2)
set.seed(2000)
aa <- sample(1:4000, 200)
for(ii in  aa){
  abline(regress[[1]][ii], regress[[2]][ii], col=rgb(0, 0, 0, alph=0.2), lty=1)
  if(regress[[2]][ii]>0) points((55-regress[[1]][ii])/regress[[2]][ii], 55, pch=16, cex=0.5)
  else   points((35-regress[[1]][ii])/regress[[2]][ii], 35, pch=16, cex=0.5)
}

#points((35-regress[[1]][aa])/regress[[2]][aa], rep(35, each=250), pch=16, cex=0.5)

hist(fit2$sl[aa], xlim=c(0, 80), breaks=c(seq(0, 355695, by=5)), main="Histogram of shelf-lives", xlab="shelf-life (mounth)")
abline(v=quantile(fit2$sl, 0.05, na.rm=TRUE))
text(23.5, 19, "5th percentile\n(shelf-life)", adj=c(1, 0.5))
dev.off()

################
### figure 1 ###
################
tiff("Figure1.tiff", width = 8, height = 6, units = 'in', res = 600)
###1 Explain the shelf-life of ICH 
tempDT <- subset(DT, batch=="62016")
fit <- lm(resp~month, data=tempDT)

par(mar=c(4, 4, 2, 1))
plot(resp~month, data=tempDT, pch=16, cex=0.7, ylim=c(28, 57), xlim=c(0, 40), type="n", 
     ylab="Assey", xlab="Time Points (Months)")
# title("The shelf-life estimated by ICH approach")
abline(h=c(LL, UU),lty=2, lwd=2)

newx<-seq(-1, 41)

prd<-predict(fit,newdata=data.frame(month=newx),interval = c("confidence"), 
             level = 0.95,type="response")
abline(fit)
lines(newx,prd[,2],lty=3)
lines(newx,prd[,3],lty=3)
legend(20, 54, c("acceptance criteria", "mean response(regression line)", "95% CI of mean response"), lty=c(2, 1, 3), lwd=c(2, 1, 1), col=1)
points(30.3, 35, pch=16)
text(30.3, 34.5, "shelf-life", adj=c(0.5, 1))
dev.off()
###############
#############
plot(Result~Month, data=stdata, pch=16, col=as.numeric(Batch), cex=0.7, ylim=c(32, 55))
for(ii in 1:nlevels(stdata$Batch)){
  temp <- subset(stdata, Batch==levels(stdata$Batch)[ii])
  lines(lowess(temp$Result ~ temp$Month), col=ii, lty=ii)
}
lines(lowess(stdata$Result ~ stdata$Month), col = 2, lwd=2)

abline(h=c(35, 55), col="gray50")
#####################

layout(matrix(c(1:16, rep(17, 16)), nrow=4))
par(mar=c(0.01, 0.01, 0.01, 0.01))
for(ii in 1:nlevels(stdata$Batch)){
  temp <- subset(stdata, Batch==levels(stdata$Batch)[ii])
  plot(Result~Month, data=temp, pch=16, cex=0.7, ylim=c(32, 55), xlim=c(0, 36), yaxt="n", xaxt="n")
  lines(lowess(temp$Result ~ temp$Month), col=ii, lty=ii)
  abline(h=c(35, 55), col="gray50")
}

par(mar=c(5, 5, 3, 1))
plot(Result~Month, data=stdata, pch=16, col=as.numeric(Batch), cex=0.7, ylim=c(32, 55), type="n", ylab="")
for(ii in 1:nlevels(stdata$Batch)){
  temp <- subset(stdata, Batch==levels(stdata$Batch)[ii])
  lines(lowess(temp$Result ~ temp$Month), col=ii, lty=ii)
}
lines(lowess(stdata$Result ~ stdata$Month), col = 2, lwd=2)

abline(h=c(35, 55), col="gray50")

########################
### regression lines ###
########################

plot(Result~Month, data=stdata, pch=16, col=as.numeric(Batch), cex=0.7, ylim=c(32, 55), xlim=c(0, 40), type="n", ylab="CLO")
for(bb in 1:nrow(tb_osm)) abline(tb_osm[bb, 2], tb_osm[bb, 3], col="blue")
abline(h= c(35, 55))

########
fit3$sl
bayOutput <- fit3

diagPlot(bayOutput, "mu.int")
diagPlot(bayOutput, "mu.slope")
diagPlot(bayOutput, "sigma.e")

sims <- bayOutput$st.out$BUGSoutput$sims.array
sims <- aperm(sims, c(1, 3, 2))

plot(Result~Month, data=stdata, pch=16, col=as.numeric(Batch), cex=0.7, ylim=c(32, 55), xlim=c(0, 100), type="n", ylab="")
#for(kk in 1:3) {
    for(bb in 1:1000) abline(sims[bb, "mu.int", 1], sims[bb, "mu.slope", 1], col=rgb(0, 0, 0, alph=0.1))
#}
#for(bb in 1:nrow(tb_osm))abline(tb_osm[bb, 2], tb_osm[bb, 3], col=rgb(1, 0, 0))
#abline(h= c(35, 55), v=0)

temp <- fit3$sl[fit3$sl<200]
hist(temp, breaks=seq(0, 200, by=2), xaxt="n", col=c("yellow", "orange", "orange3", "red3", "red"))
axis(1, seq(0, 150, by=10), seq(0, 150, by=10))
abline(v=24.2)
