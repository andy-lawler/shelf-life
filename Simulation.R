x <- matrix(rnorm(20000, 0, 1), nrow=2)
mean(apply(x, 2, sd)<1)

#####################
### generate data ###
#####################
setwd("C:/Users/alawle/Documents/P_ABC_DrugStability/Rfigure")

source("shelflife_function.R")
response <- "Response"

if(FALSE){
  ### Setting 0:Lem
  UU <- 260; LL <- 200
  alphaDist <- c(230, sqrt(5))
  betaDist <- c(1, sqrt(0.25))
  sigma <- sqrt(5)
}


if(FALSE){
  ### Setting 0:Lem2
  UU <- 260; LL <- 200
  alphaDist <- c(230, sqrt(5))
  betaDist <- c(1, 0.25)
  sigma <- sqrt(5)
}

if(FALSE){
  ### Setting 0:Lem3
  UU <- 260; LL <- 200
  alphaDist <- c(230, 5)
  betaDist <- c(1, sqrt(0.25))
  sigma <- sqrt(5)
}

if(FALSE){
  ### setting 2: RefQ
  UU <- 110; LL <- 90
  alphaDist <- c(101, sqrt(1.5))
  betaDist <- c(-0.33, sqrt(0.0015))
  sigma <- sqrt(0.5)
}

if(FALSE){
  ### setting 1:CLO   ### AAA
  UU <- 55; LL <- 35
  alphaDist <- c(45, 2)
  betaDist <- c(-0.2, 0.12)
  sigma <- 2.5
}

if(FALSE){
  ### setting 4: BUN
  UU <- 0.55*100; LL <- 0.45*100
  alphaDist <- c(0.51, 0.0033)*100
  betaDist <- c(0.0003, 0.0008)*100
  sigma <- 0.00875*100
}

#### get true 5% batch shelflife
alpha <- rnorm(10^3, alphaDist[1], alphaDist[2])
beta <- rnorm(10^3, betaDist[1], betaDist[2])
batchlife <- ifelse(beta>0, (UU-alpha)/beta, (LL-alpha)/beta)
quantile(batchlife, c(0.5, .05))


if(FALSE){
  ### setting: OSM  ###BBB
  library(MASS)
  UU <- 300; LL <- 250
  Mu <- c(275, 0.7)
  Cov <- diag(c(8, 0.25))
  Cov[1, 2] <- Cov[2, 1] <- sqrt(Cov[1, 1])*sqrt(Cov[2, 2])*(-0.50)
  sigma <- 4
  lines <- mvrnorm(1000, Mu, Cov)
  alpha <- lines[, 1]
  beta <- lines[, 2]
}

batchlife <- ifelse(beta>0, (UU-alpha)/beta, (LL-alpha)/beta)
#batchlife <- ifelse(batchlife>60, 60, batchlife)
quantile(batchlife, c(0.5, .05))

################
yobs <- NULL
n.size <- 2
month <- rep(c(0, 3, 6, 9, 12, 18, 24, 36), each=n.size)
for(batch in 1:10^3){
  y <- alpha[batch]+beta[batch]*month+rnorm(length(month), 0, sigma)
  yobs <- c(yobs, y)
}
DT <- data.frame(Y=yobs, Batch=rep(1:1000, each=length(month)), Month=month)
################
### OBS. shelf life (use mean)
Temp <- sapply(1:10^3, function(ii){
  dt <- DT[DT$Batch==ii, ]
  coef(lm(Y~Month, data=dt))
})
obs.batchlife <- ifelse(Temp[2, ]>0, (UU-Temp[1, ])/Temp[2, ], (LL-Temp[1, ])/Temp[2, ])

### shelf life (use 95%)
DT_ich <- as.data.table(DT)
colnames(DT_ich) <- c("resp", "batch", "month")
obs.batchlife2 <- ich_shelf(DT_ich, byvar="batch", spec.up=UU, spec.lo=LL, xmax=200, xmin=-70)$res$shelf_life2
###################
# c(3, 6, 10, 15, 20, 30, 40)
for(n.batch in c(3, 6, 10, 15, 20, 30, 40)){
  bayresult0 <- bayresult1 <- resultich <- resultichall <- ichresult <- NULL
  for(rr in 1:100){
    stdata <- DT[DT$Batch %in% sample(1:1000, n.batch), ]
    stdata_ich <- as.data.table(stdata)
    
    ###### ICH method
    colnames(stdata_ich) <- c("resp", "batch", "month")
    ## poolability ##
    pool2 <- FALSE
    pool1 <- coef(summary(lm(resp~month*batch, data=stdata_ich)))[4, 4] > 0.25
    if(pool1){
      pool2 <- coef(summary(lm(resp~month+batch, data=stdata_ich)))[3, 4] > 0.25
    }
    ## not pool ##
    ich <- ich_shelf(stdata_ich, byvar="batch", spec.up=UU, spec.lo=LL, xmax=200, xmin=-70)
    tb_ich <- ich$res[,c(1,2,3,5,11,12,13,14,19)]
    colnames(tb_ich)[9] <- "shelf_life"
    ## pool ##
    stdata_ich$batch2 <- "All"
    ichall <- ich_shelf(stdata_ich, byvar="batch2", spec.up=UU, spec.lo=LL, xmax=200, xmin=-70)
    tb_ichall <- ichall$res[,c(1,2,3,5,11,12,13,14,19)]
    colnames(tb_ichall)[9] <- "shelf_life"
    
    ichresult[rr] <- ifelse(pool2, tb_ichall[, "shelf_life"], min(tb_ich[, "shelf_life"]))
    
    
    stdata$Batch <- as.factor(stdata$Batch)
    fit0 <- bay_shelflife(stdata, response="Y", spec_up=UU, spec_lo=LL)  
    fit1 <- bay_shelflife(stdata, response="Y", spec_up=UU, spec_lo=LL, sigmab=sd(stdata$Y)) 
    # fit2 <- bay_shelflife(stdata, response="Y", spec_up=UU, spec_lo=LL, var2=((UU-LL)/4)^2)
    # fit3 <- bay_shelflife(stdata, response="Y", spec_up=UU, spec_lo=LL, sigmab=sd(stdata$Y), var2=((UU-LL)/4/4)^2)  ##apply different tc
    
    bayresult0 <- cbind(bayresult0, as.numeric(fit0$sl))
    bayresult1 <- cbind(bayresult1, as.numeric(fit1$sl))
    resultich[[rr]] <- tb_ich
    resultichall[[rr]] <- tb_ichall
  }
  
  result <- cbind(BAY_SAS=as.numeric(fit0$sl), BAY_R=as.numeric(fit1$sl))
  BAY <- t(apply(bayresult0, 2, function(x) quantile(x, c(0.05, 0.25, 0.5), na.rm=TRUE)))
  BAYR <- t(apply(bayresult1, 2, function(x) quantile(x, c(0.05, 0.25, 0.5), na.rm=TRUE)))
  ICH1<- apply(sapply(resultich, function(x)x[, 9]), 2, min)
  ICH2<- sapply(resultichall, function(x)x[, 9])
  if(FALSE){
  layout(matrix(c(2, 2, 1, 3, 3, 3), ncol=6))
  boxplot(cbind(batchlife, obs.batchlife, obs.batchlife2),  ylim=c(0, 150), main=paste0("true shelf-lives(", length(batchlife), ")"))
  points(c(quantile(batchlife, c(0.05)), quantile(obs.batchlife, c(0.05)), quantile(obs.batchlife2, c(0.05))), col=2, pch=16)
  
  boxplot(result,  ylim=c(0, 150), main=paste0("n.batch = ", n.batch, ", n.size = ", n.size, ", n.sample = ", nrow(result)))
  points(apply(result, 2, function(x) quantile(x, c(0.05), na.rm=TRUE)), col=2, pch=16)
  boxplot(cbind(BAY, worst=ICH1, pool=ICH2, ichresult),  ylim=c(0, 150),
          main=paste("replication:", ncol(bayresult0), "times"))}
  save(list=ls(), file=paste0("simulation/OSM_n.batch", n.batch , ".Rdata"))
}
