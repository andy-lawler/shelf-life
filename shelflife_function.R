####################################################
### Shelf Life Using Bayesian Hierarchical Model ###
####################################################
library(R2jags)
library(ggplot2)
# Define the Function to Calculate Bayesian Shelf Life
bay_shelflife <- function(stdata, response, spec_up, spec_lo, sigmab=NULL, var1=NULL, var2=NULL, nchain=2){
  # Preprocess stdata to Remove NA
  # stdata <- stdata2540
  # response="POL80"
  # sigmab=10; var1=0.0056; var2=0.00046; spec_up=0.65; spec_lo=0.35; nchain=3
  stdata0 <- stdata
  stdata0 <- stdata0[!is.na(stdata0[,response]),]
  if(is.null(sigmab)){
    ### SAS option
    data <- na.omit(stdata[, c("Batch", "Month", response)])
    fit <- lm(paste0(response, " ~ Month*Batch"), data=data)
    sigmab <- summary(fit)$sigma*3
    ### R option
    #sigmab <- sd(stdata0[,response])
  }
  if(is.null(var1)){
    var1 <-((spec_up-spec_lo)/4)^2
  }
  if(is.null(var2)){
    n.batch <- length(unique(stdata$Batch))
    fittedtc <- c(30, 23.2, 17.5, 13, 9.9, 8.1, 7, 6, 5, 4.4, 4, 3.8, 3.5, 3.2, 3, 2.9, 2.8, 
      2.7, 2.7, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2, 1.9, 1.8, 1.7, 1.7, 1.6, 
      1.5, 1.5, 1.4, 1.4, 1.4, 1.3, 1.3) 
    names(fittedtc) <- 3:40
    tc <- fittedtc[paste(ifelse(n.batch>40, 40, n.batch))]
    var2 <-((spec_up-spec_lo)/(4*tc))^2
  }
  # Define the BUGS Model  ### generate to "lmeSW.txt
  cat("
      model{
      # Random intercept and slope
      mu.int ~ dnorm(0, 1.0E-5)  
      mu.slope ~ dnorm(0, 1.0E-5) 
      # sigma.e ~ dgamma(2, 1.0)
      sigma.e ~ dunif(0,", sigmab,")
      tau.e <- pow(sigma.e,-2)
      
      # variance-covariance matrix of subject ranefs
      invSigma.u ~ dwish(R.u , 3)
      R.u[1,1] <- ", var1*3,"
      R.u[2,2] <- ", var2*3,"
      R.u[1,2] <- 0
      R.u[2,1] <- 0
      
      # Intercept and slope for each obs, including random effects
      for(j in 1:ngroups){
      b.hat[j, 1] <- mu.int
      b.hat[j, 2] <- mu.slope
      u[j,1:2] ~ dmnorm(b.hat[j, ], invSigma.u)   ###invSigma.u or Sigma.u ### check input of dmnorm
      b0[j] <- u[j, 1] # random intercept
      b1[j] <- u[j, 2] # random slope
      }
      
      # Likelihood
      for(i in 1:n){
      mu[i] <- b0[batch[i]] + b1[batch[i]] * x[i]
      y[i] ~ dnorm(mu[i], tau.e)
      }    
      }", file="lmeSW.txt")

  # Parameters to Save
  params <- c("mu.int", "mu.slope", "sigma.e", "invSigma.u")
  
  # Define the Data for the JAGS function
  st.data <- list(y=as.numeric(stdata0[,response]), batch=as.numeric(stdata0$Batch), x=as.numeric(stdata0$Month),
                  ngroups=length(levels(stdata0$Batch)), n=nrow(stdata0))
  
  start.time <- proc.time() 
  st.out <- jags(st.data, parameters.to.save=params, model.file="lmeSW.txt",
                 n.thin=5, n.chains=nchain, n.burnin=3000, n.iter=13000)
  end.time <- proc.time()
  elapse_time <- end.time - start.time
  
  theta1 <- st.out$BUGSoutput$sims.array[,,"mu.int"]
  theta2 <- st.out$BUGSoutput$sims.array[,,"mu.slope"]
  # Get Indeces where Intercept is below Lower Limit or Over Upper Limit
  oul_idx <- which(theta1 > spec_up, arr.ind = TRUE)
  bll_idx <- which(theta1 < spec_lo, arr.ind = TRUE)
  ob_idx <- rbind(oul_idx, bll_idx)
  print(head(ob_idx))
  x1 <- (spec_up - theta1) / theta2
  x2 <- (spec_lo - theta1) / theta2
  # Intercept Out of Bounds Cases are Excluded
  x1[ob_idx] <- NA
  x2[ob_idx] <- NA
  sl <- ifelse(theta2 > 0, x1, x2)
  if (length(dim(sl)) == 0) baysl <- quantile(sl, c(.05, .50, .95), na.rm=TRUE)
  else baysl <- apply(sl, 2, function(x) quantile(x, c(.05, .50, .95), na.rm=TRUE))
  list(baysl=baysl, st.out=st.out, elapse_time=elapse_time, ob_idx=ob_idx, sl=sl)
}

diagPlot <- function(bayOutput, myparam, mm=matrix(c(1, 2, 1,3), 2, 2)){
  if(myparam=="sl"){
     param1<- log(bayOutput$sl)
     n.chain <- ncol(param1)
  } else{
    sims <- bayOutput$st.out$BUGSoutput$sims.array
    sims <- aperm(sims, c(1, 3, 2))
    n.chain <- dim(sims)[3]
    if(is.numeric(myparam)) myparam <- dimnames(sims)[[2]][myparam]
    param1 <- data.frame(sims[,myparam,])
  }
  if(!is.null(mm)) layout(mm)
  par(mar=c(4, 4, 3, 2))  
  Color <- c("red3", "skyblue2", "green2", "gray70", "orange")
  ColorLine <- c("darkred", "blue3", "green4", "gray50", "orange3")
  LW <- sapply(1:n.chain, function(x)lowess(1:nrow(param1), param1[, x], f=1/200)[[2]]) 
  matplot(param1, type='l', col=Color, xlab="Iteration", ylab=gsub("theta", "mu", myparam))
  matplot(LW, type='l', col=ColorLine, add=TRUE, lwd=2, lty=1)
  title(paste("Diagnostics for", gsub("theta", "mu", myparam), "based on", n.chain, "chains"))
  
  par(mar=c(4, 4, 0.5, 2)) 
  bacf <- acf(param1[, 1], main="", ylab="Autocorrelation")
  
  ymax <- sapply(1:n.chain, function(x) max(density(param1[, x], adjust=2)$y))
  plot(density(param1[, which.max(ymax)], adjust=2), col=Color[which.max(ymax)], 
       main="", type="n", ylab="Posterior Density", xlab=gsub("theta", "mu", myparam))
  for(ii in 1:n.chain)
    lines(density(param1[,ii], adjust=2), col=ColorLine[ii])
}


# Define the Function to Generate Diagnostic Plots
diagPlot2 <- function(bayOutput, myparam){
  if(myparam=="sl"){
    param1 <- data.frame(log(bayOutput$sl))
  } else{
    st.out <- bayOutput$st.out
    param1 <- data.frame(st.out$BUGSoutput$sims.array[,,myparam])
  }
  colnames(param1) <- c("chain1", "chain2", "chain3")
  param1s <- stack(param1)
  param1 <- cbind(1:nrow(param1), param1)
  colnames(param1)[1] <- "order"
  ciline <- qnorm((1 - 0.95)/2)/sqrt(nrow(param1))
  
  bacf <- acf(param1[,2], plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  
  dplot <- ggplot(param1s, aes(values, colour=ind)) +
    geom_density(adjust=2) + 
    theme(legend.position="none") + 
    xlab(myparam)
  
  tplot <- ggplot() + 
    geom_line(data = param1, aes(x = order, y = chain2, color = "blue")) + 
    geom_line(data = param1, aes(x = order, y = chain3, color = "green")) + 
    geom_line(data = param1, aes(x = order, y = chain1, color = "red")) + 
    theme(legend.position="none") +
    xlab("Iteration") +
    ylab(myparam)
  
  acfplot <- ggplot(data = bacfdf) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(x=lag, xend=lag, y=acf, yend=0))  + 
    geom_hline(yintercept = ciline, color = "blue", size = 0.2, linetype="dashed") +
    geom_hline(yintercept = -ciline, color = "blue", size = 0.2, linetype="dashed") +
    xlab("Lag") +
    ylab("ACF")
  
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  
  multiplot(tplot, acfplot, dplot, layout=matrix(c(1,2,1,3),2,2))
}



##########################
### ICH/FDA Guidelines ###
##########################

# Load Required Packages
library(data.table)
library(ggplot2)

# Define the Function to Calculate ICH Shlef Life
ich_shelf <- function(DT, byvar="batch", spec.up=360, spec.lo=250, xmax, xmin){
  options(digits=6)
  # setNumericRounding(0)
  #  DT <- DT_pol80
  # byvar="batch"; spec.up=0.65; spec.lo=0.35; xmax=120; xmin=-100
  res <- DT[, {
    ux <- mean(month)
    uy <- mean(resp)
    Syy <- sum((resp - uy) ^ 2)
    Sxx <- sum((month - ux) ^ 2)
    Sxy <- sum((month - ux) * (resp - uy))
    
    slope <- Sxy / Sxx
    intercept <- uy - slope * ux
    n <- length(resp)
    MS <- (Syy - Sxy ^ 2 / Sxx) / (n-2)
    df <- n - 2
    t <- qt(0.975, df)
    sighat <- sqrt((Syy-Sxy^2/Sxx) / df)
    
    m <- MS * t^2
    A <- m / Sxx - slope^2
    B1 <- 2*slope*(spec.up-intercept) - 2*m*ux/Sxx
    B2 <- 2*slope*(spec.lo-intercept) - 2*m*ux/Sxx
    C1 <- m*(1/n+ux^2/Sxx) - (spec.up-intercept)^2
    C2 <- m*(1/n+ux^2/Sxx) - (spec.lo-intercept)^2
    x11 <- (-B1+sqrt(B1^2-4*A*C1))/(2*A)
    x12 <- (-B1-sqrt(B1^2-4*A*C1))/(2*A)
    x21 <- (-B2+sqrt(B2^2-4*A*C2))/(2*A)
    x22 <- (-B2-sqrt(B2^2-4*A*C2))/(2*A)
    
    xu <- ifelse(slope >= 0, min(x11, x12), max(x11, x12)) 
    xl <- ifelse(slope >= 0, max(x21, x22), min(x21, x22))
    if (!is.na(xu) & xu < 0) xu <- 999
    if (!is.na(xl) & xl < 0) xl <- 999
    shelf_life1 <- as.numeric(round(min(xu,xl), 2))
    
    if (sum(is.na(c(x11,x21))) == 0) {
      if (x11*x21 > 0 & intercept <= spec.up & intercept >= spec.lo & x11 > 0) {shelf_life2 <- min(x11, x21)
      }else if (x11*x21 > 0) {shelf_life2 <- 0
      }else shelf_life2 <- max(x11, x21)
    } else shelf_life2 <- NA
    shelf_life2 <- as.numeric(shelf_life2)
    
    list(intercept=intercept, slope=slope, MS=MS, n=n, Sxx=Sxx, xbar=ux, df=df, t=t, sighat=sighat,
         x11=x11, x12=x12, x21=x21, x22=x22, xl=xl, xu=xu, A=A,
         shelf_life1=shelf_life1, shelf_life2=shelf_life2)
  },by=byvar]
  
  # Get non-NA shelflife
  ichSL <- res$shelf_life1[res$shelf_life1 != 0]
  sl1 <- quantile(ichSL, c(0.01, 0.05, 0.5, 0.95), type=2, na.rm=T)
  sl2 <- quantile(res$shelf_life2, c(0.01, 0.05, 0.5, 0.95), type=2, na.rm=T)
  sl3 <- quantile(res$shelf_life2[res$shelf_life2 != 0], c(0.01, 0.05, 0.5, 0.95), type=2, na.rm=T)
  
  if (byvar!="") {
    int_lm <- c(range(res$intercept), mean(res$intercept), sd(res$intercept))
    slope_lm <- c(range(res$slope), mean(res$slope), sd(res$slope))
    param_lm <- rbind(int_lm, slope_lm)
    colnames(param_lm) <- c("lower", "upper", "mean", "sd")
    
    df <- rbind()
    indlist <- c(which.max(res$slope), which.min(res$slope))
    #xmax <- floor(max(res$x12, res$x22)) - 50
    #xmin <- ceiling(min(res$x11, res$x21)) + 50
    for (i in indlist){
      xpt <- seq(xmin, xmax, 0.05)
      ul <- res$intercept[i] + res$slope[i] * xpt + 
        res$sighat[i] * res$t[i] * sqrt(1/res$n[i] + (xpt-res$xbar[i])^2/res$Sxx[i])
      ll <- res$intercept[i] + res$slope[i] * xpt - 
        res$sighat[i] * res$t[i] * sqrt(1/res$n[i] + (xpt-res$xbar[i])^2/res$Sxx[i])
      rl <- res$intercept[i] + res$slope[i] * xpt
      batch <- rep(res$batch[i], length(xpt))
      dftmp <- data.frame(x=xpt, ul=ul, ll=ll, rl=rl, batch=batch)
      df <- rbind(df, dftmp)
    }
    
    res0 <- data.frame(res)
    res1 <- cbind(res0[,1], round(res0[,-1],6))
    colnames(res1)[1] <- "batch"
    
    ymax = max(ul)
    ymin = min(ll)
    p1 <- ggplot(data=df, aes(x=x, group=batch, colour=batch)) +
      geom_line(aes(y=ul), linetype=2) + 
      geom_line(aes(y=ll), linetype=2) + 
      geom_line(aes(y=rl)) + 
      scale_y_continuous(limits = c(ymin, ymax)) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = spec.up) + 
      geom_hline(yintercept = spec.lo)
  }
  
  if (byvar!="") return(list(res=res1, sl1=sl1, sl2=sl2, sl3=sl3, param_lm=param_lm, p1=p1))
  else return(list(res=res))
}
