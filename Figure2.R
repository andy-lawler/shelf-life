tc30 <- seq(0.2, 4, by=0.2)
MAE30 <- c(13.4, 10.2, 7.9, 6, 4.5, 3.1, 1.9, 0.9, 0.6, 1, 1.8, 2.3, 2.8, 3.2, 3.8, 4.2, 4.5, 4.9, 5.1, 5.2) 


tc20 <- seq(0, 5, length=21)[-1]
MAE20 <- c(13.85, 11.3, 9, 7.25, 5.8, 4.4, 3.25, 2.3, 1.5, 1.05, 0.74, 0.95, 1.3, 1.75, 2.2, 2.8, 3.3, 3.6, 3.8, 4) 

tc15 <- seq(0, 10, length=21)[-1]
MAE15 <- c(12.4, 9.2, 6.6, 4.6, 2.8, 1.7, 0.92, 0.95, 1.45, 2, 2.55, 3, 3.5, 3.9, 4.2, 4.65, 5, 5.15, 5.5, 5.7)

tc10 <- seq(0, 20, length=21)[-1]
MAE10 <- c(11.7, 7.65, 4.8, 2.8, 1.5, 1.08, 1.45, 1.95, 2.5, 3, 3.4 ,3.8, 4.2, 4.6, 4.8, 5, 5.2, 5.45, 5.55, 5.7)
 

tc6 <- seq(1, length=20, by=3)
MAE6 <- c(13.8, 7, 3.2, 1.7, 1.67, 2.2, 2.55, 2.9, 3.2, 3.4, 3.55, 3.7, 3.95, 4.05, 4.1, 4.2, 4.25, 4.3, 4.35, 4.45)

tc3 <- seq(1, length=20, by=5)
MAE3 <- c(15.35 ,10, 6.3, 4.15, 3.15, 2.7, 2.7 ,2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.45, 3.5, 3.55, 3.55, 3.6, 3.6, 3.65) 

mintc <- rbind(c(30, 2.70), c(13, 1.67), c(6, 1.08), c(3.5, 0.92), c(2.75, 0.74), c(1.8, 0.6))
rownames(mintc) <- c(3, 6, 10, 15, 20, 30)


tiff("Figure2.tiff", width = 8, height = 8, units = 'in', res = 600)
par(mfrow=c(3, 2), mar=c(4, 4, 2, 1))
for(ss in c(3, 6, 10, 15, 20, 30)){
  assign("tc", get(paste0("tc", ss)))
  assign("MAE", get(paste0("MAE", ss)))
  plot(tc, MAE, type="l", ylim=c(0, max(MAE)), xlab="tc", panel.first={
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray90")
    abline(v=axis(1), lwd=2 , col="white")
    abline(h=axis(2),lwd=2 ,  col="white")
  })
  title(paste0("number of batches=", ss))
  points(tc, MAE,pch=16)
  segments(c(0, mintc[paste(ss), 1]), c(mintc[paste(ss), 2], 0), mintc[paste(ss), 1], mintc[paste(ss), 2], lty=3, col=1)
  text(mintc[paste(ss), 1], mintc[paste(ss), 2], paste0("tc=", mintc[paste(ss), 1], ", MAE=", mintc[paste(ss), 2]), adj=c(-0.1, 1), col=1)
}
dev.off()





par(mar=c(0, 0, 0, 0), mfrow=c(1, 1))
plot(0:10, 0:10,pch=16, cex=0.5, type="p", bty="n")
abline(h=seq(0, 10, by=0.5), col="gray50")
abline(h=seq(0, 10, by=1), col=c(2, 1))


par(mar=c(0, 0, 0, 0))
mypng = readPNG('Figure2.png')
plot(0:25, 0:25,pch=16, cex=0.5, type="p", bty="n")
grid.raster(mypng, 0.5, 0.55, width=1)
abline(h=0:25)
