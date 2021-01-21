nbatch <- c(3, 6, 10, 15, 20, 30, 40)
mintc <- c(30, 13, 6, 3.5, 2.75, 1.8, 1.3)
lo <- loess(mintc ~ nbatch)
pred.lo <- predict(lo, data.frame(nbatch=3:40))

#library(gridExtra)
#library(grid)

#####
tiff("Figure3.tiff", width = 6, height = 4, units = 'in', res = 800)
#layout(matrix(c(1, 1, 1, 2), ncol=4))
par(mar=c(4, 4, 2, 1))
matplot(nbatch, mintc, type="p", pch=16, ylim=c(0, 35), main="Loess Curve Fitted",
        xlab="the number of batches", ylab="tc", cex=0.8,
        panel.first={
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray90")
        abline(v=seq(0, 40, 5), lwd=2 , col="white")
        abline(h=seq(0, 35, by=5),lwd=2 ,  col="white")
})
matplot(3:40, pred.lo, type="l", add=TRUE, pch=4)
# grid.table(matrix(pred.lo[1:5], ncol=1),  show.rownames = TRUE)
dev.off()
#######


lo2 <- loess(mintc ~ nbatch,
                  control = loess.control(surface = "direct"))
pred.lo2 <- predict(lo2, data.frame(nbatch = 3:40), se = TRUE)
temp <- matrix(rep(pred.lo2$fit, each=3)+c(-3, 0, 3)*rep(pred.lo2$se.fit, each=3), nrow=3)
matplot(nbatch, mintc, type="p", pch=16, ylim=c(-20, 45))
matplot(3:40, t(temp), mintc, type="l", add=TRUE)
