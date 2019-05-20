
library(coda)
library(rjags)
library(R2jags)
library(portalr)
library(forecast)
library(scoringRules)
source("functions.R")

download_observations()


LTCP <- c(4, 11, 14, 17)
nLTCP <- length(LTCP)
abunds_C <- summarize_rodent_data(level = "plot", plots = LTCP, min_plots = 1,
                                  clean = FALSE)
traps <- load_trapping_data(clean = FALSE)
moons <- traps$newmoons_table
abunds_C$newmoon <- NA
for(i in 1:nrow(abunds_C)){
  temp <- moons$newmoonnumber[moons$period == abunds_C$period[i]]
  abunds_C$newmoon[i] <- na.omit(temp)
}

matt <- matrix(NA, nrow = max(abunds_C$newmoon), ncol = nLTCP)
abunds_CDM <- data.frame(newmoon = 1:max(abunds_C$newmoon), matt)
for(i in 1:nLTCP){
  temp <- abunds_C$DM[abunds_C$plot == LTCP[i]]
  nms <- abunds_C$newmoon[abunds_C$plot == LTCP[i]]
  for(j in 1:length(temp)){
    abunds_CDM[nms[j], i + 1] <- temp[j]
  }
  colnames(abunds_CDM)[i + 1] <- paste0("plot ", LTCP[i])
}

CDM <- apply(abunds_CDM[,2:5], 1, mean, na.rm = TRUE)

newmoon <- 1:max(abunds_CDM$newmoon)




xx <- mod1(data = abunds_CDM[,2], in_timeseriesl = list(1:100, 1:101))

  

aa <- auto.arima(Y[1:500])
fc <- forecast(aa, 16)
fsd <- (fc$upper[,1] - fc$lower[,1]) / (2 * qnorm(.5 + fc$level[1] / 200))
fsd <- as.numeric(fsd)
fm <- fc$mean
s_ra<-rep(NA, 16)
s_la<-rep(NA, 16)
for(i in 1:16){
  yyy<-allY[500+i]
  if(!is.na(yyy)){
    s_ra[i]<-crps_norm(y = yyy, fm[i], fsd[i])
    s_la[i]<-logs_norm(y = yyy, fm[i], fsd[i])
  }
}


xx <- as.Date(moons$censusdate)

par(mar = c(3, 4.5, 1, 1), mfrow = c(2,1))
plot(xx, Y, type = "n", ylim = c(-3.1, 17),
     bty = "L", las = 1, xaxt = "n", ylab = "", xlab = "", cex.axis = 1.25)
mtext(side = 2, line = 2.5, "DM Abundance", cex = 2)
axis(2, seq(-3, 17, 1), labels = FALSE, tck = -0.01)
xds <- seq(1980, 2020, 5)
nys <- paste0(xds, "-01-01")
axis(1, at = as.Date(nys), labels = xds)
xds <- seq(1977, 2020, 1)
nys <- paste0(xds, "-01-01")
axis(1, at = as.Date(nys), labels = FALSE, tck = -0.01)


Y2 <- Y[-which(is.na(xx))]
xx2 <- xx[-which(is.na(xx))]
mp2 <- apply(ys, 2, mean)[-which(is.na(xx))]
lp2 <- apply(ys, 2, quantile, probs = 0.025)[-which(is.na(xx))]
up2 <- apply(ys, 2, quantile, probs = 0.975)[-which(is.na(xx))]
xx3 <- xx2[1:(length(xx2) - 1)]

points(xx2, Y2, cex = 1.25, lwd = 2, pch = 1, col = rgb(0.2,0.2,0.2,0.5))
points(xx2, mp2, type = 'l', lwd = 2)
points(xx2, lp2, type = "l", lty = 3)
points(xx2, up2, type = "l", lty = 3)

xx4 <- xx[501:length(xx)]
yy4 <- allY[501:length(xx)]
points(xx4, yy4, cex = 1.25, lwd = 2, pch = 0, col = rgb(0.2,0.2,0.2,1))

points(xx4, fc$mean, type = 'l', lwd = 2, col = 3)
points(xx4, fc$upper[,2], type = "l", lty = 3, col = 3)
points(xx4, fc$lower[,2], type = "l", lty = 3, col = 3)


xxaa <- xx[1:500][-which(is.na(aa$fitted))]
aaf <- aa$fitted[-which(is.na(aa$fitted))]

upper <- aaf + 1.96*sqrt(aa$sigma2)
lower <- aaf - 1.96*sqrt(aa$sigma2)

points(xxaa, aaf, type = "l", lwd = 2, col = 3)
points(xxaa, upper, type = "l", lty = 3, col = 3)
points(xxaa, lower, type = "l", lty = 3, col = 3)




plot(xx, Y, type = "n", ylim = c(-3.1, 17), 
     xlim = as.Date(c("2010-01-01", "2020-01-01")),
     bty = "L", las = 1, xaxt = "n", ylab = "", xlab = "", cex.axis = 1.25)
mtext(side = 2, line = 2.5, "DM Abundance", cex = 2)
axis(2, seq(-3, 17, 1), labels = FALSE, tck = -0.01)
xds <- seq(1980, 2020, 5)
nys <- paste0(xds, "-01-01")
axis(1, at = as.Date(nys), labels = xds)
xds <- seq(1977, 2020, 1)
nys <- paste0(xds, "-01-01")
axis(1, at = as.Date(nys), labels = FALSE, tck = -0.01)

Y2 <- Y[-which(is.na(xx))]
xx2 <- xx[-which(is.na(xx))]
mp2 <- apply(ys, 2, mean)[-which(is.na(xx))]
lp2 <- apply(ys, 2, quantile, probs = 0.025)[-which(is.na(xx))]
up2 <- apply(ys, 2, quantile, probs = 0.975)[-which(is.na(xx))]

points(xx2, Y2, cex = 1.25, lwd = 2, pch = 1, col = rgb(0.2,0.2,0.2,0.5))
points(xx2, mp2, type = 'l', lwd = 2)
points(xx2, lp2, type = "l", lty = 3)
points(xx2, up2, type = "l", lty = 3)

xx4 <- xx[501:length(xx)]
yy4 <- allY[501:length(xx)]
points(xx4, yy4, cex = 1.25, lwd = 2, pch = 0, col = rgb(0.2,0.2,0.2,1))








