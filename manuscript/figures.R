
blank <- function(...){
  plot(1, 1, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", ...)
}


tiff("fig1.tiff", width = 6.5, height = 3.5, units = "in", res = 200)

par(mar = c(2.5, 1.75, 0, 0.5), fig = c(0, 1, 0.01, 0.7))


blank(bty = "L", xlim = c(0,1.25), ylim = c(0,1))
box(bty = "L", lwd = 2)
set.seed(123)
x1 <- seq(0, 0.85, 0.025)
offs <- c(0, -0.005, 0.005)
n <- length(x1)
x2 <- x1 + sample(offs, n, replace = TRUE, prob = c(0.8, 0.1, 0.1))
NAs <- sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.2, 0.8))
x3 <- x2
x3[NAs] <- NA
x3[n] <- x3[n] + 0.02


set.seed(123)
y1 <- 0.35 + 0.175 * sin(x2/0.75 * 2 * pi)
y2 <- y1 + rnorm(n, 0, 0.1)
x4 <- c(x3, 0.95)
n <- length(x4)
y3 <- c(y2, y2[n-1] - 0.025)
points(x4, y3, pch = 1, cex = 0.85, lwd = 2)
set.seed(123)
fx <- c(0.995, 1.08, 1.15, 1.225)
fy1 <- 0.35 + 0.175 * sin(fx/0.75 * 2 * pi)
fy2 <- fy1 + rnorm(fn, 0, 0.1)
fn <- length(fx)



mtext(side = 2, "y ~ G", line = 0.75)
mtext(side = 1, "Time", line = 0.6)
axis(side = 1, at = seq(0, 1.25, 0.05), tck = -0.02, labels = FALSE, lwd = 2)
axis(side = 1, at = seq(0, 1.25, 0.25), tck = -0.04, labels = FALSE, lwd = 2)
axis(side = 2, at = seq(0, 1, 0.25), tck = -0.02, labels = FALSE, lwd = 2)
axis(side = 2, at = seq(0, 1, 0.5), tck = -0.04, labels = FALSE, lwd = 2)


arrows(x4[n-1] - 0.05, y3[n-1] + 0.225, x4[n-1] - 0.005, y3[n-1] + 0.03, 
       length = 0.05, lwd = 2)
text(x4[n-1] - 0.0625, y3[n-1] + 0.275, cex = 0.8, 
      expression(italic("y"[N-1])))
mtext(side = 1, at = x4[n-1], cex = 0.8, line = 0.5, 
      expression(italic("t"[N-1])))
points(rep(x4[n-1], 2), c(-0.1, y3[n-1] - 0.02), type = "l", lty = 2, lwd = 2, 
       xpd = TRUE)
points(rep(x4[n], 2), c(-0.1, y3[n] - 0.02), type = "l", lty = 2, lwd = 2, 
       xpd = TRUE)
mtext(side = 1, at = x4[n], cex = 0.8, line = 0.5, 
      expression(italic("t"[N])))


points(fx, fy2, pch = 1, cex = 0.85, lwd = 2)
#mtext(side = 1, at = fx[1], cex = 0.8, line = 0.5, 
#      expression(italic("t"[N+1])))
mtext(side = 1, at = fx[fn], cex = 0.8, line = 0.5, 
      expression(italic("t"[N+P])))

#points(rep(fx[1], 2), c(-0.1, 0.035), type = "l", lty = 2, lwd = 2, 
#       xpd = TRUE)
points(rep(fx[fn], 2), c(-0.1, fy2[fn] - 0.02), type = "l", lty = 2, lwd = 2, 
       xpd = TRUE)

points(c(x4[n], fx[fn]), rep(-0.26, 2), type = "l", lwd = 2, xpd = TRUE)
points(rep(fx[fn], 2), c(-0.31, -0.2), type = "l", lwd = 2, xpd = TRUE)
points(rep(x4[n], 2), c(-0.31, -0.2), type = "l", lwd = 2, xpd = TRUE)
text(x = mean(c(x4[n], fx[fn])), y = -0.22, "forecast horizon", 
     font = 3, cex = 0.55, xpd = TRUE)
text(x = mean(c(x4[n], fx[fn])), y = -0.29, "or lead time", 
     font = 3, cex = 0.55, xpd = TRUE)

text(x = x4[n] - 0.1375, y = -0.29, "forecast origin", 
     font = 3, cex = 0.55, xpd = TRUE)

arrows(x4[n] - 0.065, -0.28, x4[n] - 0.02, -0.17, length = 0.05, xpd = TRUE,
       lwd = 2)

arrows(x4[1] + 0.02, y3[1] + 0.5, x4[1] + 0.002, y3[1] + 0.03, 
       length = 0.05, lwd = 2)
text(x4[1] + 0.025, y3[1] + 0.55, cex = 0.8, expression(italic("y"[1])))
mtext(side = 1, at = x4[1], cex = 0.8, line = 0.5, expression(italic("t"[1])))
points(rep(x4[1], 2), c(-0.1, y3[1] - 0.02), type = "l", lty = 2, lwd = 2, 
       xpd = TRUE)

u <- par("usr") 
arrows(u[1], u[3], u[2], u[3], length = 0.1, lwd = 2, xpd = TRUE) 
#rect(0.3925, -0.1, 0.48, 1, xpd = NA, col = 0, border = NA)

for(i in 1:fn){
  set.seed(123)
  Gn <- c(rnorm(1e4, fy2[i] - 0.05, 0.1), rnorm(i*9e4, fy2[i] - 0.1, 0.25))
  Gn2 <- Gn[-which(Gn < 0.02)]
  dens <- density(Gn2)
  incl <- which(dens$x < 0.85)

  xxa <- fx[i] + dens$y[incl] * 0.01
  xxb <- fx[i] - dens$y[incl] * 0.01
  xx <- c(xxa, xxb[length(xxb):1])
  xx2 <- c(xx, xx[1])
  yy <- c(dens$x[incl], dens$x[incl][length(dens$x[incl]):1])
  yy2 <- c(yy, yy[1])
  polygon(xx2, yy2, col = rc)

}
points(fx, fy2, pch = 16, cex = 0.85, lwd = 2, col = 0)
points(fx, fy2, pch = 1, cex = 0.85, lwd = 2)


text(fx[1] + 0.04, fy2[1] + 0.45, cex = 0.8, 
      expression(italic("y"[N+1])))
arrows(fx[1] + 0.02, fy2[1] + 0.4, fx[1] + 0.002, fy2[1] + 0.03, 
       length = 0.05, lwd = 2)

par(mar = c(0.5, 1.25, 0.5, 0), fig = c(0.475, 0.675, 0.69, 0.99), 
    new = TRUE)


set.seed(123)
Gn <- c(rnorm(1e4, y3[n-1] - 0.05, 0.1), rnorm(7e3, y3[n-1] - 0.1, 0.25))
Gn2 <- Gn[-which(Gn < 0.05)]
dens <- density(Gn2)
blank(bty = "L", xlim = c(0,1), ylim = c(0,max(dens$y)*1.01))
box(bty = "L", lwd = 2)
points(dens$x, dens$y, type = "l", lwd = 2)
axis(side = 1, at = seq(0, 1.25, 0.25), tck = -0.025, labels = FALSE, lwd = 2)
axis(side = 1, at = seq(0, 1, 0.5), tck = -0.05, labels = FALSE, lwd = 2)
mtext(side = 1, at = 0.4, cex = 0.8, line = 0.25, expression("G"[N-1]))
points(rep(y3[n-1], 2), c(0, max(dens$y)), type = "l", lwd = 2, lty = 3)
points(c(y3[n-1], y3[n-1] + 0.13), c(-0.1, -1.6), type = "l", lwd = 2, 
       lty = 3, xpd = NA)
mtext(side = 2, cex = 0.7, line = 0.25, "Density") 

par(mar = c(0.5, 1.25, 0.5, 0), fig = c(0.05, 0.25, 0.69, 0.99), 
    new = TRUE)

set.seed(123)
Gn <- c(rnorm(1e4, y3[1] - 0.05, 0.1), rnorm(7e3, y3[1] - 0.1, 0.25))
Gn2 <- Gn[-which(Gn < 0.05)]
dens <- density(Gn2)
blank(bty = "L", xlim = c(0,1), ylim = c(0,max(dens$y)*1.01))
box(bty = "L", lwd = 2)
points(dens$x, dens$y, type = "l", lwd = 2)
axis(side = 1, at = seq(0, 1.25, 0.25), tck = -0.025, labels = FALSE, lwd = 2)
axis(side = 1, at = seq(0, 1, 0.5), tck = -0.05, labels = FALSE, lwd = 2)
mtext(side = 1, cex = 0.8, line = 0.25, expression("G"[1]))
points(rep(y3[1], 2), c(0, max(dens$y)), type = "l", lwd = 2, lty = 3)
points(c(y3[1], y2[1] - 0.2), c(-0.1, -1.6), type = "l", lwd = 2, lty = 3, 
       xpd = NA)
mtext(side = 2, cex = 0.7, line = 0.25, "Density") 


par(mar = c(0.5, 1.25, 0.5, 0), fig = c(0.75, 0.9755, 0.69, 0.99), 
    new = TRUE)

set.seed(123)
Gn <- c(rnorm(1e4, fy2[1] - 0.05, 0.1), rnorm(7e3, fy2[1] - 0.1, 0.25))
Gn2 <- Gn[-which(Gn < 0.02)]
dens <- density(Gn2)
Hn <- c(rnorm(1e4, fy2[1] - 0.05, 0.1), rnorm(9e4, fy2[1] - 0.1, 0.25))
Hn2 <- Hn[-which(Hn < 0.02)]
densH <- density(Hn2)

blank(bty = "L", xlim = c(0,1.1), ylim = c(0,max(dens$y)*1.01))
box(bty = "L", lwd = 2)
points(densH$x, densH$y, type = "l", lwd = 2, col = rc)
points(dens$x, dens$y, type = "l", lwd = 2)
axis(side = 1, at = seq(0, 1.25, 0.25), tck = -0.025, labels = FALSE, lwd = 2)
axis(side = 1, at = seq(0, 1, 0.5), tck = -0.05, labels = FALSE, lwd = 2)
mtext(side = 1, cex = 0.8, at = 0.75, line = 0.25, 
      expression("H"[N+1]*", "*"G"[N+1]))
points(rep(fy2[1], 2), c(0, max(dens$y)), type = "l", lwd = 2, lty = 3)
points(c(fy2[1], fy2[1] - 0.4), c(-0.1, -1), type = "l", lwd = 2, 
       lty = 3, xpd = NA)
mtext(side = 2, cex = 0.7, line = 0.25, "Density") 

dev.off()