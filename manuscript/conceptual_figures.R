source("functions.R")

rc <- rgb(0.8, 0.8, 0.8)

# figure 1

tiff("fig1.tiff", width = 6, height = 3.5, units = "in", res = 200)

par(mar = c(2.5, 1.75, 0, 0.75), fig = c(0, 1, 0.01, 0.7))


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
fn <- length(fx)
fy1 <- 0.35 + 0.175 * sin(fx/0.75 * 2 * pi)
fy2 <- fy1 + rnorm(fn, 0, 0.1)



mtext(side = 2, "y ~ G", line = 0.75)
mtext(side = 1, "Time", line = 0.6)
axis(side = 1, at = seq(0, 1.25, 0.05), tck = -0.02, labels = FALSE, lwd = 2)
axis(side = 1, at = seq(0, 1.25, 0.25), tck = -0.04, labels = FALSE, lwd = 2)
axis(side = 2, at = seq(0, 1, 0.25), tck = -0.02, labels = FALSE, lwd = 2)
axis(side = 2, at = seq(0, 1, 0.5), tck = -0.04, labels = FALSE, lwd = 2)


arrows(x4[n-1] - 0.05, y3[n-1] + 0.225, x4[n-1] - 0.005, y3[n-1] + 0.03, 
       length = 0.05, lwd = 2)
text(x4[n-1] - 0.0625, y3[n-1] + 0.275, cex = 0.8, 
      expression(italic("y"[o-1])))
mtext(side = 1, at = x4[n-1], cex = 0.8, line = 0.5, 
      expression(italic("t"[o-1])))
points(rep(x4[n-1], 2), c(-0.1, y3[n-1] - 0.02), type = "l", lty = 2, lwd = 2, 
       xpd = TRUE)
points(rep(x4[n], 2), c(-0.1, y3[n] - 0.02), type = "l", lty = 2, lwd = 2, 
       xpd = TRUE)
mtext(side = 1, at = x4[n], cex = 0.8, line = 0.5, 
      expression(italic("t"[o])))


points(fx, fy2, pch = 1, cex = 0.85, lwd = 2)
#mtext(side = 1, at = fx[1], cex = 0.8, line = 0.5, 
#      expression(italic("t"[o+1])))
mtext(side = 1, at = fx[fn] + 0.03, cex = 0.8, line = 0.5, 
      expression(italic("t"[o+P]*" = "*"t"[N])))

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
      expression(italic("y"[o+1])))
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
mtext(side = 1, at = 0.4, cex = 0.8, line = 0.25, expression("G"[o-1]))

u <- par("usr") 
yatx <- min(dens$y[which(round(dens$x, 2) == round(y3[n-1], 2))])
points(rep(y3[n-1], 2), c(u[3], yatx), type = "l", lwd = 2, lty = 3)
points(c(y3[n-1], y3[n-1] + 0.13), c(u[3], -1.6), type = "l", lwd = 2, 
       lty = 3, xpd = NA)
mtext(side = 2, cex = 0.7, line = 0.25, "Density") 

par(mar = c(0.5, 1.25, 0.5, 0), fig = c(0.05, 0.25, 0.69, 0.99), 
    new = TRUE)

set.seed(123)
Gn <- c(rnorm(1e4, y3[1] + 0.05, 0.1), rnorm(7e3, y3[1] + 0.1, 0.25))
Gn2 <- Gn[-which(Gn < 0.05)]
dens <- density(Gn2)
blank(bty = "L", xlim = c(0,1), ylim = c(0,max(dens$y)*1.01))
box(bty = "L", lwd = 2)
points(dens$x, dens$y, type = "l", lwd = 2)
axis(side = 1, at = seq(0, 1.25, 0.25), tck = -0.025, labels = FALSE, lwd = 2)
axis(side = 1, at = seq(0, 1, 0.5), tck = -0.05, labels = FALSE, lwd = 2)
mtext(side = 1, cex = 0.8, line = 0.25, expression("G"[1]))

u <- par("usr") 
yatx <- min(dens$y[which(round(dens$x, 2) == round(y3[1], 2))])

points(rep(y3[1], 2), c(u[3], yatx), type = "l", lwd = 2, lty = 3)
points(c(y3[1], y2[1] - 0.2), c(u[3], -1.6), type = "l", lwd = 2, lty = 3, 
       xpd = NA)
mtext(side = 2, cex = 0.7, line = 0.25, "Density") 


par(mar = c(0.5, 1.25, 0.5, 0), fig = c(0.75, 0.9755, 0.69, 0.99), 
    new = TRUE)

set.seed(123)
Gn <- c(rnorm(1e4, fy2[1] - 0.05, 0.1), rnorm(9e3, fy2[1] - 0.1, 0.25))
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
      expression("H"[o+1]*", "*"G"[o+1]))
yatx <- min(dens$y[which(round(dens$x, 2) == round(fy2[1], 2))])

u <- par("usr") 
points(rep(fy2[1], 2), c(u[3], yatx), type = "l", lwd = 2, lty = 3)
points(c(fy2[1], fy2[1] - 0.475), c(u[3], -1.1), type = "l", lwd = 2, 
       lty = 3, xpd = NA)
mtext(side = 2, cex = 0.7, line = 0.25, "Density") 

dev.off()


# figure 2

tiff("fig2.tiff", width = 6, height = 4, units = "in", res = 200)

par(mar = c(0, 0, 0, 0), fig = c(0, 1, 0.9, 1))
blank(bty = "n", xlim = c(-1, 20.5), ylim = c(0, 2))

rect(2, 0.7, 2.75, 1.7)
text(x = 3, y = 1.2, "Training datum", adj = 0, cex = 0.6)
rect(5.75, 0.7, 6.5, 1.7, col = rc)
rect(5.75, 0.7, 6.5, 1.7)
text(x = 6.75, y = 1.2, "Test datum", adj = 0, cex = 0.6)
rect(9.25, 0.7, 10, 1.7, lty = 2)
text(x = 10.25, y = 1.2, "Not yet observed datum", adj = 0, cex = 0.6)
rect(14.75, 0.7, 15.5, 1.7, lwd = 2)
text(x = 15.75, y = 1.2, "Forecast origin", adj = 0, cex = 0.6)


par(mar = c(0, 0, 0, 0), fig = c(0, 1, 0.61, 0.91), new = TRUE)
blank(bty = "n", xlim = c(-1, 20.5), ylim = c(0, 2))
text(-1.5, 2.07, "a", cex = 1.125, font = 2, xpd = NA)
text(-1.1, 2.06, "Single origin end-sample testing", cex = 0.7, font = 3, 
     xpd = NA, adj = 0)

lbox <- c(rep(1, 17), rep(2, 3))

trains <- 1:14
tests <- 15:17
oos <- 18:20

rect(-1.5, 1.4, 1, 1.87)
text(x = -0.25, y = 1.635, "Sample", cex = 0.7, xpd = NA)
rect(1, 1.4, max(tests) + 1, 1.87)
rect(max(tests) + 1, 1.4, max(oos) + 1, 1.87, lty = 2)

for(i in 1:max(oos)){
  ltyy <- 0
  if (i < max(tests)){
    ltyy <- 1
  } else if(i < max(oos)){
    ltyy <- 2
  } 
  text(x = i + 0.5, y = 1.635, i, cex = 0.7)
  points(rep(i + 1, 2), c(1.4, 1.87), type = "l", lty = ltyy)
}

rect(-1.5, 0.82, 1, 1.29)
text(x = -0.25, y = 1.133, "Test origin = 14", cex = 0.5, xpd = NA)
text(x = -0.25, y = 0.9763, "True origin = 17", cex = 0.5, xpd = NA)
rect(max(trains) + 1, 0.82, max(tests) + 1, 1.29, lty = 0, col = rc)
rect(1, 0.82, max(tests) + 1, 1.29)
rect(max(tests) + 1, 0.82, max(oos) + 1, 1.29, lty = 2)

for(i in 1:max(oos)){
  ltyy <- 0
  if (i < max(tests)){
    ltyy <- 1
  } else if(i < max(oos)){
    ltyy <- 2
  }
  points(rep(i + 1, 2), c(0.82, 1.29), type = "l", lty = ltyy)
  if (i == max(trains)){
    rect(i, 0.82, i + 1, 1.29, lty = 1, lwd = 2)
    text(x = i + 0.5, y = 1.055, expression("n"["o"[test]]), cex = 0.7)
  }
  if (i == max(tests)){
    rect(i, 0.82, i + 1, 1.29, lty = 1, lwd = 2)
    text(x = i + 0.5, y = 1.055, expression("n"["o"[true]]), cex = 0.7)
  }
}

par(mar = c(0, 0, 0, 0), fig = c(0, 1, 0.0, 0.7), new = TRUE)
blank(bty = "n", xlim = c(-1, 20.5), ylim = c(0, 2))


lbox <- c(rep(1, 17), rep(2, 3))

rect(-1.5, 1.77, 1, 1.97)
text(x = -0.25, y = 1.87, "Sample", cex = 0.7, xpd = NA)
for(i in 1:20){
  rect(i, 1.77, i + 1, 1.97, lty = lbox[i])
  text(x = i + 0.5, y = 1.87, i, cex = 0.7)
}

y10 <- 1.52
y20 <- 1.72
no <- 10:17
p1 <- no + 1
p2 <- no + 2
p3 <- no + 3

for(j in 1:8){
  y1 <- y10 - (j - 1) * 0.2 - (j) * 0.02
  y2 <- y20 - (j - 1) * 0.2 - (j) * 0.02
  y12 <- mean(c(y1, y2))
  rect(-1.5, y1, 1, y2)
  rect(1, y1, j + 10, y2)
  for(i in 1:(j+9)){
    points(c(i, i), c(y1, y2), type = "l")
  }
  if (j < 8){
    otext <- paste0("Test origin = ", 9 + j)
    text(x = -0.25, y = y12,otext, cex = 0.5, xpd = NA)
    t1 <- j + 9
    t2 <- min(max(tests), j + 12)
    rect(t1 + 1, y1, t2 + 1, y2, lty = 0, col = rc)
    rect(t1 + 1, y1, t2 + 1, y2, lty = 1)
    for(i in t1:t2){
      points(c(i, i), c(y1, y2), type = "l")
    }
    rect(j + 9, y1, j + 10, y2, lty = 1, lwd = 2)
    text(x = j + 9.5, y = y12, expression("n"["o"[test]]), cex = 0.7)
  } else{
    otext <- paste0("True origin = ", 9 + j)
    text(x = -0.25, y = y12,otext, cex = 0.5, xpd = NA)
    rect(max(tests) + 1, y1, max(oos) + 1, y2, lty = 2)
    for(i in j:11){
      points(c(i+9, i+9), c(y1, y2), type = "l", lty = 2)
    }
    rect(17, y1, 18, y2, lty = 1, lwd = 2)
    text(x = 17.5, y = y12, expression("n"["o"[true]]), cex = 0.7)
  }
 
}

text(-1.5, 2.07, "b", cex = 1.125, font = 2, xpd = NA)
text(-1.1, 2.06, "Rolling origin end-sample testing", cex = 0.7, font = 3, 
     xpd = NA, adj = 0)

dev.off()


# figure 3 


set.seed(321)

xs <- seq(1, 50, length.out = 30)
lams <- 8 + 0.25 * xs + 3 * sin(2 * pi * xs / 15)
nlams <- length(lams)
x <- rpois(nlams, lambda = lams)
N <- 1e4

tiff("fig3.tiff", width = 6, height = 7, units = "in", res = 200)

topoffset <- 0.94
par(fig= c(0, 0.6, topoffset, 1), mar = c(0, 0, 0, 0))
blank()
text(x = 1, y = 1.15, "Time series of predictive distributions and", 
     cex = 0.8)
text(x = 1, y = 0.825, "true observations", cex = 0.8, font = 1)
par(fig= c(0.6, 0.8, topoffset, 1), new = TRUE)
blank()
text(x = 1, y = 1.15, "Observed (y) vs.", cex = 0.8, font = 1)
text(x = 1, y = 0.825, "Predicted (x)", cex = 0.8, font = 1)
par(fig= c(0.8, 1, topoffset, 1), new = TRUE)
blank()
text(x = 1, y = 1.15, "Probability Integral", cex = 0.8, font = 1)
text(x = 1, y = 0.825, "Transform (PIT)", cex = 0.8, font = 1)

draws1 <- sapply(lams, rpois, n = N)
fig3_row(x, draws1, 1, title = "True generating distribution")

draws2 <- sapply(lams + 2, rpois, n = N)
fig3_row(x, draws2, 2, title = "Positively biased")

draws3 <- sapply(lams - 2, rpois, n = N)
fig3_row(x, draws3, 3, title = "Negatively biased")

draws4 <- sapply(x, rpois, n = N)
fig3_row(x, draws4, 4, title = "Too accurate")



draws5 <- matrix(NA, nrow = N, ncol = nlams)
for(i in 1:nlams){
  draws5[,i] <- round(rnorm(N, lams[i], sd = sqrt(lams[i]) / 1.6))
}
fig3_row(x, draws5, 5, title = "Too precise")

draws6 <- matrix(NA, nrow = N, ncol = nlams)
for(i in 1:nlams){
  draws6[,i] <- rnbinom(N, mu = lams[i], size = 1)
}
fig3_row(x, draws6, 6, title = "Too imprecise")

draws7 <- matrix(NA, nrow = N, ncol = nlams)
for(i in 1:nlams){
  draws7[,i] <-  c(rpois(N/2, max(c(0, lams[i] - 5))), 
                   rpois(N/2, lams[i] + 5))
}
fig3_row(x, draws7, 7, title = "Bimodal with proper median")



dev.off()





