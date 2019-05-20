#################################################################################
# R code for Czado, Gneiting and Held, Biometrics 
#################################################################################

#################################################################################
# Functions and preparations
#################################################################################

### rm(list=ls())
### source("countdata.R")

### function for nonrandomized PIT histogram 
###
### input: 
###   x    observed data 
###   Px   CDF at x 
###   Px1  CDF at x-1 

pit <- function(x, Px, Px1, n.bins=10, y.max=2.75, my.title="PIT Histogram")
  {
  a.mat <- matrix(0,n.bins,length(x))
  k.vec <- pmax(ceiling(n.bins*Px1),1)
  m.vec <- ceiling(n.bins*Px)
  d.vec <- Px-Px1
  for (i in 1:length(x))
      {
      if (k.vec[i]==m.vec[i]) {a.mat[k.vec[i],i]=1}
      else 
        { 
        a.mat[k.vec[i],i]=((k.vec[i]/n.bins)-Px1[i])/d.vec[i]
        if ((k.vec[i]+1)<=(m.vec[i]-1))
           {for (j in ((k.vec[i]+1):(m.vec[i]-1))) {a.mat[j,i]=(1/(n.bins*d.vec[i]))}}
        a.mat[m.vec[i],i]=(Px[i]-((m.vec[i]-1)/n.bins))/d.vec[i]     
        }
      }
  a <- apply(a.mat,1,sum)
  a <- (n.bins*a)/(length(x))
  p <- (0:n.bins)/n.bins
  PIT <- "Probability Integral Transform"
  RF <- "Relative Frequency"
  plot(p, p, ylim=c(0,y.max), type="n", xlab=PIT, ylab=RF, main=my.title) 
  temp1 <- ((1:n.bins)-1)/n.bins
  temp2 <- ((1:n.bins)/n.bins)
  o.vec <- rep(0,n.bins)
  segments(temp1,o.vec,temp1,a)
  segments(temp1,a,temp2,a)
  segments(temp2,o.vec,temp2,a)
  segments(0,0,1,0)
  }
  
### function to perform the Vuong Test
###
### input: 
###   nb    output of glm.nb fit
###   pois  output of glm Poisson fit
###   alpha significance level

vuong.test <- function(nb,pois,alpha=.05)
  {
  mu.pois <- fitted(pois)
  mu.nb <- fitted(nb)
  n <- length(mu.pois)
  l.pois <- rep(0,n)
  l.nb <- rep(0,n)
  theta <- nb$theta
  y <- nb$y
  for(i in 1:n)
    {
    l.pois[i] <- dpois(y[i],mu.pois[i],log=TRUE)
    l.nb[i] <- dnbinom(y[i],size=theta,mu=mu.nb[i],log=TRUE)
    }
  temp1 <- mean(l.pois-l.nb)^2
  temp2 <- mean((l.pois-l.nb)^2)
  w2 <- temp2 - temp1
  lr <- sum(l.pois-l.nb)
  vuong.stat <- (lr/sqrt(w2*n))
  c.alpha <- qnorm(1-(alpha/2))
  if (vuong.stat < (-c.alpha)) {decision <- "NB better than POIS"}
  if (vuong.stat > c.alpha) {decision <- "POIS better than NB"}
  if (abs(vuong.stat)<c.alpha) {decision <- "Vuong test inclusive"}
  out <- list(nb=nb,pois=pois,vuong.stat=vuong.stat,alpha=alpha,c.alpha=c.alpha,decision=decision)
  out
  }

### parameter settings for computing scores

kk <- 100000                            ### cut-off for summations 
my.k <- (0:kk) - 1                      ### to handle ranked probability score

#################################################################################
# Simulation study (Sections 2.4 and 3.4) 
#################################################################################

n <- 200                                ### 200 counts 

lambda <- 5           
theta <- 2                              ### parameterization in R: 
theta.2 <- 1                            ### theta = 1/lambda

set.seed(100)                           ### set random seed 
x <- rnbinom(n,mu=lambda,size=theta)    ### 200 simulated counts from NB(5,1/2) 

### Poisson(5) prediction 

pois.Px <- ppois(x,lambda)                        ### cumulative probabilities
pois.Px1 <- ppois(x-1,lambda)
pois.px <- dpois(x,lambda)                        ### probabilities 

### NB(5,1/2) prediction

nb.Px <- pnbinom(x,mu=lambda,size=theta)          ### cumulative probabilities
nb.Px1 <- pnbinom(x-1,mu=lambda,size=theta)
nb.px <- dnbinom(x,mu=lambda,size=theta)          ### probabilities 

### NB(5,1) prediction

nb.2.Px <- pnbinom(x,mu=lambda,size=theta.2)      ### cumulative probabilities
nb.2.Px1 <- pnbinom(x-1,mu=lambda,size=theta.2)
nb.2.px <- dnbinom(x,mu=lambda,size=theta.2)      ### probabilities 

### plot nonrandomized PIT histograms 
### reproduce Figure 1 

par(mfrow=c(1,3))

pit(x, pois.Px, pois.Px1, n.bins=10, y.max=2.75, my.title="NB(5,0) Prediction")
pit(x, nb.Px, nb.Px1, n.bins=10, y.max=2.75, my.title="NB(5,1/2) Prediction") 
pit(x, nb.2.Px, nb.2.Px1, n.bins=10, y.max=2.75, my.title="NB(5,1) Prediction")  

### compute scores

pois.logs <- - log(pois.px)
  pois.norm <- sum(dpois(my.k,lambda)^2) 
pois.qs <- - 2*pois.px + exp(-2*lambda)*besselI(2*lambda,0)
pois.sphs <- - pois.px / sqrt(exp(-2*lambda)*besselI(2*lambda,0))
  i.cumsum <- cumsum(ppois(my.k,lambda)^2)
  ii.sum <- sum((ppois(my.k,lambda)-1)^2)
  ii.cumsum <- cumsum((ppois(my.k,lambda)-1)^2)
pois.rps <- (i.cumsum[x+1] + ii.sum - ii.cumsum[x+1]) 
pois.dss <- (x-lambda)^2/lambda + log(lambda)
pois.ses <- (x-lambda)^2

nb.logs <- - log(nb.px)
  nb.norm <- sum(dnbinom(my.k,mu=lambda,size=theta)^2) 
nb.qs <- - 2*nb.px + nb.norm
nb.sphs <- - nb.px / sqrt(nb.norm)
  i.cumsum <- cumsum(pnbinom(my.k,mu=lambda,size=theta)^2)
  ii.sum <- sum((pnbinom(my.k,mu=lambda,size=theta)-1)^2)
  ii.cumsum <- cumsum((pnbinom(my.k,mu=lambda,size=theta)-1)^2)
nb.rps <- (i.cumsum[x+1] + ii.sum - ii.cumsum[x+1]) 
nb.dss <- (x-lambda)^2/(lambda*(1+lambda/theta)) + log(lambda*(1+lambda/theta))
nb.ses <- (x-lambda)^2

nb.2.logs <- - log(nb.2.px)
  nb.2.norm <- sum(dnbinom(my.k,mu=lambda,size=theta.2)^2) 
nb.2.qs <- - 2*nb.2.px + nb.2.norm
nb.2.sphs <- - nb.2.px / sqrt(nb.2.norm)
  i.cumsum <- cumsum(pnbinom(my.k,mu=lambda,size=theta.2)^2)
  ii.sum <- sum((pnbinom(my.k,mu=lambda,size=theta.2)-1)^2)
  ii.cumsum <- cumsum((pnbinom(my.k,mu=lambda,size=theta.2)-1)^2)
nb.2.rps <- (i.cumsum[x+1] + ii.sum - ii.cumsum[x+1]) 
nb.2.dss <- (x-lambda)^2/(lambda*(1+lambda/theta.2)) + log(lambda*(1+lambda/theta.2))
nb.2.ses <- (x-lambda)^2

### reproduce Table 1 column by column 

round(c(mean(pois.logs),mean(nb.logs),mean(nb.2.logs)),2)   ### logarithmic score 
round(c(mean(pois.qs),mean(nb.qs),mean(nb.2.qs)),3)         ### quadratic score 
round(c(mean(pois.sphs),mean(nb.sphs),mean(nb.2.sphs)),3)   ### spherical score
round(c(mean(pois.rps),mean(nb.rps),mean(nb.2.rps)),2)      ### ranked probability score 
round(c(mean(pois.dss),mean(nb.dss),mean(nb.2.dss)),2)      ### Dawid-Sebastiani score
round(c(mean(pois.ses),mean(nb.ses),mean(nb.2.ses)),1)      ### squared error score

#################################################################################  
# Count regression for patent data (Section 4) 
#################################################################################

### read patent data
### columns: patent, rdsales, rdw5 

patent.data <- read.table("patentdata.txt", sep="", header=T)

patent <- patent.data[,1]
rdsales <- patent.data[,2]
rdw5 <- patent.data[,3]

n <- length(patent)

### count regressions in cross validation mode 

### Poisson regression 

p.pois.lambda <- rep(0,n)
p.pois.Px <- rep(0,n)
p.pois.Px1 <- rep(0,n) 
p.pois.px <- rep(0,n) 

for (i in 1:n)
  {
  temp <- glm(patent[-i] ~ rdsales[-i] + rdw5[-i], family=poisson)
  beta <- coef(temp)
  p.pois.lambda[i] <- exp(t(beta)%*%c(1,rdsales[i],rdw5[i]))
  p.pois.Px[i] <- ppois(patent[i],p.pois.lambda[i])
  p.pois.Px1[i] <- ppois(patent[i]-1,p.pois.lambda[i])
  p.pois.px[i] <- dpois(patent[i],p.pois.lambda[i])
  }

### negative binomial regression 

p.nb.lambda <- rep(0,n)
p.nb.theta <- rep(0,n)
p.nb.Px <- rep(0,n)
p.nb.Px1 <- rep(0,n) 
p.nb.px <- rep(0,n) 
 
library(MASS)    

for (i in 1:n)
  {
  temp <- glm.nb(patent[-i] ~ rdsales[-i] + rdw5[-i])
  beta <- coef(temp)
  p.nb.lambda[i] <- exp(t(beta)%*%c(1,rdsales[i],rdw5[i]))
  p.nb.theta[i] <- temp$theta
  p.nb.Px[i] <- pnbinom(patent[i],size=p.nb.theta[i],mu=p.nb.lambda[i])
  p.nb.Px1[i] <- pnbinom(patent[i]-1,size=p.nb.theta[i],mu=p.nb.lambda[i])
  p.nb.px[i] <- dnbinom(patent[i],size=p.nb.theta[i],mu=p.nb.lambda[i])
  }

### plot nonrandomized PIT histograms 
### reproduce top panel in Figure 2 

par(mfrow=c(1,2))

pit(patent, p.pois.Px, p.pois.Px1, y.max=3.5, my.title="Poisson Regression")
pit(patent, p.nb.Px, p.nb.Px1, y.max=3.5, my.title="Negative Binomial Regression")

### auxiliary calculations for marginal calibration diagram 

x <- patent 
breaks <- c(5,10,25,50,100,500)

x.uniq <- 0:breaks[1]
n.uniq <- length(x.uniq)
temp <- table(x)[1:(breaks[1]+1)]
f.pois <- matrix(0,n,n.uniq)
f.nb <- matrix(0,n,n.uniq)
s <- temp
for (i in 1:n)
  {
  f.pois[i,] <- dpois(x.uniq, p.pois.lambda[i])
  f.nb[i,] <- dnbinom(x.uniq, mu=p.nb.lambda[i], size=p.nb.theta[i])
  }
fm.pois <- apply(f.pois,2,mean)
fm.nb <- apply(f.nb,2,mean)
tab <- rbind(s,fm.pois,fm.nb)
dimnames(tab)[[1]] <- c("Obs","Pois","NB")
dimnames(tab)[[2]] <- names(temp)

x.group <- cut(x[x>breaks[1]],breaks=breaks)
s.group <- table(x.group)
n.group <- length(breaks)-1
f.pois.group <- matrix(0,n,n.group)
f.nb.group <- matrix(0,n,n.group)
k <- length(breaks)
b.low <- breaks[1:k-1]
b.up <- breaks[2:k]
for (i in 1:n)
  {
  f.pois.group[i,] <- ppois(b.up, p.pois.lambda[i])-ppois(b.low,p.pois.lambda[i])
  f.nb.group[i,] <- pnbinom(b.up, mu=p.nb.lambda[i], size=p.nb.theta[i])-pnbinom(b.low,mu=p.nb.lambda[i], size=p.nb.theta[i])
    }
fm.pois.group <- apply(f.pois.group,2,mean)
fm.nb.group <- apply(f.nb.group,2,mean)

### table with observed and predicted frequencies (numbers of such counts) 

tab.total <- rbind(c(s,s.group),n*c(fm.pois,fm.pois.group),n*c(fm.nb,fm.nb.group))
dimnames(tab.total)[[1]] <- dimnames(tab)[[1]]

round(tab.total,1)

### plot marginal calibration diagram 
### reproduce lower panel in Figure 2 

par(mfrow=c(1,1))

plot(c(s,s.group), axes=F, xlab="Count", ylab="Frequency", type="b", pch="+", ylim=c(0,n*.22))
points(c(n*fm.pois,n*fm.pois.group), type="b", lty=2)
points(c(n*fm.nb,n*fm.nb.group), type="b", lty=3)
legend(2.5, n*.21, legend=c("Observations","Poisson","Negative Binomial"), lty=c(1,2,3))
axis(1, at=1:length(n*c(fm.nb,fm.nb.group)), labels=dimnames(tab.total)[[2]])
axis(2)

### compute scores

p.pois.logs <- - log(p.pois.px) 
  p.pois.norm <- 1:n
  for (i in 1:n) {p.pois.norm[i] <- sum(dpois(my.k,p.pois.lambda[i])^2)} 
p.pois.qs <- - 2*p.pois.px + p.pois.norm
p.pois.sphs <- - p.pois.px / sqrt(p.pois.norm)
p.pois.rps <- 1:n 
  for (i in 1:n) 
    {p.pois.rps[i] <- sum(ppois((-1):(patent[i]-1),p.pois.lambda[i])^2) + sum((ppois(patent[i]:kk,p.pois.lambda[i])-1)^2)}
p.pois.dss <- (patent-p.pois.lambda)^2/p.pois.lambda + log(p.pois.lambda)
p.pois.ses <- (patent-p.pois.lambda)^2

p.nb.px <- dnbinom(patent,mu=p.nb.lambda,size=p.nb.theta)
p.nb.logs <- - log(p.nb.px)
  p.nb.norm <- 1:n
  for (i in 1:n) 
    {p.nb.norm[i] <- sum(dnbinom(my.k,mu=p.nb.lambda[i],size=p.nb.theta[i])^2)} 
p.nb.qs <- - 2*p.nb.px + p.nb.norm
p.nb.sphs <- - p.nb.px / sqrt(p.nb.norm)
p.nb.rps <- 1:n 
  for (i in 1:n) 
    {
    p.nb.rps[i] <- sum(pnbinom((-1):(patent[i]-1),mu=p.nb.lambda[i],size=p.nb.theta[i])^2) 
    p.nb.rps[i] <- p.nb.rps[i] + sum((pnbinom(patent[i]:kk,mu=p.nb.lambda[i],size=p.nb.theta[i])-1)^2)
    }
p.nb.dss <- (patent-p.nb.lambda)^2/(p.nb.lambda*(1+p.nb.lambda/p.nb.theta)) + log(p.nb.lambda*(1+p.nb.lambda/p.nb.theta))
p.nb.ses <- (patent-p.nb.lambda)^2

### reproduce Table 2 column by column 

round(c(mean(p.pois.logs),mean(p.nb.logs)),2)    ### logarithmic score
round(c(mean(p.pois.qs),mean(p.nb.qs)),2)        ### quadratic score
round(c(mean(p.pois.sphs),mean(p.nb.sphs)),2)    ### spherical score  
round(c(mean(p.pois.rps),mean(p.nb.rps)),1)      ### ranked probability score
round(c(mean(p.pois.dss),mean(p.nb.dss)),1)      ### Dawid-Sebastiani score
round(c(mean(p.pois.ses),mean(p.nb.ses)),0)      ### squared error score   

### AIC and Vuong test  

attach(patent.data)
library(MASS)

patent.pois <- glm(patent ~ rdsales + rdw5, family=poisson)
patent.nb <- glm.nb(patent~rdsales+rdw5)
round(patent.pois$aic,1)
round(patent.nb$aic,1)

vuong.test(patent.nb,patent.pois,alpha=.01)
