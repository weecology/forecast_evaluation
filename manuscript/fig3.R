
# this is the foundation of fig 3



set.seed(1234)
minlam <- 1
maxlam <- 7
nlams <- 1e3
lams <- seq(minlam, maxlam, length.out = nlams)
x <- rpois(nlams, lambda = lams)
plot(x)
n_bins <- 10

par(mfrow=c(6,1), mar=c(4,4,1,1))

# each of these cases really just needs a predictive cdf, so simplify
# and function it!

# also code something to get violin layovers for plots more easily!


# on point

pits1 <- matrix(NA, nrow = nlams, ncol = n_bins)
set.seed(1234)
for(i in 1:nlams){
  lam <- lams[i]
  cdf <- function(x){ppois(x, lambda = lam)}
  
  pits1[i,] <- nrPIT(x[i], cdf, n_bins)
}
spits1 <- apply(pits1, 2, mean)
plot(spits1, type = "h", ylim = c(0, max(spits1) * 1.25))


# upward bias

pits2 <- matrix(NA, nrow = nlams, ncol = n_bins)
set.seed(1234)
for(i in 1:nlams){
  lam <- lams[i]
  cdf <- function(x){ppois(x, lambda = lam + 1)}
  
  pits2[i,] <- nrPIT(x[i], cdf, n_bins)
}
spits2 <- apply(pits2, 2, mean)
plot(spits2, type = "h", ylim = c(0, max(spits2) * 1.25))



# downward bias

pits3 <- matrix(NA, nrow = nlams, ncol = n_bins)
set.seed(1234)
for(i in 1:nlams){
  lam <- lams[i]
  cdf <- function(x){ppois(x, lambda = lam - 1)}
  
  pits3[i,] <- nrPIT(x[i], cdf, n_bins)
}
spits3 <- apply(pits3, 2, mean)
plot(spits3, type = "h", ylim = c(0, max(spits3) * 1.25))


# bimodal


pits4 <- matrix(NA, nrow = nlams, ncol = n_bins)
set.seed(1234)
for(i in 1:nlams){
  lam <- lams[i]
  cdf <- function(x){vals <- c(rpois(5e4, lam-1), rpois(5e4, lam+1))
                     ecdf(vals)(x)}
  
  pits4[i,] <- nrPIT(x[i], cdf, n_bins)
}
spits4 <- apply(pits4, 2, mean)
plot(spits4, type = "h", ylim = c(0, max(spits4) * 1.25))



# too precise

pits5 <- matrix(NA, nrow = nlams, ncol = n_bins)
#set.seed(1234)
for(i in 1:nlams){
  lam <- lams[i]
  cdf <- function(x){pnorm(x, mean = lam, sd = sqrt(lam)/2)}
  pits5[i,] <- nrPIT(x[i], cdf, n_bins)
}
spits5 <- apply(pits5, 2, mean)
plot(spits5, type = "h", ylim = c(0, max(spits5) * 1.25))



# too imprecise

pits6 <- matrix(NA, nrow = nlams, ncol = n_bins)
set.seed(123)
for(i in 1:nlams){
  lam <- lams[i]
  cdf <- function(x){pnbinom(x, mu = lam, size = 1)}
  pits6[i,] <- nrPIT(x[i], cdf, n_bins)
}
spits6 <- apply(pits6, 2, mean)
plot(spits6, type = "h", ylim = c(0, max(spits6) * 1.25))

