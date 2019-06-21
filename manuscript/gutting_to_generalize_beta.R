
# will become a function
# prep_predictive_dists <- function(){
#
  stack_origin <- 300
  list_spot <- stack_origin - 299
  max_lead_time <- 12
  lead_time <- length(which((stack_origin + 1:max_lead_time) <= 500))

  ini <- length(mod1[[list_spot]]$meta$in_timeseries) + 1
  fini <- ini + lead_time
  obs <- abunds[stack_origin + 1:lead_time]

  m1 <- as.mcmc(combine.mcmc(as.mcmc.list(mod1[[list_spot]]$model), 
                            collapse.chains = TRUE))
  m2 <- as.mcmc(combine.mcmc(as.mcmc.list(mod2[[list_spot]]$model), 
                            collapse.chains = TRUE))
  m3 <- as.mcmc(combine.mcmc(as.mcmc.list(mod3[[list_spot]]$model), 
                            collapse.chains = TRUE))

  yyy <- array(NA, dim = c(nrow(m3), 3, lead_time))

  set.seed(123)
  for(i in 1:lead_time){
    yyy[,1,i] <- rpois(nrow(m1), m1[,paste0("predY[", (ini:fini)[i], "]")])
    yyy[,2,i] <- rpois(nrow(m2), m2[,paste0("predY[", (ini:fini)[i], "]")])
    yyy[,3,i] <- rpois(nrow(m3), m3[,paste0("predY[", (ini:fini)[i], "]")])
  }
  yyy



transform_dist <- function(dist, link_fn = "identity", ...){
  do.call(link_fn, list(dist, ...))
}

score_dist <- function(dist, obs, score_fn = "logscore", ...){
  do.call(score_fn, list(dist = dist, obs = obs, ...))
}


combine_dists

estimate_stack_weights <- function(dists, obs, score_fn = "logscore",
                                   link_fn = "identity"){



}


















stack_origins <- 300:498
n_stack_origins <- length(stack_origins)
stack_fit_b <- matrix(NA, nrow = n_stack_origins, ncol = 10)

  llik <- function(probs, fwts, b1, b2){
    nobs <- ncol(probs)
    obsNA <- apply(probs, 2, sum) == 0
    nobsnNA <- sum(1 - obsNA)
    lliks <- rep(0, nobs)
    for(i in 1:nobs){
      lliks[i] <- dbeta(sum(probs[,i] * fwts), b1[1], 1, log = TRUE)
    }
    lliks[obsNA] <- 0
    (1/nobsnNA) * sum(lliks)
  }


    fn <- function(x, wts, probs){
      liks <- dbeta((wts %*% probs), x[1], x[2], log = TRUE)
      liks[liks == Inf] <- 0
      -sum(liks) 
    }


  nobs <- length(obs)
  nmods <- 3
  probs <- matrix(NA, nrow = nmods, ncol = nobs)
  for(k in 1:nobs){
    probs[1, k] <- length(which(yy1[,k] == obs[k])) / nrow(yy1)
    probs[2, k] <- length(which(yy2[,k] == obs[k])) / nrow(yy2)
    probs[3, k] <- length(which(yy3[,k] == obs[k])) / nrow(yy3)
  }



for(l in 1:n_stack_origins){
  stack_origin <- stack_origins[l]
  list_spot <- stack_origin - 299
  max_lead_time <- 12
  lead_time <- length(which((stack_origin + 1:max_lead_time) <= 500))

  ini <- length(mod1[[list_spot]]$meta$in_timeseries) + 1
  fini <- ini + lead_time
  obs <- abunds[stack_origin + 1:lead_time]

  m1 <- as.mcmc(combine.mcmc(as.mcmc.list(mod1[[list_spot]]$model), 
                            collapse.chains = TRUE))
  m2 <- as.mcmc(combine.mcmc(as.mcmc.list(mod2[[list_spot]]$model), 
                            collapse.chains = TRUE))
  m3 <- as.mcmc(combine.mcmc(as.mcmc.list(mod3[[list_spot]]$model), 
                            collapse.chains = TRUE))

  yyy <- array(NA, dim = c(nsamples, nmods, lead_time))

  set.seed(123)
  for(i in 1:lead_time){
    yyy[,1,i] <- rpois(nrow(m1), m1[,paste0("predY[", (ini:fini)[i], "]")])
    yyy[,2,i] <- rpois(nrow(m2), m2[,paste0("predY[", (ini:fini)[i], "]")])
    yyy[,3,i] <- rpois(nrow(m3), m3[,paste0("predY[", (ini:fini)[i], "]")])
  }

ppars<-c(wt1=0.001, wt2=0.001, wt3=0.998)
  sfit <- constrOptim(ppars, stack_fn, grad = NULL, ui = ui, 
                      ci = ci, obs = obs, dists = yyy)
  sfit$par

ui %*% ppars - ci


cds<-rep(NA, 12)
for(i in 1:12){  
combs <- combine_dists(yyy[,,i], sfit$par)
cds[i] <- score_dist(combs, obs[i], orientation ="pos")
}

cds2<-rep(NA, 12)
for(i in 1:12){  
combs <- combine_dists(yyy[,,i], stack_fit[1,4:6])
cds2[i] <- score_dist(combs, obs[i], orientation ="pos")
}

sum(cds)
sum(cds2)

dists <- yyy
dim(dists)

combine_dists <- function(dists, wts){
  wts[wts<0] <- 0
  wts <- wts/sum(wts)
  ndists <- length(wts)
  dist <- 1:ndists    
  nsamples <- nrow(dists)
  dist_samps <- sample(dist, nsamples, TRUE, wts)
  samps <- NULL
  for(i in 1:ndists){
    nsamplesi <- length(which(dist_samps == i))
    dist_sampsi <- sample(dists[,i], nsamplesi)
    samps <- c(samps, dist_sampsi)
  }
  samps
}



logscore <- function(dist, obs, orientation = "pos", thresh_digits = 10){
  ssign <- switch(orientation, "pos" = 1, "neg" = -1)
  dist <- round(dist, thresh_digits)
  obs <- round(obs, thresh_digits)
  ssign * length(which(dist == obs)) / length(dist)
}








stack_fn <- function(pars, obs, dists, orientation = "pos", seed = 123){
  set.seed(seed) # allows for verification of stable solution
  ssign <- switch(orientation, "pos" = -1, "neg" = 1)
  wts <- pars[grepl("wt", names(pars))]
  dd <- length(dim(dists))
  if (dd == 2){
    dists2 <- array(NA, dim = c(dim(dists), 1))
    dists2[,,1] <- dists
    dists <- dists2
  }
  lead_time <- dim(dists)[3]
  scores <- rep(NA, lead_time)
  for(i in 1:lead_time){ 
    dist <- combine_dists(dists[,,i], wts)
    t_dist <- transform_dist(dist)
    scores[i] <- score_dist(t_dist, obs)
  }
  sum(ssign * scores)
}


ppars[grepl("wt", names(ppars))]


ppars<-c("wt1" = 0, "wt2" = 0, "wt3" = 1)

  optsign <- switch(orientation, "pos" = -1, "neg" = 1)

ui <- diag(nmods)
ci <- rep(0, nmods)
constrOptim(ppars, stack_fn, grad = NULL, ui = ui, 
            ci = ci, obs = obs, dists = yyy)

stack_fn(c("wt1" = 0.0205, "wt2" = 0.5626, "wt3" = 0.505), obs, yyy)






obsX <- obs
obsX[which(is.na(obs))] <- 0
crpsv <- -crps_sample(y = obsX, dat = cc)
crpsv[which(is.na(obs))] <- NA



  optim(pars, scorefun, obs)

  wts <- rep(1/3, 3)
  scorefun <- function(){
    liks <- 
    liks[liks == Inf] <- 0
    -sum(liks)     
  }









  nsteps <- 50
  wts <- matrix(NA, nrow = nsteps + 1, ncol = 3)
  wts[1,] <- 1/3
  delta <- rep(NA, nrow = nsteps + 1)
  delta[1] <- 0
  eps <- 0.001
  b1 <- rep(NA, length = nsteps + 1)
  b1[1] <- 1
  b2 <- rep(NA, length = nsteps + 1)
  b2[1] <- 1
  loglik <- rep(NA, length(nsteps))
  obsNA <- is.na(obs)
  nobsnNA <- sum(1 - obsNA) 

  for(i in 1:nsteps){

    # optimize b given wts



    bbs <- optim(c(b1[i], b2[i]), fn, wts = wts[i,], probs = probs)    
    b1[i + 1] <- bbs$par[1]
    b2[i + 1] <- bbs$par[2]


    for(j in 1:nmods){

      # optimize wts given bs
      wt <- 0
      for(k in 1:nobs){
        datum <- pbeta((wts[i, j] * probs[j, k]), b1[i + 1], b2[i + 1]) / 
                 sum(pbeta((wts[i, ] * probs[, k]), b1[i + 1], b2[i + 1]))
        datum[obsNA[k]] <- 0
        wt <- wt + datum
      }
      wts[i + 1, j] <- (1 / nobsnNA) * wt
    }
    lli <- llik(probs, wts[i, ], b1[i], b2[i])
    lli1 <- llik(probs, wts[i + 1, ], b1[i], b2[i])
    loglik[i] <- lli1
    delta[i] <- (lli1 - lli) / abs(lli1)
  }
  
  stack_fit_b[l, ] <- c(stack_origin, nobs, nobsnNA, round(wts[nsteps,], 5), 
                      delta[nsteps], lli1,
                      beta1 = b1[nsteps], beta2 = b2[nsteps])    
}

colnames(stack_fit_b) <- c("origin", "possible_N", "actual_N", 
                         "wt_m1", "wt_m2", "wt_m3", "delta", "llik", 
                         "b1", "b2")

stack_fit_b<-data.frame(stack_fit_b)
plot(stack_fit_b$origin, stack_fit_b$wt_m1, ylim=c(0,1), type = "l", lty = 3)
points(stack_fit_b$origin, stack_fit_b$wt_m2, type = "l", col=2, lty = 3)
points(stack_fit_b$origin, stack_fit_b$wt_m3, type = "l", col=3, lty = 3)





gg<-matrix(NA, nrow = 10, ncol = 3)
for(i in 1:10){
gg[i,1] <- runif(1, 0.1, 0.3)
gg[i,2] <- runif(1, 0.1, 0.3)
gg[i,3] <- 1 - gg[i,1] - gg[i,2]
}
gg <- as.matrix(stack_fit_b[,4:6])
tt <- stack_fit_b[,1]
sm <- multinom(gg~tt, weights=stack_fit_b[,3])

points(tt,sm$fitted.values[,1],type="l", lwd=2,lty=2)
points(tt,sm$fitted.values[,2],type="l", lwd=2,lty=2,col=2)
points(tt,sm$fitted.values[,3],type="l", lwd=2,lty=2,col=3)
