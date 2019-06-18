stack_origins <- 300:499
n_stack_origins <- length(stack_origins)
stack_fit <- matrix(NA, nrow = n_stack_origins, ncol = 7)

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

  yy1 <- matrix(NA, nrow = nrow(m1), ncol = lead_time)
  yy2 <- matrix(NA, nrow = nrow(m2), ncol = lead_time)
  yy3 <- matrix(NA, nrow = nrow(m3), ncol = lead_time)

  set.seed(123)
  for(i in 1:lead_time){
    yy1[,i] <- rpois(nrow(m1), m1[,paste0("predY[", (ini:fini)[i], "]")])
    yy2[,i] <- rpois(nrow(m2), m2[,paste0("predY[", (ini:fini)[i], "]")])
    yy3[,i] <- rpois(nrow(m3), m3[,paste0("predY[", (ini:fini)[i], "]")])
  }

  nobs <- length(obs)
  nmods <- 3
  probs <- matrix(NA, nrow = nmods, ncol = nobs)
  for(k in 1:nobs){
    probs[1, k] <- length(which(yy1[,k] == obs[k])) / nrow(yy1)
    probs[2, k] <- length(which(yy2[,k] == obs[k])) / nrow(yy2)
    probs[3, k] <- length(which(yy3[,k] == obs[k])) / nrow(yy3)
  }

  llik <- function(probs, fwts){
    nobs <- ncol(probs)
    obsNA <- apply(probs, 2, sum) == 0
    nobsnNA <- sum(1 - obsNA)
    lliks <- rep(0, nobs)
    for(i in 1:nobs){
      lliks[i] <- log(sum(probs[,i] * fwts))
    }
    lliks[obsNA] <- 0
    (1/nobsnNA) * sum(lliks)
  }

  nsteps <- 50
  wts <- matrix(NA, nrow = nsteps + 1, ncol = 3)
  wts[1,] <- 1/3
  delta <- rep(NA, nrow = nsteps + 1)
  delta[1] <- 0
  eps <- 0.001

  obsNA <- is.na(obs)
  nobsnNA <- sum(1 - obsNA) 

  for(i in 1:nsteps){
    for(j in 1:nmods){

      wt <- 0
      for(k in 1:nobs){
        datum <- (wts[i, j] * probs[j, k]) / sum(probs[ , k] * wts[i, ])
        datum[obsNA[k]] <- 0
        wt <- wt + datum
      }
      wts[i + 1, j] <- (1 / nobsnNA) * wt
    }
    lli <- llik(probs, wts[i, ])
    lli1 <- llik(probs, wts[i + 1, ])
    delta[i] <- (lli1 - lli) / abs(lli1)
  }

  stack_fit[l, ] <- c(stack_origin, nobs, nobsnNA, round(wts[nsteps,], 5), 
                      round(delta[nsteps], 8))    
}

colnames(stack_fit) <- c("origin", "possible_N", "actual_N", 
                         "wt_m1", "wt_m2", "wt_m3", "delta")

stack_fit<-data.frame(stack_fit)
plot(stack_fit$origin, stack_fit$wt_m1, ylim=c(0,1), type = "l")
points(stack_fit$origin, stack_fit$wt_m2, ylim=c(0,1), type = "l", col=2)
points(stack_fit$origin, stack_fit$wt_m3, ylim=c(0,1), type = "l", col=3)