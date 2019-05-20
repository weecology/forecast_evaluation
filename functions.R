

run_forecasts <- function(models = NULL,
                          data = NULL,
                          forecasts = NULL,
                          options = forecast_options()){
  nforecasts <- length(forecasts) * length(models)
  out <- vector("list", length = nforecasts)
  a <- rep(seq(1, length(forecasts), 1), length(models))
  b <- rep(seq(1, length(models), 1), each = length(forecasts))
  for (i in 1:nforecasts){
    out[[i]] <- run_forecast(models[a[i]], data, forecasts[b[i]], options)
  }
  out
}

forecast_options <- function(validation_method = ){
  list(validation_method = validation_method, 
}




### function for nonrandomized PIT histogram 
###
### input: 
###   x    observed data 
###   cdf   CDF
###   n_bins number of bins to use

# adapted from Czado, Gneiting and Held, Biometrics

nrPIT <- function(x, cdf, n_bins = 10){
  a.mat <- matrix(0,n_bins,length(x))
  for (i in 1:length(x)){
    Px <- cdf(x[i])
    Px1 <- cdf(x[i] - 1)
    k.vec <- pmax(ceiling(n_bins*Px1),1)
    m.vec <- ceiling(n_bins*Px)
    d.vec <- Px-Px1
    if (k.vec == m.vec){
      a.mat[k.vec,i]=1
    } else { 
      a.mat[k.vec,i] <- ((k.vec/n_bins)-Px1)/d.vec
      if ((k.vec+1)<=(m.vec-1)){
        for (j in ((k.vec+1):(m.vec-1))){ 
          a.mat[j,i] <- (1/(n_bins*d.vec))
        }
      }
      a.mat[m.vec,i] <- (Px-((m.vec-1)/n_bins))/d.vec     
    }
  }
  a <- apply(a.mat, 1, sum)
  (n_bins*a)/(length(x))
}


mod1 <- function(data = NULL, in_timeseriesl = list(1), lead_time = 12, 
                  n_chains = 3, n_burnin = 5000, n_iter = 10000, n_thin = 1,
                  n_bins = 10, quiet = FALSE){
  out <- vector("list", length = length(in_timeseriesl))
  for(i in 1:length(in_timeseriesl)){
    in_times <- in_timeseriesl[[i]]
    n_in_times <- length(in_times)
    name <- iter_namer(in_times, lead_time)
    messageq(paste0("running mod1 with ", name), quiet)
    out[[i]] <- mod1i(data = data, in_timeseries = in_times,
                      lead_time = lead_time, n_chains = n_chains,
                      n_burnin = n_burnin, n_iter = n_iter, n_thin = n_thin,
                      n_bins = n_bins)
    names(out)[i] <- name
  }
  out 
}

messageq <- function(txt = NULL, quiet = FALSE){
  if (!quiet){
    message(txt)
  }
}

iter_namer <- function(in_times = 1, lead_time = 1){
  n_in_times <- length(in_times)
  first1 <- in_times[1]
  last1 <- in_times[n_in_times]
  first2 <- in_times[n_in_times] + 1
  last2 <- in_times[n_in_times] + lead_time
  name1 <- paste0("train: ", first1, " to ", last1)
  name2 <- paste0("test: ",  first2, " to ", last2)
  paste0(name1, "; ", name2)
}

mod1i <- function(data = NULL, in_timeseries = 1, lead_time = 12, 
                  n_chains = 3, n_burnin = 5000, n_iter = 10000, n_thin = 1,
                  n_bins = 10){


  n_keep <- n_iter - n_burnin
  jags_params <- c("sd.q", "mu", "phi", "predY")
  nparam <- 3

  out <- vector("list", length = length(in_timeseries) + lead_time)
  Y <- c(data[in_timeseries], rep(NA, lead_time))
  nY <- length(Y)
  jags_data <- list(Y = Y, N = length(Y))

  mod <- jags(jags_data, parameters.to.save = jags_params, 
              model.file = "mod1.txt", n.chains = n_chains, 
              n.burnin = n_burnin,
              n.thin = n_thin, n.iter = n_iter, DIC = TRUE)  

  lams <- matrix(NA, nrow = n_chains * n_keep, ncol = nY)
  ys <- matrix(NA, nrow = n_chains * n_keep, ncol = nY)

  for(i in 1:nY){
    lams[,i] <- as.vector(mod[[2]]$sims.array[,,i+nparam])
    ys[,i] <- rpois(n_chains * n_keep, lams[,i])
  }
  pred <- list(density = lams, observation = ys)
  meta <- list(in_timeseries = in_timeseries, lead_time = lead_time)


  PIT <- matrix(NA, nrow = lead_time, ncol = n_bins)
  s_r <- rep(NA, lead_time)
  s_l <- rep(NA, lead_time)
  for(i in 1:lead_time){
    time_spot <- max(in_timeseries)+i
    yyy <- data[time_spot]
    if(!is.na(yyy)){
      dist <- ys[ , time_spot]
      s_r[i] <- crps_sample(y = yyy, dat = dist)
      s_l[i] <- -log(length(dist[dist==yyy])/length(dist))
      PIT[i, ] <- nrPIT(yyy, ecdf(dist), n_bins)
    }
  }

  list(model = mod, data = data, jags_data = jags_data, prediction = pred, 
       meta = meta, PIT = PIT, crps = s_r, logscore = s_l)
}




