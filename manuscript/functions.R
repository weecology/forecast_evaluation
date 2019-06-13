
#
# figure functions to be pulled out to bbplot are in other scripts


density.default


mass <- function(x, min = NULL, max = NULL){
  xin <- na.omit(x)
  N <- length(xin)
  min <- if_null(min, min(xin))
  max <- if_null(max, max(xin))
  x <- seq.int(min, max)
  nx <- length(x)
  yraw <- rep(NA, nx)
  for(i in 1:nx){
    yraw[i] <- length(xin[xin == x[i]])
  }
  y <- yraw / N
  out <- list(x = x, y = y, n = N)
  class(out) <- c("mass", "list")
  out
}


#
#
make_eval_tab <- function(mod1, mod2, mod3, mod4){
  npitc <- ncol(mod1[[1]]$eval$PIT)
  evalmat <- matrix(NA, nrow = 4 * length(mod1) * 12, ncol = 6 + npitc)
  spots <- 0

  for(i in 1:length(mod1)){
    for(j in 1:4){
      obj <- eval(parse(text = paste0("mod", j)))
      if (i <= length(obj)){
        lead <- obj[[i]]$meta$lead_time
        leads <- 1:lead
        origin <- rep(max(obj[[i]]$meta$in_timeseries), lead)
        mod <- rep(j, lead)
        crps <- obj[[i]]$eval$crps
        logs <- obj[[i]]$eval$logscore
        PIT <- obj[[i]]$eval$PIT
        spots <- max(spots) + 1:lead
        evalmat[spots, 1] <- mod
        evalmat[spots, 2] <- origin
        evalmat[spots, 3] <- origin + leads
        evalmat[spots, 4] <- leads
        evalmat[spots, 5] <- crps
        evalmat[spots, 6] <- logs
        evalmat[spots, 7:(6+npitc)] <- PIT
      }
    }
  }

  colnames(evalmat) <- c("model", "origin", "destin", "lead", "crps", "logs",
                         paste0("PITint", 1:npitc))
  evalmat <- data.frame(evalmat)
  evalmat
}


plot_sp_abunds <- function(plot = NULL, species = NULL){

  abunds_C <- summarize_rodent_data(level = "plot", plots = plot, 
                                    min_plots = 1, clean = FALSE)

  moons <- load_trapping_data(clean = FALSE)$newmoons_table
  abunds <- rep(NA, nrow(moons))
  for(i in 1:nrow(moons)){
    matchrow <- which(abunds_C$period == moons$period[i])
    if (length(matchrow) == 1){
      abunds[i] <- as.numeric(abunds_C[matchrow, species])
    }
  }
  abunds
}

training_ts <- function(starts = 1, ends = 100){
  if(length(starts) == 1){
    starts <- rep(starts, length(ends))
  } else if (length(starts) != length(ends)){
    stop("starts not length 1 or length of ends")
  }
  n_ts <- length(ends)
  in_ts <- vector("list", length = n_ts)
  for(i in 1:n_ts){
    in_ts[[i]] <- starts[i]:ends[i]
  }
  in_ts
}


save_without_models <- function(modobj, file){
  mods <- vector("list", length = length(modobj))
  for(i in 1:length(modobj)){
    mods[[i]] <- list(summary = modobj[[i]]$summary, data = modobj[[i]]$data,
                      meta = modobj[[i]]$meta, eval = modobj[[i]]$eval)
  }
  names(mods) <- names(modobj)
  assign(deparse(substitute(modobj)), mods)
  save(list = deparse(substitute(modobj)), file = file)

}

### function for nonrandomized PIT histogram 
###
### input: 
###   x    observed data 
###   cdf   CDF
###   n_bins number of bins to use

# adapted from Czado, Gneiting and Held, Biometrics

nrPIT <- function(x, cdf, n_bins = 10, ...){
  a.mat <- matrix(0,n_bins,length(x))
  for (i in 1:length(x)){
    Px <- cdf(x[i], ...)
    Px1 <- cdf(x[i] - 1, ...)
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


inits_fun <- function(model){
  rngs <- c("base::Wichmann-Hill", "base::Marsaglia-Multicarry",
            "base::Super-Duper", "base::Mersenne-Twister")
  out <- function(chain = chain){NULL}
  if (model == 1){
    out <- function(chain = chain){
             list(.RNG.name = sample(rngs, 1),
                  .RNG.seed = sample(1:1e+06, 1),
                  mu = rnorm(1, 0, 1),
                  tau.pro = rgamma(1, shape = 0.1, rate = 0.1))
           }
  }
  if (model == 2){
    out <- function(chain = chain){
             list(.RNG.name = sample(rngs, 1),
                  .RNG.seed = sample(1:1e+06, 1),
                  mu = rnorm(1, 0, 1),
                  tau.pro = rgamma(1, shape = 0.1, rate = 0.1),
                  phi = runif(1, -0.95, 0.95))
           }
  }
  if (model == 3){
    out <- function(chain = chain){
             list(.RNG.name = sample(rngs, 1),
                  .RNG.seed = sample(1:1e+06, 1),
                  mu = rnorm(1, 0, 1),
                  tau.pro = rgamma(1, shape = 0.1, rate = 0.1),
                  phi = runif(1, -0.95, 0.95),
                  beta1 = rnorm(1, 0, 2),
                  beta2 = rnorm(1, 0, 2))
           }
  }
  if (model == 4){
    out <- function(chain = chain){
             list(.RNG.name = sample(rngs, 1),
                  .RNG.seed = sample(1:1e+06, 1),
                  mu = rnorm(1, 0, 1),
                  tau.pro = rgamma(1, shape = 0.1, rate = 0.1),
                  phi = runif(1, -0.95, 0.95),
                  theta = runif(1, -0.95, 0.95))
           }
  }
  out
}

monitor_fun <- function(model, meta = NULL, obs = TRUE, obstype = "lead"){
  if (model == 1){
    out <- c("sd.q", "mu")
  }
  if (model == 2){
    out <- c("sd.q", "mu", "phi")
  }
  if (model == 3){
    out <- c("sd.q", "mu", "phi", "beta1", "beta2")
  }
  if (model == 4){
    out <- c("sd.q", "mu", "phi", "theta")
  }
  if (obs){
    if (obstype == "lead"){
      ptimes <- max(meta$in_timeseries) + 1:meta$lead_time
      preds <- paste0("predY[", ptimes, "]")
      out <- c(out, preds)
    }
    if (obstype == "all"){
      out <- c(out, "predY")
    }
  }
  out
}

fod <- function(dates){
  dates <- as.Date(dates)
  yr <- format(dates, "%Y")
  yr0 <- paste0(substr(yr, 1, 3), "0")
  yr9 <- paste0(substr(yr, 1, 3), "9")
  nyd <- as.Date(paste0(yr0, "-01-01"))
  nye <- as.Date(paste0(yr9, "-12-31"))
  round(as.numeric((dates - nyd))/as.numeric((nye - nyd)), 4)
}

data_fun <- function(model, abunds, in_timeseries, lead_time, 
                     moon_dates = NULL){
  if (model == 1){
    Y <- c(abunds[in_timeseries], rep(NA, lead_time))
    meanY <- mean(Y, na.rm = TRUE)
    out <- list(Y = Y, N = length(Y), meanY = meanY)
  }
  if (model == 2){
    Y <- c(abunds[in_timeseries], rep(NA, lead_time))
    meanY <- mean(Y, na.rm = TRUE)
    out <- list(Y = Y, N = length(Y), meanY = meanY)
  }
  if (model == 3){
    all_in <- (min(in_timeseries)):(max(in_timeseries) + lead_time)    
    frde <- fod(moon_dates[all_in])
    Y <- c(abunds[in_timeseries], rep(NA, lead_time))
    meanY <- mean(Y, na.rm = TRUE)
    out <- list(Y = Y, N = length(Y), meanY = meanY,
                sinde = sin(frde), cosde = cos(frde))
  }
  if (model == 4){
    Y <- c(abunds[in_timeseries], rep(NA, lead_time))
    meanY <- mean(Y, na.rm = TRUE)
    out <- list(Y = Y, N = length(Y), meanY = meanY)
  }
  out
}

name_fun <- function(model){
  if (model == 1){
    out <- "model_scripts/mod1.txt"
  }
  if (model == 2){
    out <- "model_scripts/mod2.txt"
  }
  if (model == 3){
    out <- "model_scripts/mod3.txt"
  }
  if (model == 4){
    out <- "model_scripts/mod4.txt"
  }
  out
}


modi <- function(model = NULL, abunds = NULL, moon_dates = NULL,
                 in_timeseries = 1:500, lead_time = 12,
                 n_chains = 3, n_burnin = 5000, n_sample = 10000, n_thin = 1,
                 n_bins = 10, n_adapt = 1000, keeploc = FALSE, quiet = FALSE,
                 save_preds = FALSE){

  run_name <- run_namer(model, in_timeseries, lead_time)
  messageq(run_name, quiet)
  jags_data <- data_fun(model, abunds, in_timeseries, lead_time, moon_dates)
  mod_name <- name_fun(model)

  meta <- list(in_timeseries = in_timeseries, lead_time = lead_time,
               n_chains = n_chains, n_burnin = n_burnin, 
               n_sample = n_sample, n_thin = n_thin, n_bins = n_bins,
               n_adapt = n_adapt, keeploc = keeploc)
  data <- list(abunds = abunds, moon_dates = moon_dates, jags = jags_data)

  mod <- tryCatch(
           run.jags(model = mod_name, monitor = monitor_fun(model, meta),
                    inits = inits_fun(model), data = jags_data, 
                    n.chains = n_chains, adapt = n_adapt, burnin = n_burnin, 
                    sample = n_sample, thin = n_thin, summarise = FALSE, 
                    plots = FALSE, method = "parallel", 
                    keep.jags.files = keeploc),
                  error = function(x){NA})
  eval <- tryCatch(
            eval_mod(model = model, mod = mod, abunds = abunds, 
                          moon_dates = moon_dates,
                          in_timeseries = in_timeseries, 
                          lead_time = lead_time,
                          n_chains = n_chains, n_burnin = n_burnin, 
                          n_sample = n_sample, n_thin = n_thin, 
                          n_bins = n_bins, n_adapt = n_adapt,
                          save_preds = save_preds),
                  error = function(x){NA})
  summary <- summary(mod)
  list(model = mod, summary = summary, data = data, meta = meta, eval = eval)
}


eval_mod <- function(model = NULL, mod = NULL, abunds = NULL, 
                          moon_dates = NULL,
                          in_timeseries = list(1:500), lead_time = 12, 
                          n_chains = 3, n_burnin = 5000, n_sample = 10000, 
                          n_thin = 1, n_bins = 10, n_adapt = 1000, 
                          save_preds = FALSE){
  if (!all(is.na(mod))){
    jags_data <- data_fun(model, abunds, in_timeseries, lead_time, moon_dates)

    ptimes <- max(in_timeseries) + 1:lead_time
    preds <- paste0("predY[", ptimes, "]")

    m <- as.mcmc(combine.mcmc(as.mcmc.list(mod), collapse.chains = TRUE))
    lams <- matrix(NA, nrow = n_chains * n_sample, ncol = lead_time)
    ys <- matrix(NA, nrow = n_chains * n_sample, ncol = lead_time)

    for(i in 1:lead_time){
      matchcol <- which(colnames(m) == preds[i])
      lamst <- m[,matchcol]
      lams[,i] <- lamst
      lamst[which(lamst > 2e9)] <- 2e9
      lamst[which(lamst < 0)] <- 0
      ys[,i] <- rpois(n_chains * n_sample, lamst)
    }

    PIT <- matrix(NA, nrow = lead_time, ncol = n_bins)
    s_r <- rep(NA, lead_time)
    s_l <- rep(NA, lead_time)

    for(i in 1:lead_time){
      time_spot <- max(in_timeseries)+i
      yyy <- abunds[time_spot]
      if(!is.na(yyy)){
        dist <- ys[ , i]
        s_r[i] <- crps_sample(y = yyy, dat = dist)
        s_l[i] <- -log(length(dist[dist==yyy])/length(dist))
        PIT[i, ] <- nrPIT(yyy, ecdf(dist), n_bins)
      }
    }
    preds <- NULL
    if (save_preds){
      preds <- list(y = ys, lambda = lams) 
    }
    out <- list(PIT = PIT, crps = s_r, logscore = s_l, prediction = preds)
  } else{
    out <- list(PIT = NA, crps = NA, logscore = NA, prediction = NULL)
  }
  out
}


mods <- function(model = NULL, abunds = NULL, moon_dates = NULL, 
                 in_timeseries = list(1:500), lead_time = 12, 
                 n_chains = 3, n_burnin = 5000, n_sample = 10000, n_thin = 1,
                 n_bins = 10, n_adapt = 1000, keeploc = FALSE, quiet = FALSE,
                 save_preds = FALSE){
  out <- vector("list", length = length(in_timeseries))
  for(i in 1:length(in_timeseries)){
    out[[i]] <- tryCatch(
                  modi(model = model, abunds = abunds, 
                       moon_dates = moon_dates,
                       in_timeseries = in_timeseries[[i]], 
                       lead_time = lead_time,
                       n_chains = n_chains, n_burnin = n_burnin, 
                       n_sample = n_sample, n_thin = n_thin, n_bins = n_bins,
                       n_adapt = n_adapt, keeploc = keeploc, quiet = quiet,
                       save_preds = save_preds),
                  error = function(x){NA})
    names(out)[i] <- run_namer(model, in_timeseries[[i]], lead_time)
  }
  out 
}

messageq <- function(txt = NULL, quiet = FALSE){
  if (!quiet){
    message(txt)
  }
}

run_namer <- function(model, in_times = 1, lead_time = 1){
  modname <- paste0("model ", model)
  n_in_times <- length(in_times)
  first1 <- in_times[1]
  last1 <- in_times[n_in_times]
  first2 <- in_times[n_in_times] + 1
  last2 <- in_times[n_in_times] + lead_time
  name1 <- paste0("train: ", first1, " to ", last1)
  name2 <- paste0("test: ",  first2, " to ", last2)
  paste0(modname, ": ", name1, "; ", name2)
}




