# working space
# optim works alright to get weights fit now
# things to do
#  1. make combine work with a series or singularly
#  2. make score work with a series or singularly
#  3. add in a link function step (dist conversion)
#  3. make it so you don't have to tell it what the format is


dists <- draw_predictive_dists(list(mod1, mod2, mod3), 300)
dists4 <- dists[,4,]
obs <- abunds[300 + 1:12]
wts <- c(1/3, 1/3, 1/3)



combined_cdf <- combine(dists1, wts, "cdf")
combined_pdf <- combine(dists1, wts, "pdf")
combined_draw <- combine(dists1, wts, "draw")

score(combined_cdf, obs[1], "crps", "cdf")
score(combined_draw, obs[1], "crps", "draw")
score(combined_pdf, obs[1], "logscore", "pdf")
score(combined_draw, obs[1], "logscore", "draw")

combined_series_cdf <- combine_series(dists, wts, "cdf")
combined_series_pdf <- combine_series(dists, wts, "pdf")
combined_series_draw <- combine_series(dists, wts, "draw")

score_series(combined_series_cdf, obs, "crps", "cdf")
score_series(combined_series_draw, obs, "crps", "draw")
score_series(combined_series_pdf, obs, "logscore", "pdf")
score_series(combined_series_draw, obs, "logscore", "draw")



twts <- alr(c(1/3, 1/3, 1/3))
fit <- optim(twts, ofn, dists = dists, obs = obs, 
             control = list(fnscale = -1))
ialr(fit$par)



############### functions




combine_series <- function(dists, wts = NULL, method = "cdf", nsamps = NULL, 
                           seed = NULL){
  nsteps <- dim(dists)[2]
  out <- vector("list", length = nsteps)
  nmodels <- dim(dists)[3]
  if(is.null(wts)){
    wts <- rep(1/nmodels, nmodels)
  }
  if(is.null(dim(wts))){
    wts <- matrix(wts, nrow = nsteps, ncol = nmodels, byrow = TRUE)
  }
  for(i in 1:nsteps){
    out[[i]] <- combine(dists[,i,], wts[i,], method, nsamps, seed)
  }
  out
}





combine <- function(dists, wts = NULL, method = "cdf", nsamps = NULL, 
                    seed = NULL, digits = 4){
  wts <- enforce_wts(wts, nmods = ncol(dists))
  if(method == "cdf"){
    out <- combine_to_cdf(dists, wts, digits)
  }
  if(method == "pdf"){
    out <- combine_to_pdf(dists, wts, digits)
  }
  if(method == "draw"){
    out <- combine_to_draw(dists, wts, nsamps, seed)
  }
  out
}


combine_to_cdf <- function(dists, wts, digits = 4){
  if(is_int_conf(dists)){
    combine_int_cdfs(dists, wts, digits)
  } else{
    stop("not yet coded for continuous distributions, sorry!")
  }
}

combine_to_pdf <- function(dists, wts, digits = 4){
  if(is_int_conf(dists)){
    combine_int_pmfs(dists, wts, digits)
  } else{
    stop("not yet coded for continuous distributions, sorry!")
  }
}

combine_to_draw <- function(dists, wts, nsamps = NULL, seed = NULL){
  set.seed(seed)
  if(is.null(nsamps)){
    nsamps <- nrow(dists)
  }
  ndists <- ncol(dists)
  dist_choice <- sample(1:ndists, nsamps, TRUE, wts)
  samp_choice <- sample(1:nsamps, nsamps, TRUE)
  draw <- rep(NA, nsamps)
  for(i in 1:nsamps){
    draw[i] <- dists[samp_choice[i], dist_choice[i]]
  }
  draw
}


is_int_conf <- function(x){
  all(x %% 1 == 0)
}


combine_int_cdfs <- function(dists, wts, digits = 4){
  drange <- range(dists, na.rm = TRUE)
  nmods <- ncol(dists)
  x <- seq(drange[1], drange[2], 1)
  y <- rep(0, length(x))
  for(i in 1:nmods){
    y <- y + wts[i] * ecdf(dists[,i])(x)
  }
  y <- round(y, digits)
  incl <- 1:which(y == 1)[1]
  data.frame(x = x[incl], y = y[incl])
}


combine_int_pmfs <- function(dists, wts, digits = 4){
  drange <- range(dists, na.rm = TRUE)
  nmods <- ncol(dists)
  x <- seq(drange[1], drange[2], 1)
  nvals <- length(x)
  y <- rep(0, nvals)
  for(i in 1:nmods){
    nsamps <- nrow(dists)
    tab <- table(dists[,i])
    tabentries <- as.numeric(names(tab))
    ylocations <- match(tabentries, x)
    adder <- wts[i] * tab / nsamps
    y[ylocations] <- y[ylocations] + adder
  }
  cumsumy <- round(cumsum(y), digits)
  incl <- 1:which(cumsumy == 1)[1]
  data.frame(x = x[incl], y = y[incl])
}




ofn <- function(twts, dists, obs){
  wts <- ialr(twts)
  combined_cdf <- combine_series(dists, wts, "cdf")
  score_series(combined_cdf, obs, "crps", "cdf")
}

alr <- function(wts){
  nwts <- length(wts)
  log(wts[1:(nwts - 1)]/wts[nwts])
}

ialr <- function(twts){
  allbut <- exp(twts)/(1 + sum(exp(twts)))
  c(allbut, 1 - sum(allbut))
}


score_series <- function(dists, obs, rule, type = "draw", 
                         orientation = "pos"){
  nsteps <- length(dists)
  scores <- rep(NA, nsteps)
  for(i in 1:nsteps){
    if(!is.na(obs[i])){
      scores[i] <- score(dists[[i]], obs[i], rule, type, orientation)
    }
  }
  sum(scores, na.rm = TRUE)
}




score <- function(dist, obs, rule, type = "draw", orientation = "pos"){
  if(!is_int_conf(obs) | is.vector(dist) && !is_int_conf(dist)){
    stop("not yet coded for continuous distributions, sorry!")
  }
  if(rule == "crps"){
    if(type == "cdf"){
      val <- crps_cdf(dist, obs, orientation)
    } else if (type == "draw"){
      val <- crps_draw(dist, obs, orientation)
    }
  } 
  if(rule == "logscore"){
    if(type == "pdf"){
      val <- logs_pdf(dist, obs, orientation)
    } else if (type == "draw"){
      val <- logs_draw(dist, obs, orientation)
    }
  }
  val
}

    
logs_draw <- function(dist, obs, orientation = "pos"){
  signtoggle <- switch(orientation, "pos" = 1, "neg" = -1)
  log(length(which(dist == obs)) / length(dist)) * signtoggle  
}

logs_pdf <- function(pdf, obs, orientation = "pos"){
  signtoggle <- switch(orientation, "pos" = 1, "neg" = -1)
  log(pdf$y[pdf$x == obs]) * signtoggle
}


crps_cdf <- function(cdf, obs, orientation = "pos"){
  signtoggle <- switch(orientation, "pos" = -1, "neg" = 1)
  k <- seq(min(cdf$x), max(cdf$x), 1)
  nk <- length(k)
  klocations <- match(k, cdf$x)
  cp <- cdf$y[klocations]
  cf <- rep(1, nk)
  cf[which(obs > k)] <- 0
  rpsv <- (cp - cf)^2
  sum(rpsv) * signtoggle
}


crps_draw <- function(draw, obs, orientation = "pos"){
  signtoggle <- switch(orientation, "pos" = -1, "neg" = 1)
  crps_sample(obs, draw) * signtoggle
}



enforce_wts <- function(wts = NULL, nmods = NULL){
  if(is.null(wts)){
    if(is.null(nmods)){
      stop("wts or nmods must be specified")
    }
    wts <- rep(1/nmods, nmods)
  } 
  if(any(wts < 0) || round(sum(wts), 5) != 1){
    stop("wts must be non-negative and sum to 1")
  }
  wts
}