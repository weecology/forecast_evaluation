# working space

dists <- draw_predictive_dists(list(mod1, mod2, mod3), 300)
obs <- abunds[300 + 1:12]
wts <- c(1/3, 1/3, 1/3)

dists<-dists[,1:2,]
obs<-obs[1:2]

combined_cdf <- combine(dists1, wts, "cdf")
combined_pdf <- combine(dists1, wts, "pdf")
combined_draw <- combine(dists1, wts, "draw")

score(combined_cdf, obs[1], "crps")
score(combined_draw, obs[1], "crps")
score(combined_pdf, obs[1], "logscore")
score(combined_draw, obs[1], "logscore")

combined_series_cdf <- combine(dists, wts, "cdf")
combined_series_pdf <- combine(dists, wts, "pdf")
combined_series_draw <- combine(dists, wts, "draw")

score(combined_series_cdf, obs, "crps")
score(combined_series_draw, obs, "crps")
score(combined_series_pdf, obs, "logscore")
score(combined_series_draw, obs, "logscore")



pars <- c(alr(c(1/3, 1/3, 1/3)), logp1(c(2, 2)))
names(pars) <- c(paste0("twt", 1:2), paste0("beta", 1:2))
fit <- optim(pars, ofn, dists = dists, obs = obs, rule = "crps", 
             linkname = "beta", control = list(fnscale = -1))
ialr(fit$par[1:2])
ilogp1(fit$par[3:4])
fit$val

############### functions

ofn <- function(pars, dists, obs, rule, linkname){
  wts <- ialr(pars[grepl("wt", names(pars))])
  type <- switch(rule, "crps" = "cdf", "logscore" = "pdf")
  combined <- combine(dists, wts, type)
  linked <- link(combined, linkname, pars, raw = dists)
  score(linked, obs, rule)
}


link <- function(dists, linkname = "identity", pars = NULL, raw = NULL){
  if(linkname == "identity"){
    ldists <- dists
  } else if(linkname == "beta"){
    nsteps <- length(dists)
    ldists <- dists
    if (class(dists[[1]])[1] == "cdf"){
      for(i in 1:nsteps){
        ldists[[i]]$y <- pbeta(dists[[i]]$y, 
                               ilogp1(pars["beta1"]), 
                               ilogp1(pars["beta2"]))
      }
    }
    if (class(dists[[1]])[1] == "pdf"){
      wts <- ialr(pars[grepl("wt", names(pars))])
      combined_cdfs <- combine(raw, wts, "cdf")
      for(i in 1:nsteps){
        ldists[[i]]$y <- dists[[i]]$y * 
                         dbeta(combined_cdfs[[i]]$y, 
                               ilogp1(pars["beta1"]), 
                               ilogp1(pars["beta2"]))
      }
    }

  }
  ldists
}

score <- function(dists, obs, rule, orientation = "pos"){
  if(!"list" %in% class(dists)){
    dists <- list(dists)
  }
  nsteps <- length(dists)
  scores <- rep(NA, nsteps)
  for(i in 1:nsteps){
    if(!is.na(obs[i])){
      scores[i] <- score_dist(dists[[i]], obs[i], rule, orientation)
    }
  }
  sum(scores, na.rm = TRUE)
}

score_dist <- function(dist, obs, rule, orientation = "pos"){
  type <- class(dist)[1]
  if(!is_int_conf(obs) | type == "draw" && !is_int_conf(dist)){
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

    
logs_draw <- function(draw, obs, orientation = "pos"){
  signtoggle <- switch(orientation, "pos" = 1, "neg" = -1)
  log(length(which(draw == obs)) / length(draw)) * signtoggle  
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
  crps_sample(obs, as.vector(draw)) * signtoggle
}





combine <- function(dists, wts = NULL, method = "cdf", nsamps = NULL, 
                           seed = NULL){
  ndims <- length(dim(dists))
  if(ndims == 2){
    dists <- array(dists, dim = c(dim(dists)[1], 1, dim(dists)[2]))
  }
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
    out[[i]] <- combine_dists(dists[,i,], wts[i,], method, nsamps, seed)
  }
  if(ndims == 2){
    out[[1]]
  } else{
    out
  }
}





combine_dists <- function(dists, wts = NULL, method = "cdf", nsamps = NULL, 
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
    cdf <- combine_int_cdfs(dists, wts, digits)
  } else{
    stop("not yet coded for continuous distributions, sorry!")
  }
  class(cdf) <- c("cdf", class(cdf))
  cdf
}

combine_to_pdf <- function(dists, wts, digits = 4){
  if(is_int_conf(dists)){
    pdf <- combine_int_pmfs(dists, wts, digits)
  } else{
    stop("not yet coded for continuous distributions, sorry!")
  }
  class(pdf) <- c("pdf", class(pdf))
  pdf
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
  class(draw) <- c("draw", class(draw))
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

alr <- function(wts){
  nwts <- length(wts)
  log(wts[1:(nwts - 1)]/wts[nwts])
}

ialr <- function(twts){
  allbut <- exp(twts)/(1 + sum(exp(twts)))
  c(allbut, 1 - sum(allbut))
}

logp1 <- function(x){
  log(x - 1.0001)
}

ilogp1 <- function(x){
  exp(x) + 1.0001
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