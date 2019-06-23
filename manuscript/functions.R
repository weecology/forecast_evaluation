
load_stack_df <- function(){
  load("model_output/stack_df.RData", envir = parent.frame(n = 2))
}

stack_mods <- function(mods, test_origins, abunds, save = TRUE){
  ntos <- length(test_origins)
  pars1 <- c(alr(c(1/3, 1/3, 1/3)))
  names(pars1) <- c(paste0("twt", 1:2))
  pars2 <- c(alr(c(1/3, 1/3, 1/3)), logp1(2))
  names(pars2) <- c(paste0("twt", 1:2), "beta")
  pars3 <- c(alr(c(1/3, 1/3, 1/3)), logp1(c(2, 2)))
  names(pars3) <- c(paste0("twt", 1:2), paste0("beta", 1:2))


  stack_tab <- matrix(NA, nrow = ntos * 6, ncol = 9)
  rid <- 1
  for(i in 1:ntos){

    test_origin <- test_origins[i]
    print("#######################")
    print(paste0("test origin ", test_origin))
    dists <- draw_predictive_dists(mods, test_origin, last_in = 500)
    obs <- abunds[test_origin + 1:12]

    print("1")
    fit1 <- tryCatch(
              optim(pars1, ofn, dists = dists, obs = obs, rule = "crps", 
                    linkname = "identity", control = list(fnscale = -1)),
              error = function(x){list()})
    print("2")
    fit2 <- tryCatch(
               optim(pars2, ofn, dists = dists, obs = obs, rule = "crps", 
                     linkname = "beta1", control = list(fnscale = -1)),
              error = function(x){list()})
    print("3")
    fit3 <- tryCatch(
               optim(pars3, ofn, dists = dists, obs = obs, rule = "crps", 
                     linkname = "beta2", control = list(fnscale = -1)),
              error = function(x){list()})

    print("4")
    fit4 <- tryCatch(
              optim(pars1, ofn, dists = dists, obs = obs, rule = "logscore", 
                    linkname = "identity", control = list(fnscale = -1)),
              error = function(x){list()})
    print("5")
    fit5 <- tryCatch(
               optim(pars2, ofn, dists = dists, obs = obs, rule = "logscore", 
                     linkname = "beta1", control = list(fnscale = -1)),
              error = function(x){list()})
    print("6")
    fit6 <- tryCatch(
               optim(pars3, ofn, dists = dists, obs = obs, rule = "logscore", 
                     linkname = "beta2", control = list(fnscale = -1)),
              error = function(x){list()})

    stack_tab[rid + 0, ] <- c(test_origin, "crps", "id", 
                              ialr(fit1$par), NA, NA, 
                              fit1$val)
    stack_tab[rid + 1, ] <- c(test_origin, "crps", "beta1", 
                              ialr(fit2$par[1:2]), ilogp1(fit2$par[3]), NA,  
                              fit2$val)
    stack_tab[rid + 2, ] <- c(test_origin, "crps", "beta2", 
                              ialr(fit3$par[1:2]), ilogp1(fit3$par[3:4]), 
                              fit3$val)
    stack_tab[rid + 3, ] <- c(test_origin, "logscore", "id", 
                              ialr(fit4$par), NA, NA, 
                              fit4$val)
    stack_tab[rid + 4, ] <- c(test_origin, "logscore", "beta1", 
                              ialr(fit5$par[1:2]), ilogp1(fit5$par[3]), NA,  
                              fit5$val)
    stack_tab[rid + 5, ] <- c(test_origin, "logscore", "beta2", 
                              ialr(fit6$par[1:2]), ilogp1(fit6$par[3:4]), 
                              fit6$val)
    rid <- rid + 6
  }

  stack_df <- data.frame(stack_tab, stringsAsFactors = FALSE)
  colnames(stack_df)<-c("origin", "rule", "link", paste0("wt", 1:3), 
                        "b1", "b2", "val")

  for(i in c(1,4:9)){
    stack_df[,i] <- as.numeric(stack_df[,i])
  }

  if(save){
    save(stack_df, file = "model_output/stack_df.RData")
  }
  stack_df
}

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
  } else if(linkname == "beta1"){
    nsteps <- length(dists)
    ldists <- dists
    if (class(dists[[1]])[1] == "cdf"){
      for(i in 1:nsteps){
        ldists[[i]]$y <- pbeta(dists[[i]]$y, 
                               ilogp1(pars["beta"]), 
                               ilogp1(pars["beta"]))
      }
    }
    if (class(dists[[1]])[1] == "pdf"){
      wts <- ialr(pars[grepl("wt", names(pars))])
      combined_cdfs <- combine(raw, wts, "cdf")
      for(i in 1:nsteps){
        ldists[[i]]$y <- dists[[i]]$y * 
                         dbeta(combined_cdfs[[i]]$y, 
                               ilogp1(pars["beta"]), 
                               ilogp1(pars["beta"]))
      }
    }
  } else if(linkname == "beta2"){
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
  yy <- sum(pdf$y)
  log(pdf$y[pdf$x == obs] / yy) * signtoggle
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

draw_predictive_dists <- function(models, origin, list_spot = origin-299, 
                                 max_lead_time = 12, last_in = Inf,
                                 seed = NULL, n = NULL){

  nmodels <- length(models)
  mds <- vector("list", length = nmodels)
  nsamps <- rep(NA, length = nmodels)
  for(i in 1:nmodels){
    mds[[i]] <- as.mcmc(
                combine.mcmc(
                as.mcmc.list(
                  models[[1]][[list_spot]]$model), collapse.chains = TRUE))
    nsamps[i] <- nrow(mds[[i]])
  }
  if(is.null(n)){
    n <- max(nsamps)
  }
  lead_time <- length(which((origin + 1:max_lead_time) <= last_in))
  yyy <- array(NA, dim = c(n, lead_time, nmodels))
  for(i in 1:nmodels){
    yyy[,,i] <- draw_predictive_dist(models[[i]], origin, list_spot,
                                     max_lead_time, last_in, seed, n)
  }
  yyy
}


draw_predictive_dist <- function(model, origin, list_spot = origin-299, 
                                 max_lead_time = 12, last_in = Inf,
                                 seed = NULL, n = NULL){
  
  lead_time <- length(which((origin + 1:max_lead_time) <= last_in))

  ini <- length(mod1[[list_spot]]$meta$in_timeseries) + 1
  fini <- ini + lead_time
  md <- as.mcmc(combine.mcmc(as.mcmc.list(model[[list_spot]]$model), 
                            collapse.chains = TRUE))

  if(is.null(n)){
    n <- nrow(md)
  }
  yyy <- matrix(NA, nrow = n, ncol = lead_time)

  set.seed(seed)
  for(i in 1:lead_time){
    lams <-  md[,paste0("predY[", (ini:fini)[i], "]")]
    if(n != nrow(md)){
      lams <- sample(lams, n, TRUE)
    } 
    yyy[,i] <- rpois(n, lams)
  }
  yyy
}


fig4 <- function(abunds, moon_dates){
  tiff("fig4.tiff", width = 6, height = 3, units = "in", res = 200)

  par(fig = c(0, 1, 0, 1), mar = c(1.25, 2.125, 1, 1))
  daterange <- as.Date(c("1993-08-01", "2019-01-01"))
  blank(bty = "L", xlim = daterange, ylim = c(-0.5, 18))

  rectx1 <- as.Date(moon_dates[300]) - 14
  rectx2 <- as.Date(moon_dates[500]) + 14
  rectcol1 <- rgb(0.7, 0.7, 0.7, 0.4)
  rect(rectx1, -0.3, rectx2, 18, col = rectcol1, border = NA)
  rectx1 <- as.Date(moon_dates[501]) - 14
  rectx2 <- as.Date(moon_dates[512]) + 14
  rectcol1 <- rgb(0.2, 0.2, 0.2, 0.4)
  rect(rectx1, -0.3, rectx2, 18, col = rectcol1, border = NA)


  x <- as.Date(moon_dates[200:512])
  y <- abunds[200:512]
  set.seed(123)
  yy <- y + runif(length(y), -0.25, 0.25)
  nas <- which(is.na(y))
  points(x[-nas], yy[-nas], type = "l", col = rgb(0.1, 0.1, 0.1, 0.5))
  points(x, yy, cex = 0.6, pch = 16, col = rgb(1, 1, 1, 1))
  points(x, yy, cex = 0.6, col = rgb(0.1, 0.1, 0.1, 0.5))
  axis(2, at = seq(0, 15, 5), labels = FALSE, tck = -0.0225)
  axis(2, at = seq(0, 15, 5), las = 1, lwd = 0, line = -0.5, cex.axis = 0.7)
  axis(2, at = seq(0, 18, 1), labels = FALSE, tck = -0.01)
  mtext(side = 2, line = 1.1, cex = 0.8, 
        expression(paste(italic(C.), " ", italic(penicillatus ), " counts")))
  datevals <- as.Date(paste0(1993:2019, "-01-01"))
  axis(1, at = datevals, labels = FALSE, tck = -0.01)
  datevals <- as.Date(paste0(seq(1995, 2019, 5), "-01-01"))
  datelabs <- seq(1995, 2019, 5)
  axis(1, at = datevals, labels = datelabs, lwd = 0, line = -0.9, 
       cex.axis = 0.7)
  axis(1, at = datevals, labels = FALSE, tck = -0.02)
  
  par(fig = c(0.1, 0.35, 0.5, 1), mar = c(1.25, 1.5, 1, 0.5), new = TRUE)
  blank(bty = "L", xlim = c(-0.5, 17.5), ylim = c(0, 0.4))
  tt <- table(abunds[200:512])
  points(tt/sum(tt))

  axis(1, at = seq(0, 15, 5), labels = FALSE, tck = -0.025)
  axis(1, at = seq(0, 15, 5), las = 1, lwd = 0, line = -1.25, cex.axis = 0.5)
  axis(1, at = seq(0, 18, 1), labels = FALSE, tck = -0.015)
  mtext(side = 1, line = 0.35, cex = 0.5, 
        expression(paste(italic(C.), " ", italic(penicillatus ), " counts")))

  axis(2, at = seq(0, 0.4, 0.1), labels = FALSE, tck = -0.025)
  axis(2, at = seq(0, 0.4, 0.1), las = 1, lwd = 0, line = -0.75,
       cex.axis = 0.5)
  axis(2, at = seq(0, 0.4, 0.05), labels = FALSE, tck = -0.01)
  mtext(side = 2, line = 1.1, cex = 0.5, "Frequency")

  dev.off()
}


load_models <- function(modns, small = FALSE){
  if (small){
    fnames <- paste0("model_output/mod", modns, "_withoutmodels.RData")
  } else{
    fnames <- paste0("model_output/mod", modns, "_withmodels.RData")
  }
  for(i in 1:length(modns)){
    load(fnames[i], envir = parent.frame(n = 2))
  }
}


summarize_abunds <- function(x, digits = 2){
  x2 <- na.omit(x)
  out <- c(n_all = length(x), n = length(x2), min = min(x2), max = max(x2), 
           median = median(x2),
           mean = mean(x2), var = var(x2), skewness = skewness(x2))
  round(out, digits)
}

#
# figure functions to be pulled out to bbplot are in other scripts




fig4_row <- function(x, draws, draws2, spits, frow, minx = 1, maxx = 30, 
                     nxs = 30, 
                     n_bins = 10, ym1 = 35, ym2 = 35, wex = 1.45, 
                     steps = 1:length(x), minstep = 490,
                     maxstep = 513, buff = 0.35, minp = 1e-3, vwidth = 0.5, 
                     xrange = 0:25, vxs = 501:512, vc = rgb(0.6, 0.6, 0.6), 
                     poc = rgb(0.4, 0.4, 0.4, 0.1),
                     nfrows = 5, jitv = 0.4, seed = 1234, title = NULL){
  set.seed(seed)
  topoffset <- 0.94
  rtop <- topoffset - ((frow - 1) * topoffset) / nfrows
  rbot <- max(c(0, rtop - topoffset / nfrows))

  par(fig = c(0, 0.5, rbot, rtop), mar = c(1, 1, 1.25, 0.5), 
       new = TRUE)

  blank(bty = "L", xlim = c(minstep - buff, maxstep + buff), ylim = c(0, ym1))
  nvxs <- length(vxs)
  for(j in 1:nvxs){
    vx <- vxs[j] - 500 
    violin(draws[, vx], vxs[vx], wex = wex, col = vc, border = NA)
  }
  points(steps[steps >= minstep], x[steps >= minstep], cex = 0.4, pch = 16)
  axis(side = 2, at = seq(0, ym1, 5), tck = -0.01, labels = FALSE)
  axis(side = 2, at = seq(0, ym1, 10), tck = -0.03, labels = FALSE)
  mtext(side = 2, at = seq(0, ym1, 10), seq(0, ym1, 10), xpd = NA, las = 1,
        line = 0.25, cex = 0.5)
  axis(side = 1, at = seq(490, 515, 5), tck = -0.01, labels = FALSE)
  axis(side = 1, at = seq(490, 510, 10), tck = -0.03, labels = FALSE)
  axis(side = 1, at = 489.1, tck = -0.08, labels = FALSE)
  mtext(side = 1, at = 489.1, line = 0., cex = 0.5, "2017")
  axis(side = 1, at = 501.4, tck = -0.08, labels = FALSE)
  mtext(side = 1, at = 501.4, line = 0., cex = 0.5, "2018")
  axis(side = 1, at = 513.8, tck = -0.08, labels = FALSE)
  mtext(side = 1, at = 513.8, line = 0., cex = 0.5, "2019")


  par(fig = c(0, 1, rtop - 0.1, rtop), new = TRUE)
  blank()
  text(x = 0.575, y = 1.65, title, cex = 0.75, xpd = NA, adj = 0)

  par(fig = c(0.5, 0.75, rbot, rtop), new = TRUE)

  blank(bty = "L", xlim = c(-1, ym2), ylim = c(-1, ym2))
  abline(a = 0 , b = 1)
  
  for(i in 1:nvxs){
    specs <- sample(1:nrow(draws2), 100)
    xx <- draws2[specs,i]
    yy <- rep(x[i+488], length(xx))
    xx2 <- xx + runif(length(yy), -jitv, jitv)
    yy2 <- yy + runif(length(yy), -jitv, jitv)
    points(xx2, yy2, col = poc, cex = 0.75)
  }

  axis(side = 1, at = seq(0, ym1, 5), tck = -0.01, labels = FALSE)
  axis(side = 1, at = seq(0, ym1, 10), tck = -0.03, labels = FALSE)
  mtext(side = 1, at = seq(0, ym1, 10), seq(0, ym1, 10), xpd = NA, las = 1,
        line = -0.25, cex = 0.5)
  axis(side = 2, at = seq(0, ym1, 5), tck = -0.01, labels = FALSE)
  axis(side = 2, at = seq(0, ym1, 10), tck = -0.03, labels = FALSE)
  mtext(side = 2, at = seq(0, ym1, 10), seq(0, ym1, 10), xpd = NA, las = 1,
        line = 0.25, cex = 0.5)

  par(fig = c(0.75, 1, rbot, rtop), new = TRUE)

  blank(xlim = c(1, 10), ylim = c(0, max(spits) * 1.25), bty = "L")
  abline(h = 1, lty = 3)
  points(spits, type = "h", lwd = 2)

}



fig3_row <- function(x, draws, frow, minx = 1, maxx = 30, nxs = 30, 
                     n_bins = 10, ym1 = 35, ym2 = 45, wex = 1.45, 
                     steps = 1:length(x), minstep = 1,
                     maxstep = 30, buff = -0.5, minp = 1e-3, vwidth = 0.5, 
                     xrange = 0:25, vxs = 1:30, vc = rgb(0.6, 0.6, 0.6), 
                     poc = rgb(0, 0, 0, 0.025),
                     nfrows = 7, jitv = 0.4, seed = 1234, title = NULL){
  set.seed(seed)
  topoffset <- 0.94
  rtop <- topoffset - ((frow - 1) * topoffset) / nfrows
  rbot <- max(c(0, rtop - topoffset / nfrows))

  par(fig = c(0, 0.6, rbot, rtop), mar = c(0.5, 0.5, 1.25, 0.5), 
       new = TRUE)

  blank(bty = "L", xlim = c(minstep - buff, maxstep + buff), ylim = c(0, ym1))
  nvxs <- length(vxs)
  for(j in 1:nvxs){
    vx <- vxs[j]
    violin(draws[, vx], vx, wex = wex, col = vc, border = NA)
  }
  points(steps[steps >= minstep], x[steps >= minstep], cex = 0.4, pch = 16)

  par(fig = c(0, 1, rtop - 0.1, rtop), new = TRUE)
  blank()
  text(x = 0.575, y = 1.65, title, cex = 0.75, xpd = NA, adj = 0)

  par(fig = c(0.6, 0.8, rbot, rtop), new = TRUE)

  blank(bty = "L", xlim = c(0, ym2), ylim = c(0, ym2))
  abline(a = 0 , b = 1)
  
  for(i in 1:nlams){
    specs <- sample(1:nrow(draws), 100)
    xx <- draws[specs,i]
    yy <- rep(x[i], length(xx))
    xx2 <- xx + runif(length(yy), -jitv, jitv)
    yy2 <- yy + runif(length(yy), -jitv, jitv)
    points(xx2, yy2, col = poc, cex = 0.75)
  }
  par(fig = c(0.8, 1, rbot, rtop), new = TRUE)

  cdf <- function(x, lam){ecdf(lam)(x)}

  nx <- length(x)
  pits <- matrix(NA, nrow = nx, ncol = n_bins)
  for(i in 1:nlams){
    pits[i,] <- nrPIT(x[i], cdf, n_bins, draws[,i])
  }
  
  spits <- apply(pits, 2, mean)
  blank(xlim = c(1, 10), ylim = c(0, max(spits) * 1.25), bty = "L")
  abline(h = 1, lty = 3)
  points(spits, type = "h", lwd = 2)

}



mass <- function(x, min = NULL, max = NULL){
  if (!all(x %% 1 == 0)){
    stop("x must be integer conformable")
  }
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
make_eval_tab <- function(mod1, mod2, mod3){
  npitc <- ncol(mod1[[1]]$eval$PIT)
  evalmat <- matrix(NA, nrow = 3 * length(mod1) * 12, ncol = 6 + npitc)
  spots <- 0

  for(i in 1:length(mod1)){
    for(j in 1:3){
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
                  mu0 = rnorm(1, 0, 1),
                  tau = rgamma(1, shape = 0.1, rate = 0.1))
           }
  }
  if (model == 2){
    out <- function(chain = chain){
             list(.RNG.name = sample(rngs, 1),
                  .RNG.seed = sample(1:1e+06, 1),
                  mu0 = rnorm(1, 0, 1),
                  tau = rgamma(1, shape = 0.1, rate = 0.1),
                  phi = runif(1, -0.95, 0.95))
           }
  }
  if (model == 3){
    out <- function(chain = chain){
             list(.RNG.name = sample(rngs, 1),
                  .RNG.seed = sample(1:1e+06, 1),
                  mu0 = rnorm(1, 0, 1),
                  tau = rgamma(1, shape = 0.1, rate = 0.1),
                  phi = runif(1, -0.95, 0.95),
                  beta1 = rnorm(1, 0, 2),
                  beta2 = rnorm(1, 0, 2))
    }
  }
  out
}

monitor_fun <- function(model, meta = NULL, obs = TRUE, obstype = "lead"){
  if (model == 1){
    out <- c("sd", "mu0")
  }
  if (model == 2){
    out <- c("sd", "mu0", "phi")
  }
  if (model == 3){
    out <- c("sd", "mu0", "phi", "beta1", "beta2")
  }
  if (obs){
    if (obstype == "lead"){
      ptimes <- max(meta$in_timeseries) + 1:meta$lead_time - 
                min(meta$in_timeseries) + 1
      preds <- paste0("predY[", ptimes, "]")
      out <- c(out, preds)
    }
    if (obstype == "all"){
      out <- c(out, "predY")
    }
  }
  out
}

foy <- function(dates){
  dates <- as.Date(dates)
  yr <- format(dates, "%Y")
  nye <- as.Date(paste0(yr, "-12-31"))
  jnye <- as.numeric(format(nye, "%j"))
  jdates <- as.numeric(format(dates, "%j"))
  round(jdates/jnye, 4)
}

data_fun <- function(model, abunds, in_timeseries, lead_time, 
                     moon_dates = NULL){
  if(model %in% 1:2){
    Y <- c(abunds[in_timeseries], rep(NA, lead_time))
    meanY <- mean(Y, na.rm = TRUE)
    out <- list(Y = Y, N = length(Y), meanY = meanY)
  }
  if (model == 3){
    all_in <- (min(in_timeseries)):(max(in_timeseries) + lead_time)    
    fryr <- foy(moon_dates[all_in])
    Y <- c(abunds[in_timeseries], rep(NA, lead_time))
    meanY <- mean(Y, na.rm = TRUE)
    out <- list(Y = Y, N = length(Y), meanY = meanY,
                sinyr = sin(fryr), cosyr = cos(fryr))
  }
  out 
}

name_fun <- function(model){
  paste0("model_scripts/mod", model, ".txt")
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

    ptimes <- max(in_timeseries) + 1:lead_time - min(in_timeseries) + 1
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
        s_r[i] <- -crps_sample(y = yyy, dat = dist)
        s_l[i] <- log(length(dist[dist==yyy])/length(dist))
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
                 save_preds = FALSE, save_raw = TRUE, save_small = TRUE){
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
  modname <- paste0("mod", model)
  assign(modname, out)
  if (save_raw){
    fname <- paste0("model_output/mod", model, "_withmodels.RData")
    save(list = modname, file = fname)
  }
  if (save_small){
    fname <- paste0("model_output/mod", model, "_withoutmodels.RData")


    out2 <- vector("list", length = length(out))
    for(i in 1:length(out)){
      out2[[i]] <- list(summary = out[[i]]$summary, data = out[[i]]$data,
                        meta = out[[i]]$meta, eval = out[[i]]$eval)
    }
    names(out2) <- names(out)
    assign(modname, out2)
    save(list = modname, file = fname)
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


# graphics functions for drawing a variety of violins
# still in development but pretty and flexible right now
# also will get pulled over into bbplot 

#' @param x The vector of values or function to be summarized
#' @param location The graphical location with respect to the axes where
#'   the violin is to be placed. Default assumption is for the x-axis, but
#'   can be a named-by-axis vector.
#' @param type Character indicating the type of plotting. Defaults to NULL
#'   which defines itself then based on the support for x to 
#'   either "l" (line) for continuous and "r" (histogram-like rectangles) for 
#'   non-continuous/integer-conformable
#' @param wex Numeric height scale that transforms the distribution
#'   evalution to the plotting axis value. now wex for width expansion
#' @param nvalues Integer number of values to use for the drawing of the 
#'   violin. If NULL, set by vtype. 
# also now have values explicitly
# rotate T/F, wrt the plotting axes
# side top bottom both


violin <- function(x, location = NULL, rotate = TRUE,
                   type = NULL, wex = 1, values = NULL, nvalues = NULL, 
                   side = "both", ...){

  vlocation <- violin_location(location, rotate)
  vtype <- if_null(type, default_violin_type(x))
  dvals <- dist_values(x, values, nvalues, vtype)
  vvals <- violin_values(dvals, vlocation, rotate, vtype, wex, side)
  draw_violin(vvals, vtype, ...)
}

draw_violin <- function(values, type = NULL, ...){
  if (type == "l"){
    polygon(values, ...)
  }
  if (type == "r"){
    nvalues <- nrow(values) / 2
    for(i in 1:nvalues){
      row1 <- 2 * i - 1
      row2 <- 2 * i 
      xleft <- values[row1, "x"]
      ybottom <- values[row1, "y"]
      xright <- values[row2, "x"]
      ytop <- values[row2, "y"]
      rect(xleft, ybottom, xright, ytop, ...)
    }
  }
}

violin_values <- function(dist_values, location = NULL, rotate = TRUE, 
                          type = "l", wex = 1, side = "both"){

  length_vals <- violin_length_values(dist_values, type, side)
  width_vals <- violin_width_values(dist_values, type, wex, side)
  if (rotate){
    y_vals <- length_vals + location["y"]
    x_vals <- width_vals + location["x"]
  } else{
    x_vals <- length_vals + location["x"]
    y_vals <- width_vals + location["y"]
  }
  data.frame(x = x_vals, y = y_vals)
}

violin_length_values <- function(dist_values, type = "l", side = "both"){
  nvalues <- nrow(dist_values)
  values <- dist_values[ , "x"]
  if (type == "l"){
    if (side == "both"){
      c(values, values[nvalues:1])
    } else {
      values
    }
  } else if (type == "r"){
    rwex <- 0.45
    value_diff <- min(diff(values))
    x_offset <- value_diff * c(-1, 1) * rwex
    rep(values, each = 2) + rep(x_offset, nvalues)
  }
}

violin_width_values <- function(dist_values, type = "l", wex = 1, 
                                side = "both"){
  nvalues <- nrow(dist_values)
  values <- dist_values[ , "y"] * wex
  if (type == "l"){
    if (side == "both"){
      c(values, -values[nvalues:1])
    } else {
      values_sign <- switch(side, "pos" = 1, "neg" = -1)
      values * values_sign
    }
  } else if (type == "r"){
    bottom_mult <- switch(side, "both" = -1, "pos" = 0, "neg" = -1)
    top_mult <- switch(side, "both" = 1, "pos" = 1, "neg" = 0)
    rep(values, each = 2) * c(bottom_mult, top_mult) * wex
  }
}

#' determines the values of the variable to evaluate (x_values) and the 
#' resulting evaluation values (y_values) of the distribution. these are 
#' the raw values
#'
#' @return data.frame of x and y  
#'
#' @export
#' 
dist_values <- function(x, values = NULL, nvalues = NULL, type = "l"){

  if (is.null(values)) {
    if (is.null(nvalues)) {
      nvalues <- default_nvalues(x, type)
    }
    xvals <- dist_x_values(x, nvalues)
  } else {
    xvals <- values
  }
  yvals <- dist_y_values(x, xvals)
  data.frame(x = xvals, y = yvals)
}


dist_x_values <- function(x, nvalues){
  minx <- min(x, na.rm = TRUE)
  maxx <- max(x, na.rm = TRUE)
  seq(minx, maxx, length.out = nvalues)
}

dist_y_values <- function(x, xvals){
  if (all(x %% 1 == 0) & all(xvals %% 1 == 0)){
    yvals <- mass(x)$y
  } else{ 
    den <- density(x)
    nvalues <- length(xvals)
    yvals <- rep(NA, nvalues)
    for(i in 1:nvalues){
      match_less <- which(den$x < xvals[i])
      match_more <- which(den$x > xvals[i])
      match_hit <- which(den$x == xvals[i])
      nmatch_less <- length(match_less)
      nmatch_more <- length(match_more)
      nmatch_hit <- length(match_hit)
      if (nmatch_hit == 1){
        yvals[i] <- mean(den$y[match_hit])
      } else if (nmatch_less == 0 | nmatch_more == 0){
        yvals[i] <- 0
      } else{
        xval_val_1 <- den$x[match_less[nmatch_less]]
        xval_val_2 <- den$x[match_more[1]]
        yval_val_1 <- den$y[match_less[nmatch_less]]
        yval_val_2 <- den$y[match_more[1]]
        xval_diff_1 <- xvals[i] - xval_val_1
        xval_diff_2 <- xval_val_2 - xvals[i] 
        xval_diff_12 <- xval_val_2 - xval_val_1
        xval_ratio_1 <- 1 - xval_diff_1 / xval_diff_12
        xval_ratio_2 <- 1 - xval_diff_2 / xval_diff_12
        yvals[i] <- (yval_val_1 * xval_ratio_1 + yval_val_2 * xval_ratio_2) /2
      }
    }
  }  
  yvals
}

if_null <- function(x = NULL, val_if_null = NULL){
  if (is.null(x)){
    val_if_null
  } else {
    x
  }
}



#' nrf n reduction factor

default_nvalues <- function(x, type = "l", nrf = NULL, minn = NULL, 
                            maxn = NULL){
  if (type == "n"){
    10
  } else {
    if (all(x %% 1 == 0)){
      length(seq(min(x), max(x), 1))
    } else{
      if (type == "l"){
        nrf <- if_null(nrf, 100)
        minn <- if_null(minn, 100)
        maxn <- if_null(maxn, 1000)
      }
      if (type == "r"){
        nrf <- if_null(nrf, 100)
        minn <- if_null(minn, 2)
        maxn <- if_null(maxn, 10)
      }

      min(c(max(c((length(x) / nrf), minn)), maxn))
    }
  }
}





default_violin_type <- function(x = NULL){
  type <- "l"
  if (!is.numeric(x)){
    stop("presently only supported for numeric values")
  } 
  if (all(x %% 1 == 0)){
    type <- "r"
  }
  type
}


violin_location <- function(location = NULL, rotate = TRUE){
  if(is.null(location)){
    out <- c(x = 0, y = 0)
  } else {
    if (is.null(names(location))){
      out <- c(x = location[1], y = location[2])
      out[which(is.na(out))] <- 0      
    } else {
      out <- c(location["x"], location["y"])
      out[which(is.na(out))] <- 0
    }
  }
  if(!rotate){
    names(out) <- c("y", "x")
  } else{
    names(out) <- c("x", "y")
  }
  out
}



blank <- function(x = 1, y = 1, type = "n", xlab = "", ylab = "", 
                  xaxt = "n", yaxt = "n", bty = "n", ...){
  plot(x = x, y = y, type = type, xlab = xlab, ylab = ylab, 
                  xaxt = xaxt, yaxt = yaxt, bty = bty, ...)
}


