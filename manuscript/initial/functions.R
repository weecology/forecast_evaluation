load_dependencies <- function(){
  depend <- c("portalr", "coda", "runjags", "scoringRules")
  ndepend <- length(depend)
  present <- installed.packages()[ , "Package"]
  needed <- depend[!(depend %in% present)] 
  nneeded <- length(needed)
  if(nneeded > 0){
    install.packages(needed)
  }
  for(dep in 1:ndepend){
    eval(bquote(library(.(depend[dep]))))
  }
}


prep_data_set <- function(remote = TRUE, update = TRUE){
  if(remote){
    abunds <- plot_sp_abunds(plot = 19, species = "PP")
    moon_dates <- load_trapping_data(clean = FALSE)$newmoons_table$newmoondate

    if(update){
      datatable <- data.frame(abunds, moon_dates)
      write.table(datatable, "data_set.txt", sep = ",", row.names = FALSE)
    }
  }
  else{
    datatable <- read.table("data_set.txt", sep = ",", header = TRUE)
    abunds <- as.numeric(datatable$abunds)
    moon_dates <- as.character(datatable$moon_dates)
  }
  list(abunds = abunds, moon_dates = moon_dates)
}

run_models <- function(data_set = NULL, models = 1:3, starts = 200, 
                       ends = 300:500, run = TRUE, small = FALSE){
  if(run == FALSE){
    out <- load_models(models, small)
    return(out)  
  }
  dir.create("model_output", showWarnings = FALSE)
  in_ts <- training_ts(starts, ends)
  nmodels <- length(models)
  out <- vector("list", length = nmodels)
  for(i in 1:nmodels){
    out[[i]] <- mods(i, data_set, in_ts)
  }
  names(out) <- paste0("mod", models)  
  out
}

load_models <- function(models, small = FALSE){
  if (small){
    fnames <- paste0("model_output/mod", models, "_withoutmodels.rds")
  } else{
    fnames <- paste0("model_output/mod", models, "_withmodels.rds")
  }
  nmodels <- length(models)
  out <- vector("list", length = nmodels)
  for(i in 1:nmodels){
    out[[i]] <- readRDS(fnames[i])
  }
  names(out) <- paste0("mod", models)  
  out
}



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
  data.frame(evalmat)
}

mods <- function(model = NULL, data_set = NULL, 
                 in_timeseries = list(1:500), lead_time = 12, 
                 n_chains = 3, n_burnin = 5000, n_sample = 10000, n_thin = 1,
                 n_bins = 10, n_adapt = 1000, keeploc = FALSE, quiet = FALSE,
                 save_preds = FALSE, save_raw = TRUE, save_small = TRUE){
  abunds <- data_set$abunds
  moon_dates <- data_set$moon_dates
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
    fname <- paste0("model_output/mod", model, "_withmodels.rds")
    save(modname, file = fname)
  }
  if (save_small){
    fname <- paste0("model_output/mod", model, "_withoutmodels.rds")


    out2 <- vector("list", length = length(out))
    for(i in 1:length(out)){
      out2[[i]] <- list(summary = out[[i]]$summary, data = out[[i]]$data,
                        meta = out[[i]]$meta, eval = out[[i]]$eval)
    }
    names(out2) <- names(out)
    assign(modname, out2)
    saveRDS(modname, file = fname)
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
      ys[which(ys[,i]>49),i] <- 49
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

name_fun <- function(model){
  paste0("mod", model, ".txt")
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



# function for nonrandomized PIT histogram 
# inputs: 
#   x    observed data 
#   cdf   CDF
#   n_bins number of bins to use
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



if_null <- function(x = NULL, val_if_null = NULL){
  if (is.null(x)){
    val_if_null
  } else {
    x
  }
}

messageq <- function(txt = NULL, quiet = FALSE){
  if (!quiet){
    message(txt)
  }
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

blank <- function(x = 1, y = 1, type = "n", xlab = "", ylab = "", 
                  xaxt = "n", yaxt = "n", bty = "n", ...){
  plot(x = x, y = y, type = type, xlab = xlab, ylab = ylab, 
                  xaxt = xaxt, yaxt = yaxt, bty = bty, ...)
}

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



fig1top <- function(){
  tiff("fig1top.tiff", width = 6, height = 4, units = "in", res = 200)
  rc <- grey(0.6)

  par(mar = c(2.5, 1.75, 0, 0.75), fig = c(0, 1, 0.007, 0.49))


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
  points(x4, y3, pch = 1, cex = 1.2, lwd = 2, col = grey(0,0.9))
  set.seed(123)
  fx <- c(0.995, 1.08, 1.15, 1.225)
  fn <- length(fx)
  fy1 <- 0.35 + 0.175 * sin(fx/0.75 * 2 * pi)
  fy2 <- fy1 + rnorm(fn, 0, 0.1)

  mtext(side = 2, "y ~ G", line = 0.75)
  mtext(side = 1, "Time", line = 0.6)
  axis(side = 1, at = seq(0, 1.25, 0.05), tck = -0.02, labels = FALSE, 
       lwd = 2)
  axis(side = 1, at = seq(0, 1.25, 0.25), tck = -0.03, labels = FALSE, 
       lwd = 2)
  axis(side = 2, at = seq(0, 1, 0.25), tck = -0.02, labels = FALSE, lwd = 2)
  axis(side = 2, at = seq(0, 1, 0.5), tck = -0.03, labels = FALSE, lwd = 2)

  arrows(x4[n-1] - 0.05, y3[n-1] + 0.2, x4[n-1] - 0.005, y3[n-1] + 0.05, 
         length = 0.05, lwd = 2)
  text(x4[n-1] - 0.0625, y3[n-1] + 0.275, cex = 0.8, 
        expression(italic("y"[o-1])))
  mtext(side = 1, at = x4[n-1], cex = 0.8, line = 0.5, 
        expression(italic("t"[o-1])))
  points(rep(x4[n-1], 2), c(-0.1, y3[n-1] - 0.02), type = "l", lty = 2, 
         lwd = 2, xpd = TRUE, cex = 1.2, col = grey(0,0.9))
  points(rep(x4[n], 2), c(-0.1, y3[n] - 0.02), type = "l", lty = 2, lwd = 2, 
         xpd = TRUE, cex = 1.2, col = grey(0,0.9))
  mtext(side = 1, at = x4[n], cex = 0.8, line = 0.5, 
        expression(italic("t"[o])))

  points(fx, fy2, pch = 1, cex = 1.2, lwd = 2, col = grey(0,0.9))
  mtext(side = 1, at = fx[fn] + 0.03, cex = 0.8, line = 0.5, 
        expression(italic("t"[o+P]*" = "*"t"[N])))
  points(rep(fx[fn], 2), c(-0.1, fy2[fn] - 0.02), type = "l", lty = 2, 
         lwd = 2, xpd = TRUE)

  points(c(x4[n], fx[fn]), rep(-0.31, 2), type = "l", lwd = 2, xpd = TRUE)
  points(rep(fx[fn], 2), c(-0.36, -0.26), type = "l", lwd = 2, xpd = TRUE)
  points(rep(x4[n], 2), c(-0.36, -0.26), type = "l", lwd = 2, xpd = TRUE)
  text(x = mean(c(x4[n], fx[fn])), y = -0.23, "forecast horizon", 
       font = 3, cex = 0.6, xpd = TRUE)

  text(x = x4[n] - 0.1375, y = -0.38, "forecast origin", 
       font = 3, cex = 0.6, xpd = TRUE)

  arrows(x4[n] - 0.065, -0.32, x4[n] - 0.02, -0.20, length = 0.05, xpd = TRUE,
         lwd = 2)

  arrows(x4[1] + 0.02, y3[1] + 0.5, x4[1] + 0.002, y3[1] + 0.05, 
         length = 0.05, lwd = 2)
  text(x4[1] + 0.025, y3[1] + 0.55, cex = 0.8, expression(italic("y"[1])))
  mtext(side = 1, at = x4[1], cex = 0.8, line = 0.5, 
        expression(italic("t"[1])))
  points(rep(x4[1], 2), c(-0.1, y3[1] - 0.02), type = "l", lty = 2, lwd = 2, 
         xpd = TRUE)

  u <- par("usr") 
  arrows(u[1], u[3], u[2], u[3], length = 0.1, lwd = 2, xpd = TRUE) 
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
  points(fx, fy2, pch = 16, cex = 1.2, lwd = 2, col = grey(1,0.9))
  points(fx, fy2, pch = 1, cex = 1.2, lwd = 2, col = grey(0,0.9))

  text(fx[1] + 0.04, fy2[1] + 0.45, cex = 0.8, 
        expression(italic("y"[o+1])))
  arrows(fx[1] + 0.02, fy2[1] + 0.4, fx[1] + 0.002, fy2[1] + 0.05, 
         length = 0.05, lwd = 2)

  par(mar = c(0, 0, 0, 0), fig = c(0, 1, 0, 0.72), new = TRUE)
  blank(ylim = c(0, 1), xlim = c(0, 1))
  text(x = 0, y = 1, "(b)", cex = 1.1, xpd = TRUE)


  par(mar = c(0.5, 1.25, 0.5, 0), fig = c(0.475, 0.675, 0.483, 0.693), 
      new = TRUE)

  set.seed(123)
  Gn <- c(rnorm(1e4, y3[n-1] - 0.05, 0.1), rnorm(7e3, y3[n-1] - 0.1, 0.25))
  Gn2 <- Gn[-which(Gn < 0.05)]
  dens <- density(Gn2)
  blank(bty = "L", xlim = c(0,1), ylim = c(0,max(dens$y)*1.01))
  box(bty = "L", lwd = 2)
  points(dens$x, dens$y, type = "l", lwd = 2)
  axis(side = 1, at = seq(0, 1.25, 0.25), tck = -0.025, labels = FALSE, 
       lwd = 2)
  axis(side = 1, at = seq(0, 1, 0.5), tck = -0.05, labels = FALSE, lwd = 2)
  mtext(side = 1, at = 0.4, cex = 0.8, line = 0.25, expression("G"[o-1]))

  u <- par("usr") 
  yatx <- min(dens$y[which(round(dens$x, 2) == round(y3[n-1], 2))])
  points(rep(y3[n-1], 2), c(u[3], yatx), type = "l", lwd = 2, lty = 3)
  points(c(y3[n-1], y3[n-1] + 0.13), c(u[3], -1.6), type = "l", lwd = 2, 
         lty = 3, xpd = NA)
  mtext(side = 2, cex = 0.7, line = 0.25, "Density") 

  par(mar = c(0.5, 1.25, 0.5, 0), fig = c(0.05, 0.25, 0.483, 0.693), 
      new = TRUE)

  set.seed(123)
  Gn <- c(rnorm(1e4, y3[1] + 0.05, 0.1), rnorm(7e3, y3[1] + 0.1, 0.25))
  Gn2 <- Gn[-which(Gn < 0.05)]
  dens <- density(Gn2)
  blank(bty = "L", xlim = c(0,1), ylim = c(0,max(dens$y)*1.01))
  box(bty = "L", lwd = 2)
  points(dens$x, dens$y, type = "l", lwd = 2)
  axis(side = 1, at = seq(0, 1.25, 0.25), tck = -0.025, labels = FALSE, 
       lwd = 2)
  axis(side = 1, at = seq(0, 1, 0.5), tck = -0.05, labels = FALSE, lwd = 2)
  mtext(side = 1, cex = 0.8, line = 0.25, expression("G"[1]))
 
  u <- par("usr") 
  yatx <- min(dens$y[which(round(dens$x, 2) == round(y3[1], 2))])

  points(rep(y3[1], 2), c(u[3], yatx), type = "l", lwd = 2, lty = 3)
  points(c(y3[1], y2[1] - 0.2), c(u[3], -1.6), type = "l", lwd = 2, lty = 3, 
         xpd = NA)
  mtext(side = 2, cex = 0.7, line = 0.25, "Density") 

  par(mar = c(0.5, 1.25, 0.5, 0), fig = c(0.75, 0.9755, 0.483, 0.693), 
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
  axis(side = 1, at = seq(0, 1.25, 0.25), tck = -0.025, labels = FALSE, 
       lwd = 2)
  axis(side = 1, at = seq(0, 1, 0.5), tck = -0.05, labels = FALSE, lwd = 2)
  mtext(side = 1, cex = 0.8, at = 0.75, line = 0.25, 
        expression("H"[o+1]*", "*"G"[o+1]))
  yatx <- min(dens$y[which(round(dens$x, 2) == round(fy2[1], 2))])

  u <- par("usr") 
  points(rep(fy2[1], 2), c(u[3], yatx), type = "l", lwd = 2, lty = 3)
  points(c(fy2[1], fy2[1] - 0.475), c(u[3], -1.1), type = "l", lwd = 2, 
         lty = 3, xpd = NA)
  mtext(side = 2, cex = 0.7, line = 0.25, "Density") 



  par(mar = c(0, 0, 0, 0), fig = c(0, 1, 0.67, 1), new = TRUE)
  blank(ylim = c(0, 1), xlim = c(0, 1))
  text(x = 0, y = 0.95, "(a)", cex = 1.1, xpd = TRUE)
  set.seed(123)
  x1 <- rnorm(1e6)
  x2 <- rnorm(1e6, 0, 1.5)
  x3 <- rt(1e6, 3.5, 0)
  x4 <- c(rnorm(3e5, 2.25, 2), rnorm(7e5, -0.965, 1))
  x5 <- c(rnorm(5e5, -3), rnorm(5e5, 3))
  x6 <- c(rnorm(3e5, -5, 0.5), rnorm(3e5, 5, 0.5), rnorm(4e5, 0))

  par(mar = c(1.5, 2.5, 1.0, 1.5), c(0, 1, 0.7, 1), new = TRUE)
  blank(bty = "L", xlim = c(-8, 8), ylim = c(0, 0.55))
  box(bty = "L", lwd = 2)

  points(density(x1), type = "l", lwd = 2, col = grey(0, 0.8))
  points(density(x2), type = "l", lwd = 2, col = grey(0.1, 0.8))
  points(density(x3), type = "l", lwd = 2, col = grey(0.2, 0.8))
  points(density(x4), type = "l", lwd = 2, col = grey(0.3, 0.8))
  points(density(x5), type = "l", lwd = 2, col = grey(0.4, 0.8))
  points(density(x6), type = "l", lwd = 2, col = grey(0.5, 0.8))
  axis(1, at = seq(-10, 10, 1), tck = -0.04, labels = FALSE)
  axis(1, at = seq(-10, 10, 5), tck = -0.08, labels = FALSE)

  mtext(side = 2, line = 0.5, "Density", cex = 0.75)


  dev.off()
}

fig1bottom <- function(){

  tiff("fig1bottom.tiff", width = 6, height = 3, units = "in", res = 200)
  rc <- rgb(0.8, 0.8, 0.8)
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

  text(-0.8, 2.55, "(c)", cex = 1.1, font = 1, xpd = NA)

  text(-1.1, 2.06, "Single Origin", cex = 0.7, font = 3, xpd = NA, adj = 0)

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

  text(-1.1, 2.06, "Rolling Origin", cex = 0.7, font = 3, xpd = NA, adj = 0)

  dev.off()
}

figA1 <- function(){
  set.seed(321)

  xs <- seq(1, 50, length.out = 30)
  lams <- 8 + 0.25 * xs + 3 * sin(2 * pi * xs / 15)
  nlams <- length(lams)
  x <- rpois(nlams, lambda = lams)
  N <- 1e4

  tiff("figA1.tiff", width = 6, height = 7, units = "in", res = 200)

  topoffset <- 0.94
  par(fig = c(0, 0.6, topoffset, 1), mar = c(0, 0, 0, 0))
  blank()
  text(x = 1, y = 1.15, "Time series of predictive distributions and", 
       cex = 0.8)
  text(x = 1, y = 0.825, "true observations", cex = 0.8, font = 1)
  par(fig = c(0.6, 0.8, topoffset, 1), new = TRUE)
  blank()
  text(x = 1, y = 1.15, "Observed (y) vs.", cex = 0.8, font = 1)
  text(x = 1, y = 0.825, "Predicted (x)", cex = 0.8, font = 1)
  par(fig = c(0.8, 1, topoffset, 1), new = TRUE)
  blank()
  text(x = 1, y = 1.15, "Probability Integral", cex = 0.8, font = 1)
  text(x = 1, y = 0.825, "Transform (PIT)", cex = 0.8, font = 1)

  draws1 <- sapply(lams, rpois, n = N)
  figA1_row(x, draws1, 1, title = "True generating distribution")

  draws2 <- sapply(lams + 2, rpois, n = N)
  figA1_row(x, draws2, 2, title = "Positively biased")

  draws3 <- sapply(lams - 2, rpois, n = N)
  figA1_row(x, draws3, 3, title = "Negatively biased")

  draws4 <- sapply(x, rpois, n = N)
  figA1_row(x, draws4, 4, title = "Too accurate")

  draws5 <- matrix(NA, nrow = N, ncol = nlams)
  for(i in 1:nlams){
    draws5[,i] <- round(rnorm(N, lams[i], sd = sqrt(lams[i]) / 1.6))
  }
  figA1_row(x, draws5, 5, title = "Too precise")

  draws6 <- matrix(NA, nrow = N, ncol = nlams)
  for(i in 1:nlams){
    draws6[,i] <- rnbinom(N, mu = lams[i], size = 1)
  }
  figA1_row(x, draws6, 6, title = "Too imprecise")

  draws7 <- matrix(NA, nrow = N, ncol = nlams)
  for(i in 1:nlams){
    draws7[,i] <-  c(rpois(N/2, max(c(0, lams[i] - 5))), 
                     rpois(N/2, lams[i] + 5))
  }
  figA1_row(x, draws7, 7, title = "Bimodal with proper median")

  dev.off()
}


figA1_row <- function(x, draws, frow, minx = 1, maxx = 30, nxs = 30, 
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
  nlams <- dim(draws)[2]
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



fig2top <- function(data_set, models){

  abunds <- data_set$abunds
  moon_dates <- data_set$moon_dates
  mod1 <- models$mod1
  mod2 <- models$mod2 
  mod3 <- models$mod3

  tiff("fig2top.tiff", width = 6, height = 4.5, units = "in", res = 200)

  par(fig = c(0, 1, 0, 0.7), mar = c(1.25, 2.125, 0.5, 1))
  daterange <- as.Date(c("1993-08-01", "2019-01-01"))
  blank(bty = "L", xlim = daterange, ylim = c(-0.5, 18))

  rectx1 <- as.Date(moon_dates[300]) - 14
  rectx2 <- as.Date(moon_dates[500]) + 14
  rectcol1 <- rgb(0.9, 0.9, 0.9, 1)
  rect(rectx1, -0.3, rectx2, 18, col = rectcol1, border = NA)
  rectx1 <- as.Date(moon_dates[501]) - 14
  rectx2 <- as.Date(moon_dates[512]) + 14
  rectcol1 <- rgb(0.8, 0.8, 0.8, 1)
  rect(rectx1, -0.3, rectx2, 18, col = rectcol1, border = NA)


  x <- as.Date(moon_dates[200:512])
  y <- abunds[200:512]
  set.seed(123)
  yy <- y
  nas <- which(is.na(y))
  points(x[-nas], yy[-nas], type = "l", col = rgb(0.1, 0.1, 0.1, 0.6))
  points(x, yy, cex = 0.45, pch = 16, col = rgb(1, 1, 1, 1))
  points(x, yy, cex = 0.45, col = rgb(0.1, 0.1, 0.1, 0.6))
  axis(2, at = seq(0, 10, 10), labels = FALSE, tck = -0.025)
  axis(2, at = seq(0, 10, 10), las = 1, lwd = 0, line = -0.5, cex.axis = 0.8)
  axis(2, at = seq(0, 18, 5), labels = FALSE, tck = -0.02)
  axis(2, at = seq(0, 18, 1), labels = FALSE, tck = -0.01)
  mtext(side = 2, line = 1.45, cex = 0.6, 
        expression(paste(italic(C.), " ", italic(penicillatus ), " counts")))
  datevals <- as.Date(paste0(1993:2019, "-01-01"))
  axis(1, at = datevals, labels = FALSE, tck = -0.02)
  datevals <- as.Date(paste0(seq(1995, 2019, 5), "-01-01"))
  axis(1, at = datevals, labels = FALSE, tck = -0.03)
  datevals <- as.Date(paste0(seq(1995, 2015, 5), "-01-01"))
  datelabs <- seq(1995, 2015, 5)
  axis(1, at = datevals, labels = datelabs, lwd = 0, line = -0.8, 
       cex.axis = 0.7)
  axis(1, at = datevals, labels = FALSE, tck = -0.02)
  
  par(fig = c(0.1, 0.4, 0.4, 0.65), mar = c(1.25, 1.5, 0, 0.5), new = TRUE)
  blank(bty = "L", xlim = c(-0.5, 17.5), ylim = c(0, 0.4))
  tt <- table(abunds[200:512])
  points(tt/sum(tt), lwd = 3)

  axis(1, at = seq(0, 15, 5), labels = FALSE, tck = -0.025)
  axis(1, at = seq(0, 15, 10), las = 1, lwd = 0, line = -1.25, cex.axis = 0.5)
  axis(1, at = seq(0, 18, 1), labels = FALSE, tck = -0.015)
  mtext(side = 1, line = 0.3, cex = 0.5, "Counts")

  axis(2, at = seq(0, 0.4, 0.1), labels = FALSE, tck = -0.025)
  axis(2, at = seq(0, 0.4, 0.2), las = 1, lwd = 0, line = -0.75,
       cex.axis = 0.5)
  axis(2, at = seq(0, 0.4, 0.05), labels = FALSE, tck = -0.01)
  mtext(side = 2, line = 1, cex = 0.5, "Frequency")

  par(fig = c(0, 1, 0.66, 0.72), mar = c(0, 0, 0, 0), new = TRUE)
  blank(xlim = c(0, 1), ylim = c(0, 1))
  text(0, 0.5, "(a)", xpd = TRUE, cex = 1.25)

  par(fig = c(0.25, 0.95, 0.65, 1), mar = c(1, 1, 1, 1), new = TRUE)
  blank(bty = "L", xlim = c(0.5, 12.5), ylim = c(-0.5, 15))


  polygon(c(0.5, 12.5, 13.26275, 13., 0.5),
          c(-1, -4, -4, -1, -1), col = rectcol1, border = NA,
         xpd = NA)
  rect(0.5, -1, 13.25, 20, col = rectcol1, border = NA)
  box(bty = "L", lwd = 2)

  for(i in 1:12){
    Hn <- unlist(mod3[[1]]$model$mcmc[,5+i])
    dens <- density(Hn)
    incl <- which(dens$x < 15 & dens$y > 0.01)
    fx <- i
    xxa <- fx + dens$y[incl] * 0.5
    xxb <- fx - dens$y[incl] * 0.5 
    xx <- c(xxa, xxb[length(xxb):1])
    xx2 <- c(xx, xx[1])
    violiny <- c(dens$x[incl], dens$x[incl][length(dens$x[incl]):1])
    violiny2 <- c(violiny, violiny[1])
    polygon(xx2, violiny2, border = grey(0.1, 1), col = grey(0.1, 1))
  }


  for(i in 1:12){
    Hn <- unlist(mod2[[1]]$model$mcmc[,3+i])
    dens <- density(Hn)
    incl <- which(dens$x < 15 & dens$y > 0.01)
    fx<-i
    xxa <- fx + dens$y[incl] * 0.5
    xxb <- fx - dens$y[incl] * 0.5 
    xx <- c(xxa, xxb[length(xxb):1])
    xx2 <- c(xx, xx[1])
    violiny <- c(dens$x[incl], dens$x[incl][length(dens$x[incl]):1])
    violiny2 <- c(violiny, violiny[1])
    polygon(xx2, violiny2, border = grey(0.4, 1), col = grey(0.4, 0.5))
}



  for(i in 1:12){
    Hn <- unlist(mod1[[1]]$model$mcmc[,2+i])
    dens <- density(Hn)
    incl <- which(dens$x < 15 & dens$y > 0.01)
    fx<-i
    xxa <- fx + dens$y[incl] * 0.5 
    xxb <- fx - dens$y[incl] * 0.5 
    xx <- c(xxa, xxb[length(xxb):1])
    xx2 <- c(xx, xx[1])
    violiny <- c(dens$x[incl], dens$x[incl][length(dens$x[incl]):1])
    violiny2 <- c(violiny, violiny[1])
    polygon(xx2, violiny2, border = grey(0.6, 1), col = grey(0.6, 0.5))
  }

  points(1:12, y[302:313], pch = 16, col = rgb(1, 1, 1, 0.8))
  points(1:12, y[302:313], lwd = 2, cex = 1.1, col = rgb(0, 0, 0, 0.9))

  mtext(side = 3, c("RW", "AR(1)", "cAR(1)"), at = c(4, 6, 8), line = 0.2,
        cex = 0.6, font = 2, col = c(grey(0.6), grey(0.4), grey(0.1)))

  axis(2, at = seq(0, 20, 5), labels = FALSE, tck = -0.05)
  axis(2, at = seq(0, 20, 1), labels = FALSE, tck = -0.03)
  axis(2, at = seq(0, 20, 10), las = 1, lwd = 0, line = -0.6,
       cex.axis = 0.6)
  mtext(side = 2, line = 1, cex = 0.55, 
        expression(paste(italic(C.), " ", italic(penicillatus ))))
  axis(1, at = seq(1, 12, 1), labels = FALSE, tck = -0.03)

  par(fig = c(0, 1, 0.9, 1), mar = c(0, 0, 0, 0), new = TRUE)
  blank(xlim = c(0, 1), ylim = c(0, 1))
  text(0.22, 0.8, "(d)", xpd = TRUE, cex = 1.25)

  dev.off()
}

fig2bottom <- function(models, nbins = 10, last_in = 500){

  mod1 <- models$mod1
  mod2 <- models$mod2 
  mod3 <- models$mod3

  tiff("fig2bottom.tiff", width = 6, height = 2.5, units = "in", res = 200)

  par(fig = c(0, 1, 0, 1))
  blank()

  eval_tab <- make_eval_tab(mod1, mod2, mod3)
  umods <- unique(eval_tab$model)
  numods <- length(umods)
  mod_ids <- umods
  mod_names <- c("RW", "AR(1)", "cAR(1)")
  mod_names2 <- gsub(" Ensemble", "", mod_names)
  PITs <- matrix(NA, nrow = numods, ncol = nbins)
  crpss <- rep(NA, numods)
  logss <- rep(NA, numods)
  for(i in 1:numods){
    incl <- which(eval_tab$model == umods[i] & eval_tab$destin <= last_in)
    eval_tabi <- eval_tab[incl,]
    PITcols <- grepl("PIT", colnames(eval_tabi))
    eval_tabiPIT <- eval_tabi[,PITcols]
    for(j in 1:ncol(eval_tabiPIT)){
      eval_tabiPIT[,j] <- as.numeric(eval_tabiPIT[,j])
    }
    PITsum <- apply(eval_tabiPIT, 2, sum, na.rm = TRUE)
    PITmean <- PITsum / nrow(eval_tabiPIT)
    PITs[i, ] <- PITmean
    crpss[i] <- mean(as.numeric(eval_tabi[,5]), na.rm = TRUE)
    logss[i] <- mean(as.numeric(eval_tabi[,6]), na.rm = TRUE)
  }

  max_crpss <- which.max(crpss)
  max_logss <- which.max(logss)

  for(i in 1:numods){
    x1 <- (i - 1) * 0.2
    x2 <- (i) * 0.2
    par(fig = c(x1, x2, 0.5, 1), mar = c(0.5, 0.5, 1.25, 0.25), new = TRUE)
    blank(bty = "L", ylim = c(0, 1.75), xlim = c(0.5, 10.5))
    points(PITs[i,], type = "h", lwd = 3)
    abline(h = 1, lty = 3, lwd = 2)
    text(1, 1.7, mod_names[i], cex = 0.75, adj = 0)
  }

  par(fig = c(0.6, 0.8, 0.5, 1), mar = c(1.5, 1.75, 1, 0.5), 
      new = TRUE)
  blank(bty = "L", ylim = c(-2.5, -1.5), xlim = c(0.4, 3.6))
  points(max_crpss, crpss[max_crpss], cex = 1.5, pch = 16, lwd = 0, 
         col = rgb(0.6, 0.6, 0.6, 0.6))
  points(max_crpss, crpss[max_crpss], cex = 1.5, pch = 1, lwd = 1, 
         col = rgb(0.6, 0.6, 0.6))
  points(crpss, pch = 16, cex = 0.9)

  text(1:3, -2.64, mod_names2, xpd = TRUE, srt = 55, cex = 0.5, adj = 1)
  mtext(side = 2, "RPS", line = 1.1, cex = 0.6)
  axis(side = 2, las = 1, cex.axis = 0.5, tck = -0.03, 
       at = seq(-2.5, -1.5, 0.1), labels = FALSE)
  axis(side = 2, las = 1, cex.axis = 0.5, tck = -0.05, 
       at = seq(-2.5, -1.5, 0.5), labels = FALSE)
  axis(side = 2, las = 1, cex.axis = 0.5, tck = -0.05, line = -0.7, 
       at = seq(-2.5, -1.5, 0.5), lwd = 0)

  par(fig = c(0.8, 1, 0.5, 1), mar = c(1.5, 1.75, 1, 0.5), new = TRUE)

  blank(bty = "L", ylim = c(-2.5, -2), xlim = c(0.4, 3.6))
  points(max_logss, logss[max_logss], cex = 1.5, pch = 16, lwd = 0, 
         col = rgb(0.6, 0.6, 0.6, 0.6))
  points(max_logss, logss[max_logss], cex = 1.5, pch = 1, lwd = 1, 
         col = rgb(0.6, 0.6, 0.6))
  points(logss, pch = 16, cex = 0.9)

  mtext(side = 2, "Log Score", line = 1.15, cex = 0.6)
  axis(side = 2, las = 1, cex.axis = 0.5, tck = -0.05, labels = FALSE)
  axis(side = 2, las = 1, at = c(-2.7, -2.5, -2.3, -2.1),
        cex.axis = 0.5, tck = -0.05, line = -0.7, lwd = 0)

  text(1:3, -2.57, mod_names2, xpd = TRUE, srt = 55, cex = 0.5, adj = 1)
  par(fig = c(0, 1, 0.92, 1), mar = c(0, 0, 0, 0), new = TRUE)
  blank(xlim = c(0, 1), ylim = c(0, 1))
  text(0, 0.45, "(b)", xpd = TRUE, cex = 1.25)

  eval_tabft <- eval_tab[eval_tab$origin == 500,]
  eval_tabft_PIT <- matrix(NA, nrow = numods, ncol = nbins)
  eval_tabft_c <- rep(NA, numods)
  eval_tabft_l <- rep(NA, numods)

  for(i in 1:numods){
    incl <- which(eval_tabft$model == umods[i])
    PITcol <- grepl("PIT", colnames(eval_tabft))
    x <- eval_tabft[incl, PITcol]
    for(j in 1:nbins){
      x[, j] <- as.numeric(x[,j])
    }

    eval_tabft_PIT[i, ] <- apply(x, 2, sum) / nrow(x)
    eval_tabft_c[i] <- mean(as.numeric(eval_tabft$crps[incl])) 
    eval_tabft_l[i] <- mean(as.numeric(eval_tabft$logs[incl])) 
  }

  max_eval_tabft_c <- which.max(eval_tabft_c)
  max_eval_tabft_l <- which.max(eval_tabft_l)

  for(i in 1:numods){
    x1 <- (i - 1) * 0.2
    x2 <- (i) * 0.2
    par(fig = c(x1, x2, 0, 0.5), mar = c(0.5, 0.5, 1.25, 0.25), new = TRUE)
    blank(bty = "L", ylim = c(0, 3.5), xlim = c(0.5, 10.5))
    points(eval_tabft_PIT[i,], type = "h", lwd = 3)
    abline(h = 1, lty = 3, lwd = 2)
    text(1, 3.4, mod_names[i], cex = 0.75, adj = 0, xpd = TRUE)
  }

  par(fig = c(0.6, 0.8, 0, 0.5), mar = c(1.5, 1.75, 1, 0.5), new = TRUE)
  blank(bty = "L", ylim = c(-3.1, -0.9), xlim = c(0.5, 3.5))
  points(max_eval_tabft_c, eval_tabft_c[max_eval_tabft_c], cex = 1.5, 
         pch = 16, lwd = 0, col = rgb(0.6, 0.6, 0.6, 0.6))
  points(max_eval_tabft_c, eval_tabft_c[max_eval_tabft_c], cex = 1.5, 
         pch = 1, lwd = 1, col = rgb(0.6, 0.6, 0.6))
  points(eval_tabft_c, pch = 16, cex = 0.9)
  text(1:3, -3.4, mod_names2, xpd = TRUE, srt = 55, cex = 0.5, adj = 1)
  mtext(side = 2, "RPS", line = 1.1, cex = 0.6)

  axis(side = 2, las = 1, cex.axis = 0.5, tck = -0.03, 
       at = seq(-3, 1, 0.5), labels = FALSE)
  axis(side = 2, las = 1, cex.axis = 0.5, tck = -0.05, 
       at = seq(-3, 1, 1), labels = FALSE)
  axis(side = 2, las = 1, cex.axis = 0.5, tck = -0.05, line = -0.7, 
       at = seq(-3, 1, 1), lwd = 0)
  
  par(fig = c(0.8, 1, 0, 0.5), mar = c(1.5, 1.75, 1, 0.5), new = TRUE)
  blank(bty = "L", ylim = c(-2.7, -2.2), xlim = c(0.5, 3.5))
  points(max_eval_tabft_l, eval_tabft_l[max_eval_tabft_l], cex = 1.5, 
         pch = 16, lwd = 0, col = rgb(0.6, 0.6, 0.6, 0.6))
  points(max_eval_tabft_l, eval_tabft_l[max_eval_tabft_l], cex = 1.5, 
         pch = 1, lwd = 1, col = rgb(0.6, 0.6, 0.6))
  points(eval_tabft_l, pch = 16, cex = 0.9)
  text(1:3, -2.76, mod_names2, xpd = TRUE, srt = 55, cex = 0.5, adj = 1)
  mtext(side = 2, "Log Score", line = 1.15, cex = 0.6)
  axis(side = 2, las = 1, cex.axis = 0.5, tck = -0.05, labels = FALSE)
  axis(side = 2, las = 1, at = c(-2.7, -2.5, -2.3),
        cex.axis = 0.5, tck = -0.05, line = -0.7, lwd = 0)

  par(fig = c(0, 1, 0.44, 0.52), mar = c(0, 0, 0, 0), new = TRUE)
  blank(xlim = c(0, 1), ylim = c(0, 1))
  text(0, 0.45, "(c)", xpd = TRUE, cex = 1.25)
  dev.off()
}



make_figures <- function(data_set, models){
  fig1top()
  fig1bottom()
  figA1()
  fig2top(data_set, models)
  fig2bottom(models)
}

