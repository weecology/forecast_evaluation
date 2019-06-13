




dens












# so this means we need to have a function that helps set the classes of
#  the functions so it's clear how the methods get dispatched

cdf <- function(x){ppois(x, 1)}
class(cdf) <- c("cdf", "function")

pmf <- function(x){dpois(x, 1)}
class(pmf) <- c("pmf", "function")

emcdf <- ecdf(rpois(1e6, 1))
class(emcdf)

exem <- rpois(1e6, 1)
class(exem)


# this is the foundation of fig 3

p <- function(pars, sym){
  pars[[sym, exact = TRUE]]
}

pchoice <- function(p, def1, def2 = NULL){
  if (length(p)){
    p 
  } else if (length(def1)){ 
    def1 
  } else{
    def2
  }
}






draw_violin <- function(location, vals, type = "l", 
                        side = "both"){




}




draw_violin_line <- function(){

}


draw_violin_hist <- function(){






}

# working in and around here... make it so we can actually start drawing
# some basic version and then iterate off it!

x <- rnorm(1e6, 0, 1)

blank(bty = "L")
points(1,1)
violin(xx, 1)

& is.null(names(location))

# we'll get to this, but first, keep it simple



violin <- function(x, ...){
  UseMethod("violin", x)
}

# start with a vector of numbers

violin.numeric <- function(x, ...){

}




violin.cdf <- function(x, ...){

}
violin.pmf <- function(x, ...){

}
violin.ecdf <- function(x, ...){

}
violin.integer <- function(x, ...){


}












# not sure...



# defines the sample space of a variable or function

support <- function(x, ...){
  UseMethod("support", x)
}

support.integer <- function(x, ...){
  as.integer(names(table(x)))
}

support.numeric <- function(x, step_size = NULL, step_density = 10, 
                            nsteps = NULL, ...){
  if(is.null(step_size) & is.null(nsteps)){
    step_size <- 0
    step_width <- dist(abs(x))
    if (length(step_width) > 0 ){
      min_step_width <- max(c(min(step_width), min(abs(x))))
      step_size <- min_step_width / step_density
    }
  } else if(!is.null(nsteps)){
    step_size <- (max(x) - min(x)) / nsteps
  }
  seq(min(x), max(x), step_size)
}



# based off of boxplot.stats
#  this works fine for an empirical distribution sample but what if
#  we have a defined distribution?
#   and a distribution defined using a cdf or a pd/mf
# decision points are: 
#    definition (empirical or function, and then cdf or pd/mf)
#    support (binomial, continuous, integer, pos/neg, categorical)

# what fivenum does needs to be generalized for sure

violin_stats <- function(x, coef = 1.5, do.conf = TRUE, do.out = TRUE){
  if (coef < 0) {
    stop("'coef' must not be negative")
  }
  nna <- !is.na(x)
  n <- sum(nna)
  stats <- fivenum(x, na.rm = TRUE)
  iqr <- diff(stats[c(2, 4)])
  if (coef == 0) {
    do.out <- FALSE
  } else {
    out <- if (!is.na(iqr)) {
    x < (stats[2L] - coef * iqr) | x > (stats[4L] + coef * 
                iqr)
        }
       else !is.finite(x)
        if (any(out[nna], na.rm = TRUE)) 
            stats[c(1, 5)] <- range(x[!out], na.rm = TRUE)
  }
  conf <- NULL
  if (do.conf){ 
    conf <- stats[3L] + c(-1.58, 1.58) * iqr/sqrt(n)
  }
  out <- NULL
  if (do.out){ 
   out <- x[out & nna]
  }
    list(stats = stats, n = n, conf = conf, out = out)
}





set.seed(1234)
minlam <- 1
maxlam <- 7
nlams <- 1e3
lams <- seq(minlam, maxlam, length.out = nlams)
steps <- 1:nlams
x <- rpois(nlams, lambda = lams)
n_bins <- 10



par(mfrow=c(6,3), mar=c(4,4,1,1))

# each of these cases really just needs a predictive cdf, so simplify
# and function it!

# also code something to get violin layovers for plots more easily!
# and to differentiate between discrete and continuous violins


# working here right now on turning the "on point" forecast into good
# basic forecasting figure stuff!


# on point

minstep <- 990
maxstep <- 1000
buff <- 0.5
minp <- 1e-3
vwidth <- 0.5
xrange <- 0:25
vxs <- 996:1000
pdfun <- function(x, lam){dpois(x, lam)}


ts_panel <- function(pdfun, lams[vx]
  blank(bty = "L", xlim = c(minstep - buff, maxstep + buff), ylim = c(0, 25))

  nvxs <- length(vxs)
  for(j in 1:nvxs){
    vx <- vxs[j]

    nxr <- length(xrange)
    mass <- pdfun(xrange, ...)
    vwidth2 <- vwidth / 2
    vwidthmult <- vwidth2 / max(mass)
    vwidths <- mass * vwidthmult
    vwidths[vwidths < minp] <- NA
    vxs1 <- vx + vwidths
    vxs2 <- vx - vwidths

    for(i in 1:nxr){
      vxs12 <- c(vxs1[i], vxs2[i])
      vys12 <- rep(xrange[i], 2)
      points(vxs12, vys12, type = "l", lwd = 4, col = rc)
    }
  }
  points(steps[steps >= minstep], x[steps >= minstep], cex = 1.5, lwd = 2)
}



cdf <- function(x, lam){ppois(x, lambda = lam)}

nrPITs(x, cdf, lams)


diff(cdf(0:20, lams[nlams]))


nrPITs <- function(x, cdf, ...){
  nx <- length(x)
  pits <- matrix(NA, nrow = nx, ncol = n_bins)
  for(i in 1:nlams){
    pits[i,] <- nrPIT(x[i], cdf, n_bins, lams[i])
  }
  pits
}



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
  cdf <- function(x){vals <- c(rpois(5e4, max(c(0,lam-2))), rpois(5e4, lam+2))
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

