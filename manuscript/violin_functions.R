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
  nvalues <- length(xvals)
  yvals <- rep(NA, nvalues)
  den <- density(x)
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
      yvals[i] <- (yval_val_1 * xval_ratio_1 + yval_val_2 * xval_ratio_2) / 2
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
