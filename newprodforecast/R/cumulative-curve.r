cumulative_curve <- function(params, time, acv) {  
  # Internal function to calculate cumulative trial curve as a function of time
  y <- if(is.null(params$lambda)) with(params, 1 - (alpha / (alpha + time))^r) else with(params, 1 - exp(-lambda*time))
  if (!is.null(params$p0)) 
    y <- params$p0 * y
  if (!is.null(params$gamma)) 
    y <- (acv ^ params$gamma) * y
  y
}
