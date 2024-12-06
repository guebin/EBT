ebt <- function(
    t,                  # Time vector
    f = NULL,           # Signal vector (optional; defaults to `t`)
    tau,                # Parameter for filtering
    mfunc = "mean",     # Method for calculating the central tendency ("mean", "median", or "volume")
    vfunc = "var",      # Method for calculating variability ("var" or "volume")
    inter_method = "linear" # Interpolation method ("linear" or "cubic")
) {
  # If no signal vector is provided, set `f` as `t` and generate `t` as an index
  if (is.null(f)) {
    f <- t
    t <- 1:length(f)
  }
  
  # Save original time and signal vectors
  tsave <- t
  fsave <- f
  
  # Extend signal with mirrored values for boundary handling
  f <- c(rep(f[1], tau * 2), f, rep(f[length(f)], tau * 2))
  t <- 1:length(f)
  len <- length(f)
  
  # Initialize variables
  sampled_index <- list()
  missing_index <- list()
  band <- matrix(0, nrow = len, ncol = tau)
  
  # Perform interpolation based on the specified method
  if (inter_method == "cubic") {
    for (eta in 1:tau) {
      sampled_index[[eta]] <- seq(from = eta, to = len, by = tau)
      if (sampled_index[[eta]][1] != 1) sampled_index[[eta]] <- c(1, sampled_index[[eta]])
      if (sampled_index[[eta]][length(sampled_index[[eta]])] != len) {
        sampled_index[[eta]] <- c(sampled_index[[eta]], len)
      }
      missing_index[[eta]] <- (1:len)[-sampled_index[[eta]]]
      result <- smooth.spline(t[sampled_index[[eta]]], f[sampled_index[[eta]]], spar = 0, all.knots = TRUE)
      band[sampled_index[[eta]], eta] <- result$y
      band[missing_index[[eta]], eta] <- predict(result, t[missing_index[[eta]]])$y
    }
  } else if (inter_method == "linear") {
    for (eta in 1:tau) {
      sampled_index[[eta]] <- seq(from = eta, to = len, by = tau)
      if (sampled_index[[eta]][1] != 1) sampled_index[[eta]] <- c(1, sampled_index[[eta]])
      if (sampled_index[[eta]][length(sampled_index[[eta]])] != len) {
        sampled_index[[eta]] <- c(sampled_index[[eta]], len)
      }
      missing_index[[eta]] <- (1:len)[-sampled_index[[eta]]]
      band[, eta] <- .lin_impute(t, f, missing_index[[eta]])
    }
  }
  
  # Compute upper and lower bounds
  U <- apply(band, 1, max)
  L <- apply(band, 1, min)
  
  # Compute mean, median, and volume-based central tendency
  M1 <- apply(band, 1, mean)
  M2 <- apply(band, 1, median)
  M3 <- (L + U) / 2
  
  # Compute variability metrics
  V1 <- U - L
  V2 <- apply(band, 1, var) * (tau - 1) / tau
  
  # Select central tendency based on `mfunc`
  if (mfunc == "mean") {
    M <- M1
  } else if (mfunc == "median") {
    M <- M2
  } else if (mfunc == "volume") {
    M <- M3
  }
  
  # Select variability metrics based on `vfunc`
  if (vfunc == "volume") {
    V <- V1
  } else if (vfunc == "var") {
    V <- V2
    L <- M - V
    U <- M + V
  }
  
  # Define index for the output
  index <- (2 * tau + 1):(length(f) - 2 * tau)
  
  # Filter sampled indices to match the length of the original signal
  sampled_index <- lapply(sampled_index, function(x) x[x <= length(fsave)])
  
  # Return results as a list
  list(
    t = tsave,            # Original time vector
    f = fsave,            # Original signal vector
    V = V[index],         # Variability vector
    L = L[index],         # Lower bounds
    U = U[index],         # Upper bounds
    M = M[index],         # Central tendency vector
    tau = tau,            # Parameter tau
    band = band[index, ], # Interpolated band matrix
    knot = sampled_index  # Sampled indices
  )
}