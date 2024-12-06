create_maps <- function(
    t,                 # Time vector or index
    f = NULL,          # Signal vector (optional; defaults to `t`)
    maxtau,            # Maximum value for tau (frequency parameter)
    M = "mean",        # Method for central tendency (default: "mean")
    V = "var",         # Method for variability (default: "var")
    inter_method = "linear" # Interpolation method (default: "linear")
) {
  # If no signal vector is provided, set `f` as `t` and generate `t` as an index
  if (is.null(f)) {
    f <- t
    t <- 1:length(f)
  }
  
  # Initialize matrices for Vmap and Mmap
  VM <- matrix(0, nrow = length(f), ncol = maxtau)
  MM <- VM
  
  # Compute Vmap and Mmap for each tau value
  for (tau in 2:maxtau) {
    out <- ebt(t, f, tau = tau, M, V, inter_method)
    VM[, tau] <- out$V
    MM[, tau] <- out$M
  }
  
  # Return the results as a list
  list(
    t = t,      # Time vector
    vmap = VM,  # Variability map
    mmap = MM   # Mean map
  )
}