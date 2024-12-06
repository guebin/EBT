extract_signal <- function(
    time_vector,       # The time vector of the input data
    signal,            # The signal vector of the input data
    tau,               # Parameter for signal filtering
    M = "mean",        # Method for estimating the central tendency (default: "mean")
    V = "sd",          # Method for estimating variability (default: "sd")
    tol = 0.05,        # Convergence threshold (stopping condition for iterations)
    iter = 100,        # Maximum number of iterations
    interpolation = "linear" # Method for interpolation (default: "linear")
) {
  # Initialize with the original signal for filtering
  prev_filtered_signal <- signal
  
  # Perform the initial filtering using EBT (Envelope Band Transformation)
  current_filtered_signal <- ebt(
    signal, tau = tau, M = M, V = V, interpolation = interpolation
  )$M
  
  # Compute the residual signal after the first filtering
  residual_signal <- signal - current_filtered_signal
  
  iteration <- 0  # Initialize the iteration counter
  
  # Continue iterating until convergence or the maximum number of iterations is reached
  while ((max(abs(prev_filtered_signal - current_filtered_signal)) > tol) & 
         (iteration < iter)) {
    prev_filtered_signal <- current_filtered_signal  # Store the signal from the previous iteration
    
    # Reapply filtering to the residual signal
    current_filtered_signal <- ebt(
      residual_signal, tau = tau, M = M, V = V, interpolation = interpolation
    )$M
    
    # Update the residual signal
    residual_signal <- residual_signal - current_filtered_signal
    
    # Increment the iteration counter
    iteration <- iteration + 1
  }
  
  # Return the final residual signal
  final_residual_signal <- residual_signal
}