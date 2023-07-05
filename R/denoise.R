#' De-noise a numeric vector
#'
#' @description De-noise a numeric vector by variable modeby Variational Mode
#'   Decomposition. Find optimum number of modes, do the decomposition (using at
#'   least three modes) and return the first mode, containing the quasi-static
#'   component.
#'
#' @noRd
denoise_vector <- function(data, verbose = FALSE) {
  default_vmd_params <- list(
    alpha = 2000, # moderate bandwidth constraint
    tau = 0, # noise-tolerance (no strict fidelity enforcement)
    DC = TRUE, # a DC part imposed
    init = 1, # initialize omegas uniformly
    tol = 1e-6
  )

  optimum_k <- VMDecomp::estimate_k_modes(
    signal_1d = data,
    cor_thresh = 0.1,
    default_vmd_params = default_vmd_params,
    min_K = 2,
    seed = 1,
    verbose = verbose
  )

  # de-noise with the optimum number of modes
  data.vmd <- VMDecomp::vmd(data,
    alpha = default_vmd_params$alpha,
    tau = default_vmd_params$tau,
    K = optimum_k,
    DC = default_vmd_params$DC,
    init = default_vmd_params$init,
    tol = default_vmd_params$tol,
    verbose = verbose
  )

  # return the first mode (DC-component) only
  return(data.vmd$u[, 1])
}
