#' De-noise a numeric vector
#'
#' @description De-noise a numeric vector by variable mode by Variational Mode
#'   Decomposition. Find optimum number of modes, do the decomposition (using at
#'   least three modes) and return the first mode, containing the quasi-static
#'   component.
#'
#' @noRd
denoise_vector.internal <- function(data,
                                    alpha = 2000,
                                    K = 100,
                                    verbose = FALSE) {
  # de-noise with the optimum number of modes
  data.vmd <- VMDecomp::vmd(
    data,
    alpha = alpha,
    tau = 0, # noise-tolerance (no strict fidelity enforcement)
    K = K,
    DC = TRUE, # a DC part imposed
    init = 1, # initialize omegas uniformly
    tol = 1e-6,
    verbose = verbose
  )

  # return all modes
  return(data.vmd)
}

#' De-noise a numeric vector
#'
#' @description De-noise a numeric vector by variable mode by Variational Mode
#'   Decomposition. Find optimum number of modes, do the decomposition (using at
#'   least three modes) and return the first mode, containing the quasi-static
#'   component.
#'
#' @noRd
denoise_vector <- function(data,
                           alpha = 2000,
                           K = 100,
                           verbose = FALSE) {
  # de-noise
  vmd <- denoise_vector.internal(data, alpha, K, verbose)

  # return the first mode (DC-component) only
  return(vmd$u[, 1])
}

# FIXME
denoise_vector.tst <- function(data,
                               alpha = 2000,
                               K = 100,
                               verbose = FALSE) {
  # de-noise
  vmd <- denoise_vector.internal(data, alpha, K, verbose)

  # return the first mode (DC-component) only
  return(vmd)
}

#' De-noise a numeric vector
#'
#' @description De-noise a numeric vector by variable mode by Variational Mode
#'   Decomposition. Find optimum number of modes, do the decomposition (using at
#'   least three modes) and return the first mode, containing the quasi-static
#'   component.
#'
#' @noRd
denoise_vector.dc <- function(data,
                              alpha = 2000,
                              K = 100,
                              verbose = FALSE) {
  # de-noise
  vmd <- denoise_vector.internal(data, alpha, K, verbose)

  # return the first mode (DC-component) only
  return(vmd$u[, 1])
}
