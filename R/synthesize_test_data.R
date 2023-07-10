#' @title Synthesize Test Data
#'
#' @description Synthesize Test Data by Ramberg-Osgood equation.
#'
#' @details Use the Ramberg-Osgood equation in the formulation of Hill to create
#'   synthetic data similar to quasi-static mechanical testing. Parameters K and
#'   n of the Ramberg-Osgood equation are calculated from the yield-strength
#'   (typically at 0.2% offset) and a second point in the post-yield region
#'   (e.g. the point of ultimate stress) as proposed by Hill.
#'
#' @note As random values are drawn, set a [`random seed`][base::RNG] beforehand
#'   to get reproducible results.
#'
#' @references Hill, H. N. (1944). Determination of stress-strain relations from
#'   "offset" yield strength values. Aluminum Co of America Pittsburgh Pa.
#'
#' @references Ramberg, W., & Osgood, W. R. (1943). Description of Stress-Strain
#'   Curves by Three Parameters; National Advisory Committee for Aeronautics
#'   Technical Note. NACA-TN-902.
#'
#' @references Rice, R. C., Jackson, J. L., Bakuckas, J., & Thompson, S. (2003).
#'   Metallic Materials Properties Development and Standardization (MMPDS)
#'   (DOT/FAA/AR-MMPDS-01).
#'
#' @param slope Slope in the linear region.
#'
#' @param yield.y y-value at yield point.
#'
#' @param yield.xp Plastic deformation at yield point. Defaults to 0.002.
#'
#' @param ultimate.y Maximum y-value in the post-linear region.
#'
#' @param ultimate.x Maximum x-value in the post-linear region.
#'
#' @param offset Value of y-offset.
#'
#' @param toe.start.y Intersection of toe-region with y-axis (before adding an
#'   offset).
#'
#' @param toe.start.slope Initial slope of toe region.
#'
#' @param toe.end.y End of toe region.
#'
#' @param enob.x Effective number of bits for the synthetic data x-range. Will
#'   determine the number of points in the returned data (i.e. 2^enob.x). Also
#'   used for adding quantization noise.
#'
#' @param enob.y Effective number of bits for the synthetic data y-range. Used
#'   for adding quantization noise.
#'
#' @param enob.x_FS Effective number of bits for the full-scale x-range. Using
#'   to determine level of quantization/data-noise in x-values.
#'
#' @param enob.y_FS Effective number of bits for the full-scale y-range. Using
#'   to determine level of quantization/data-noise in y-values.
#'
#' @param enob.x_noise Add noise to x-data. Give the effective number of bits
#'   for the full-scale x-range.
#'
#' @param enob.y_noise Add noise to y-data. Give the effective number of bits
#'   for the full-scale y-range.
#'
#' @param x.name Name for x-values. Defaults to `"strain"`.
#'
#' @param y.name Name for y-values. Defaults to `"stress"`.
#'
#' @param x.unit Unit for x-values. Can be `NULL`.
#'
#' @param y.unit Unit for y-values. Can be `NULL`. Defaults to `"MPa"`.
#'
#' @returns A data.frame with the synthetic data. Units (when provided) are
#'   given as variable labels.
#'
#' @examples
#' # Synthesize a test record resembling Al 6060 T66
#' # (Values according to Metallic Material Properties
#' # Development and Standardization (MMPDS) Handbook)
#' Al_6060_T66 <- synthesize_test_data(
#'   slope = 68000,
#'   yield.y = 160,
#'   ultimate.y = 215,
#'   ultimate.x = 0.091
#' )
#'
#' plot(
#'   x = Al_6060_T66$strain, y = Al_6060_T66$stress,
#'   type = "l",
#'   xlab = "strain", ylab = "stress (in MPa)"
#' )
#'
#' @export
synthesize_test_data <- function(slope, yield.y, yield.xp = 0.002,
                                 ultimate.y, ultimate.x, offset = 0,
                                 toe.start.y = 0, toe.start.slope = slope, toe.end.y = 0,
                                 enob.x = 14, enob.y = 14, enob.x_FS = 16, enob.y_FS = 16,
                                 enob.x_noise = 0, enob.y_noise = 0,
                                 x.name = "strain", y.name = "stress",
                                 x.unit = NULL, y.unit = "MPa") {
  # calculate plastic strain at ultimate stress
  ef <- ultimate.x - ultimate.y / slope
  # calculate Ramberg-Osgood coefficient (inverse of strain hardening exponent)
  n <- log(ef / yield.xp) / log(ultimate.y / yield.y)

  # get a splinefun to predict stress/force from sequence of strain/displacement values
  synthetic_data.sfun <- data.frame("y" = seq(0, ultimate.y, length.out = 2^enob.x)) %>%
    dplyr::mutate("x" = .data$y / slope + yield.xp * (.data$y / yield.y)^n) %>%
    {
      stats::splinefun(.$x, .$y,
        method = "monoH.FC",
        ties = list("ordered", mean)
      )
    }

  # get a splinefun for creating synthetic data with added toe region
  # (strain/displacement...) with added toe region
  synthetic_data.toed.sfun <- data.frame("x" = seq(from = 0, to = ultimate.x, length.out = 2^enob.x)) %>%
    dplyr::mutate("y" = synthetic_data.sfun(.data$x)) %>%
    {
      # add a toe region (if possible)
      synthetic_data <- .
      synthetic_data.toe <- synthetic_data %>% dplyr::filter(.data$y <= toe.end.y)
      synthetic_data.post_toe <- synthetic_data %>% dplyr::filter(.data$y > toe.end.y)

      # does it make any sense to add a toe region to the synthetic data?
      if (nrow(synthetic_data.toe) > 2) {
        toe.sfun <- data.frame(
          "x" = c(
            seq(-23 * synthetic_data.toe[2, "x"], 0,
              length.out = 24
            ),
            utils::tail(synthetic_data.toe, 1)$x,
            utils::head(synthetic_data.post_toe, 23)$x
          ),
          "y" = c(
            seq(-23 * synthetic_data.toe[2, "x"], 0,
              length.out = 24
            ) * toe.start.slope + toe.start.y,
            utils::tail(synthetic_data.toe, 1)$y,
            utils::head(synthetic_data.post_toe, 23)$y
          )
        ) %>%
          {
            stats::splinefun(.$x, .$y,
              method = "monoH.FC",
              ties = list("ordered", mean)
            )
          }

        synthetic_data.toe %>%
          dplyr::mutate("y" = toe.sfun(.data$x)) %>%
          dplyr::bind_rows(synthetic_data.post_toe)
      } else {
        if (toe.end.y > 0) {
          warning("Toe-region is not feasible. (Too little data...)")
        }
        synthetic_data
      }
    } %>% {
      stats::splinefun(.$x, .$y,
                       method = "monoH.FC",
                       ties = list("ordered", mean))
    }

  # synthesize data starting with a regular sequence of x-values
  # (strain/displacement...) with added toe region
  # start with the sequence of x-values
  data.frame("x" = seq(from = 0, to = ultimate.x, length.out = 2^enob.x)) %>% {
    # add noise to x-data
    synthetic_data <- .
    synthetic_data.nrow <- nrow(synthetic_data)
    # maximum value "recorded" x-data
    x.max <- max(synthetic_data$x)
    # hypothetical x-value at FS
    x.max.FS <- x.max * 2^(enob.x_FS - enob.x)
    # x-value worth at 1 bit of FS
    x.resolution <- x.max.FS * 2^-enob.x_FS
    # sd of noise
    x.noise.sd <- (2^enob.x_noise - 1) * x.resolution

    # add noise from a normal distribution
    synthetic_data %>%
      dplyr::mutate("x" = .data$x + stats::rnorm(synthetic_data.nrow, sd = x.noise.sd))
  } %>% {
    # add quantization noise
    synthetic_data <- .

    # maximum value "recorded" x-data
    x.max <- max(synthetic_data$x)
    # hypothetical x-value at FS
    x.max.FS <- x.max * 2^(enob.x_FS - enob.x)
    # x-value worth at 1 bit of FS
    x.resolution <- x.max.FS * 2^-enob.x_FS

    synthetic_data %>%
      dplyr::mutate("x" = round(.data$x / x.resolution, 0) * x.resolution)
  } %>%
    # order the x-data (in case adding noise mixed something up)
    dplyr::arrange(.data$x) %>%
    # use the splinefun with the toe region to synthesize y-data
    dplyr::mutate("y" = synthetic_data.toed.sfun(.data$x)) %>% {
      # add noise to y-data
      synthetic_data <- .
      synthetic_data.nrow <- nrow(synthetic_data)
      # maximum value "recorded" y-data
      y.max <- max(synthetic_data$y)
      # hypothetical y-value at FS
      y.max.FS <- y.max * 2^(enob.y_FS - enob.y)
      # y-value worth at 1 bit of FS
      y.resolution <- y.max.FS * 2^-enob.y_FS
      # sd of noise
      y.noise.sd <- (2^enob.y_noise - 1) * y.resolution

      # add noise from a normal distribution
      synthetic_data %>%
        dplyr::mutate("y" = .data$y + stats::rnorm(synthetic_data.nrow, sd = y.noise.sd))
    } %>%
    # add offset to y-data
    dplyr::mutate("y" = .data$y + offset) %>% {
      # add quantization noise
      synthetic_data <- .

      # maximum value "recorded" y-data
      y.max <- max(synthetic_data$y)
      # hypothetical y-value at FS
      y.max.FS <- y.max * 2^(enob.y_FS - enob.y)
      # y-value worth at 1 bit of FS
      y.resolution <- y.max.FS * 2^-enob.y_FS

      synthetic_data %>%
        dplyr::mutate("y" = round(.data$y / y.resolution, 0) * y.resolution)
    } %>%
    # set labels
    labelled::set_variable_labels(.labels = list(x = x.unit, y = y.unit)) %>%
    #set names
    magrittr::set_names(c(x.name, y.name)) %>%
    return()
}
