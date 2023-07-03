#' First Data Quality Metric
#'
#' @note
#' Refer to chapter 6.4 in ASTM E3076-18
#'
#' @description
#' Check for noise and resolution of the normalized data
#' pass the normalized test-data and make sure that x and y are in columns "x.normalized" and "y.normalized"
#'
#' @returns a list with noise metrics for normalized x and y
#'
#' @noRd
check_data_quality.noise <- function(data.normalized, verbose = F, warn = T) {

  # calculate noise for x and y
  normalized.noise <- data.normalized %>%
    dplyr::mutate(delta.x = dplyr::lead(x.normalized) - x.normalized) %>%
    dplyr::mutate(delta.y = dplyr::lead(y.normalized) - y.normalized) %>%
    utils::head(-1) %>%
    dplyr::mutate(delta.x.r = delta.x - mean(delta.x)) %>%
    dplyr::mutate(delta.y.r = delta.y - mean(delta.y)) %>% {
      data.frame(normalized.noise.x = stats::sd(.$delta.x.r),
                 normalized.noise.y = stats::sd(.$delta.y.r))
    }

  normalized.noise.x <- normalized.noise[1,"normalized.noise.x"] %>% unlist(T,F)
  normalized.noise.y <- normalized.noise[1,"normalized.noise.y"] %>% unlist(T,F)
  relative.noise.x <- normalized.noise.x/0.005
  relative.noise.y <- normalized.noise.y/0.005
  passed.check.noise.x <- ifelse(relative.noise.x <= 1, T, F)
  passed.check.noise.y <- ifelse(relative.noise.y <= 1, T, F)

  # print messages
  if(verbose) {
    message("  Data Quality Metric: Noise")
    message(paste0("      normalized x noise:      ", normalized.noise.x))
    message(paste0("      relative x noise:        ", relative.noise.x))
    message(paste0("      normalized y noise:      ", normalized.noise.y))
    message(paste0("      relative y noise:        ", relative.noise.y, "\n"))
  }


  # Warn in case of excessive noise
  if(warn) {
    if(!passed.check.noise.x) {
      warning(paste0("excessive relative noise for x-values (", round(relative.noise.x*100, 1), " %)\n"), call. = F)
    }

    if(!passed.check.noise.y) {
      warning(paste0("excessive relative noise for y-values (", round(relative.noise.y*100, 1), " %)\n"), call. = F)
    }
  }

  # passed all checks for noise?
  passed.check.noise <- passed.check.noise.x && passed.check.noise.y

  # return quality metrics
  list(
    "passed.check" = passed.check.noise,
    "x" = list(
      "Relative_noise" = relative.noise.x,
      "passed.check" = passed.check.noise.x
    ),
    "y" = list(
      "Relative_noise" = relative.noise.y,
      "passed.check" = passed.check.noise.y
    )
  )
}

#' Second Data Quality Metric
#'
#' @note
#' Refer to chapter 6.5 in ASTM E3076-18
#'
#' @description
#' Check for noise and resolution of the normalized data
#' pass the normalized test-data and make sure that x and y are in columns "x.normalized" and "y.normalized"
#'
#' @returns a list with resolution metrics for x and y
#'
#' @noRd
check_data_quality.resolution <- function(data.normalized,
                                          verbose = F,
                                          showPlots = F,
                                          savePlots = F,
                                          warn = T) {

  # The optimum digital resolution, δ (delta.optimal), for the normalized data set is:
  delta.optimal <- 2^-12

  # calculate resolution for x and y
  normalized.resolution <- data.normalized %>%
    dplyr::mutate(delta.x = dplyr::lead(x.normalized) - x.normalized) %>%
    dplyr::mutate(delta.y = dplyr::lead(y.normalized) - y.normalized) %>%
    utils::head(-1) %>%
    dplyr::mutate(delta.delta.x = abs(dplyr::lead(delta.x) - delta.x)) %>%
    dplyr::mutate(delta.delta.y = abs(dplyr::lead(delta.y) - delta.y)) %>%
    utils::head(-1) %>% {
      data.frame(delta.delta.x = .$delta.delta.x,
                 delta.delta.y = .$delta.delta.y)
    }

  # number of samples after calculation (N-2)
  samples <- normalized.resolution %>% nrow

  # Analyze x-resolution
  max.delta.delta.x <- normalized.resolution$delta.delta.x %>% max
  n_bins.delta.delta.x <- (max.delta.delta.x / delta.optimal) %>% ceiling
  bin.breaks.x <- seq(-.5,n_bins.delta.delta.x + 1, 1)

  sdar.hist.x <- graphics::hist(normalized.resolution$delta.delta.x / delta.optimal,
                                breaks = bin.breaks.x, plot = showPlots,
                                main = "Histogram for assessment of x-Resolution",
                                xlab = "\u2206\u2206xi (binned by \u03B4)", warn.unused = F)

  z.delta.delta.x <- sdar.hist.x[["counts"]] %>%
    utils::tail(-1) %>%
    which.max
  relativeResolution.x <- z.delta.delta.x/3
  passed.check.resolution.x <- ifelse(relativeResolution.x <= 1, T, F)
  passed.check.binz.ratio.x <- ifelse((sdar.hist.x[["counts"]][z.delta.delta.x+1])/samples > 0.25, F, T)
  passed.check.bin0.ratio.x <- ifelse((sdar.hist.x[["counts"]][1])/samples > 0.25, F, T)

  # Analyze y-resolution
  max.delta.delta.y <- normalized.resolution$delta.delta.y %>% max
  n_bins.delta.delta.y <- (max.delta.delta.y / delta.optimal) %>% ceiling
  bin.breaks.y <- seq(-.5,n_bins.delta.delta.y + 1, 1)

  sdar.hist.y <- graphics::hist(normalized.resolution$delta.delta.y / delta.optimal,
                                breaks = bin.breaks.y, plot = showPlots,
                                main = "Histogram for assessment of y-Resolution",
                                xlab = "\u2206\u2206yi (binned by \u03B4)", warn.unused = F)

  z.delta.delta.y <- sdar.hist.y[["counts"]] %>%
    utils::tail(-1) %>%
    which.max
  relativeResolution.y <- z.delta.delta.y/3
  passed.check.resolution.y <- ifelse(relativeResolution.y <= 1, T, F)
  passed.check.binz.ratio.y <- ifelse((sdar.hist.y[["counts"]][z.delta.delta.y+1])/samples > 0.25, F, T)
  passed.check.bin0.ratio.y <- ifelse((sdar.hist.y[["counts"]][1])/samples > 0.25, F, T)

  # secondary checks. for resolution if first have failed
  if(!passed.check.resolution.x) {
    passed.check.resolution.x <- passed.check.binz.ratio.x && passed.check.binz.ratio.x
  }
  if(!passed.check.resolution.y) {
    passed.check.resolution.y <- passed.check.binz.ratio.y && passed.check.binz.ratio.y
  }

  # print messages
  if(verbose) {
    message("  Data Quality Metric: Digital Resolution")
    message(paste0("      Relative x-Resolution:   ", relativeResolution.x))
    message(paste0("      Relative y-Resolution:   ", relativeResolution.y, "\n"))
  }

  # Warn in case of insufficient resolution
  if(warn) {
    if(!passed.check.resolution.x) {
      warning(paste0("relative digital resolution of x-values might be insufficient (",
                     round(relativeResolution.x*100, 1), " %)\n"), call. = F)
    }

    if(!passed.check.resolution.y) {
      warning(paste0("relative digital resolution of y-values might be insufficient (",
                     round(relativeResolution.y*100, 1), " %)\n"), call. = F)
    }
  }

  # passed all checks for resolution?
  passed.check.resolution <- passed.check.resolution.x && passed.check.resolution.y

  # assemble quality metrics
  results <- list(
    "passed.check" = passed.check.resolution,
    "x" = list(
      "Relative_resolution" = relativeResolution.x,
      "Percent_at_this_resolution" = sdar.hist.x[["counts"]][z.delta.delta.x+1]/samples*100,
      "Percent_in_zeroth_bin" = sdar.hist.x[["counts"]][1]/samples*100,
      "passed.check" = passed.check.resolution.x
    ),
    "y" = list(
      "Relative_resolution" = relativeResolution.y,
      "Percent_at_this_resolution" = sdar.hist.y[["counts"]][z.delta.delta.y+1]/samples*100,
      "Percent_in_zeroth_bin" = sdar.hist.y[["counts"]][1]/samples*100,
      "passed.check" = passed.check.resolution.y
    )
  )

  # add plots if requested
  if(savePlots) {

    hist.x <- carrier::crate(function(value) {
      graphics::hist(x,
                     breaks,
                     main = plot.main,
                     xlab = plot.xlab)
    }, breaks = bin.breaks.x,
    x = normalized.resolution$delta.delta.x / delta.optimal,
    plot.main = "Histogram for assessment of x-Resolution",
    plot.xlab = "\u2206\u2206xi (binned by \u03B4)")

    hist.y <- carrier::crate(function(value) {
      graphics::hist(x,
                     breaks,
                     main = plot.main,
                     xlab = plot.xlab)
    }, breaks = bin.breaks.y,
    x = normalized.resolution$delta.delta.y / delta.optimal,
    plot.main = "Histogram for assessment of y-Resolution",
    plot.xlab = "\u2206\u2206yi (binned by \u03B4)")

    results <- results %>% append(list(
      "plots" = list(
        "hist.x" = hist.x,
        "hist.y" = hist.y
      )
    ))
  }

  # return results
  results
}


#' Data Quality Metrics
#'
#' @note
#' Refer to chapters 6.4 & 6.5 in ASTM E3076-18
#'
#' @description
#' Check for noise and resolution of the normalized data
#'
#' @noRd
check_data_quality <- function(data.normalized,
                               verbose = F,
                               showPlots = F,
                               savePlots = F) {

  # get and check data quality metrics
  list(
    "Num.Obs.normalized" = nrow(data.normalized),
    "digitalResolution" = check_data_quality.resolution(data.normalized,
                                                        verbose = verbose,
                                                        showPlots = showPlots,
                                                        savePlots = savePlots),
    "Noise" = check_data_quality.noise(data.normalized, verbose)) %>% {
      results <- .
      list(
        "passed.check" = results[["digitalResolution"]][["passed.check"]] && results[["Noise"]][["passed.check"]]) %>% append(results)
    }
}



# simplified version of check_data_quality, will return only the (boolean) result of the checks
check_data_quality.lazy <- function(data.normalized) {
  # satisfy pipe addiction...
  `%>%` <- magrittr::`%>%`

  # get and check both data quality metrics
  digitalResolution <- check_data_quality.resolution(data.normalized, warn = F)
  Noise <- check_data_quality.noise(data.normalized, warn = F)

  list("passed.check" = digitalResolution$passed.check && Noise$passed.check,
       "passed.check.resolution.x" = digitalResolution$x$passed.check,
       "passed.check.resolution.y" = digitalResolution$y$passed.check,
       "passed.check.noise.x" = Noise$x$passed.check,
       "passed.check.noise.y" = Noise$y$passed.check)

}


#' Fit Quality Metrics
#'
#' @note
#' Refer to chapters 6.7 & 6.8 in ASTM E3076-18
#'
#' @description
#' Check for curvature in the vicinity of the fit range
#'
#' @noRd
check_fit_quality <- function(data.normalized,
                              fit,
                              verbose = F,
                              showPlots = F,
                              savePlots = F,
                              warn = T) {

  # shorthands for slope and intercept
  m <- fit[["m"]]
  b <- fit[["b"]]

  # filter data for analysis
  data.filtered <- data.normalized %>%
    dplyr::filter(dplyr::between(.data[["otr.idx"]],
                                 fit[["otr.idx.start"]],
                                 fit[["otr.idx.end"]])) %>%
    dplyr::mutate("y.residual" = .data[["y.normalized"]] - (m * .data[["x.normalized"]] + b))

  # find range of data in normalized x and y
  x.min <- data.filtered$x.normalized %>% min
  x.max <- data.filtered$x.normalized %>% max
  x.range <- x.max - x.min
  y.min <- data.filtered$y.normalized %>% min
  y.max <- data.filtered$y.normalized %>% max
  y.range <- y.max - y.min

  # check slope of residuals in outer quartiles
  allowable_residual_slope <- m * 0.05

  data.first_quartile <- data.filtered %>%
    dplyr::filter(x.normalized <= (x.min + 0.25 * x.range))

  slope.first_quartile <- stats::.lm.fit(cbind(1, data.first_quartile$x.normalized),
                                         data.first_quartile$y.residual)$coefficients[2]

  data.fourth_quartile <- data.filtered %>%
    dplyr::filter(x.normalized >= (x.min + 0.75 * x.range))

  slope.fourth_quartile <- stats::.lm.fit(cbind(1, data.fourth_quartile$x.normalized),
                                          data.fourth_quartile$y.residual)$coefficients[2]

  # Assemble first fit quality metrics
  number_of_points.first_quartile <- nrow(data.first_quartile)
  residual_slope.first_quartile <- slope.first_quartile
  relative_residual_slope.first_quartile <- slope.first_quartile/allowable_residual_slope
  check.passed.first_quartile <- ifelse(number_of_points.first_quartile >= 5,
                                        ifelse(abs(relative_residual_slope.first_quartile) <= 1, T, F),
                                        F)

  number_of_points.fourth_quartile <- nrow(data.fourth_quartile)
  residual_slope.fourth_quartile <- slope.fourth_quartile
  relative_residual_slope.fourth_quartile <- slope.fourth_quartile/allowable_residual_slope
  check.passed.fourth_quartile <- ifelse(number_of_points.fourth_quartile >= 5,
                                         ifelse(abs(relative_residual_slope.fourth_quartile) <= 1, T, F),
                                         F)

  # Assemble second fit quality metrics
  normalized_y_range <- y.range
  relative_fit_range <- 0.4/y.range
  check.passed.relative_fit_range <- ifelse(relative_fit_range > 1, F, T)

  # print messages
  if(verbose) {
    message("  Fit Quality Metric: Curvature")
    message("    1st Quartile")
    message(paste0("      Number of Points:         ", number_of_points.first_quartile))
    message(paste0("      Residual Slope:           ", residual_slope.first_quartile))
    message(paste0("      Allowable Residual Slope: ", allowable_residual_slope))
    message(paste0("      Relative Residual Slope:  ", relative_residual_slope.first_quartile))
    message("    4th Quartile")
    message(paste0("      Number of Points:         ", number_of_points.fourth_quartile))
    message(paste0("      Residual Slope:           ", residual_slope.fourth_quartile))
    message(paste0("      Allowable Residual Slope: ", allowable_residual_slope))
    message(paste0("      Relative Residual Slope:  ", relative_residual_slope.fourth_quartile, "\n"))

    message("  Second Fit Quality Metric")
    message(paste0("      normalized y-range:       ", y.min, " - ", y.max))
    message(paste0("      relative fit range:       ", relative_fit_range, "\n"))
  }

  # give warning when quality metric checks have failed
  if(warn) {
    if(!(check.passed.first_quartile && check.passed.fourth_quartile)) {
      warning("First Fit Quality Metric checks failed: Excessive curvature found in the vicinity of the fitted range!\n", call. = F)
    }

    if(!check.passed.relative_fit_range) {
      warning("Second Fit Quality Metric check failed: unacceptably small linear region!\n", call. = F)
    }
  }

  # Let's get plotty...
  if(showPlots || savePlots) {
    # generate function for plotting the residuals over normalized x
    plot.residuals <- carrier::crate(function(value) {
      # satisfy pipe addiction...
      `%>%` <- magrittr::`%>%`

      # shorthands
      y.residual.min <- unlist(plot.data$y.residual, T, F) %>% min
      y.residual.max <- unlist(plot.data$y.residual, T, F) %>% max
      y.residual.range <- y.residual.max - y.residual.min
      x.min <- unlist(plot.data$x.normalized, T, F) %>% min
      x.max <- unlist(plot.data$x.normalized, T, F) %>% max
      x.range <- x.max - x.min

      # generate plot
      plot(x = unlist(plot.data$x.normalized, T, F),
           y = unlist(plot.data$y.residual, T, F),
           type ="p",
           main = plot.main,
           xlab = plot.xlab,
           ylab = plot.ylab,
           xlim = c(x.min - 0.05 * x.range, x.max + 0.05* x.range),
           ylim = c(y.residual.min - 0.05 * y.residual.range, y.residual.max + 0.1 * y.residual.range))

      # add vertical lines, arrows, and text for the 1st and 4th quartiles
      graphics::abline(v = c(x.end.first_quartile, x.start.fourth_quartile), lty= c("dotdash", "dotdash"))
      graphics::arrows(x0 = x.end.first_quartile,
                       y0 = y.residual.max + 0.05 * y.residual.range,
                       x1 = x.min,
                       y1 = y.residual.max + 0.05 * y.residual.range,
                       length = 0.1)
      graphics::text(x.min + 0.125 * x.range, y.residual.max + 0.1 * y.residual.range, "1st Quartile")
      graphics::arrows(x0 = x.start.fourth_quartile,
                       y0 = y.residual.max + 0.05 * y.residual.range,
                       x1 = x.max,
                       y1 = y.residual.max + 0.05 * y.residual.range,
                       length = 0.1)
      graphics::text(x.max - 0.125 * x.range, y.residual.max + 0.1 * y.residual.range, "4th Quartile")
    },
    plot.data = data.frame("x.normalized" = data.filtered$x.normalized %>% unlist(T, F),
                           "y.residual" = data.filtered$y.residual %>% unlist(T, F)),
    plot.x = "x.normalized",
    plot.y = "y.residual",
    plot.xlab = "normalized x",
    plot.ylab = "normalized y residual",
    plot.main = "Residuals for Final Fit Test Record",
    x.end.first_quartile = x.min + 0.25 * x.range,
    x.start.fourth_quartile = x.min + 0.75 * x.range
    )

    # Plot Residuals
    if(showPlots) {
      plot.residuals()
    }
  }

  # assemble quality metrics. results
  results <- list(
    "passed.check" = check.passed.first_quartile && check.passed.fourth_quartile && check.passed.relative_fit_range,
    "Curvature" = list(
      "first_quartile" = list(
        "relative_residual_slope" = relative_residual_slope.first_quartile,
        "Number_of_points" = number_of_points.first_quartile,
        "passed.check" = check.passed.first_quartile
      ),

      "fourth_quartile" = list(
        "relative_residual_slope" = relative_residual_slope.fourth_quartile,
        "Number_of_points" = number_of_points.fourth_quartile,
        "passed.check" = check.passed.fourth_quartile
      )
    ),
    "Fit_range" = list(
      "relative_fit_range" = relative_fit_range,
      "passed.check" = check.passed.relative_fit_range
    )
  )

  # append plot
  if(savePlots) {
    results <- results %>% append(list("plots" = list(
      "plot.residuals" = plot.residuals
    )))
  }

  # return results
  return(results)
}

#' Fit Quality Metrics
#'
#' @note
#' lazy variant
#'
#' @description
#' check fit quality and only return only the (boolean) result
#'
#' @noRd
check_fit_quality.lazy <- function(data.normalized, fit) {
  fit_quality_metrics <- check_fit_quality(data.normalized, fit, F, F, F, F)

  list("passed.check" = fit_quality_metrics$passed.check,
       "fit.quality.passed.first_quartile" = fit_quality_metrics$Curvature$first_quartile$passed.check,
       "fit.quality.passed.fourth_quartile" = fit_quality_metrics$Curvature$fourth_quartile$passed.check,
       "fit.quality.passed.Fit_range" = fit_quality_metrics$Fit_range$passed.check
  )
}
