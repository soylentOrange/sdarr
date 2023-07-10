# # Run the SDAR algorithm with down-sampling
# # @noRd
# sdar_execute.downsample <- function(prepared_data,
#                                     otr.info,
#                                     normalized_data.hint = NULL,
#                                     denoise.x = FALSE,
#                                     denoise.y = FALSE,
#                                     vmd.alpha = 1000,
#                                     fit.candidates,
#                                     optimum.range.size,
#                                     cutoff_probability,
#                                     quality_penalty = 0.1,
#                                     verbose.all = FALSE,
#                                     verbose = TRUE,
#                                     plot.all = FALSE,
#                                     plot = TRUE,
#                                     plotFun.all = FALSE,
#                                     plotFun = FALSE,
#                                     Nmin_factor = 0.2) {
#   # Give a (long) welcome message
#   if (verbose) {
#     message("Down-sampling modification of the SDAR-algorithm\n")
#   }
#
#   # maybe the offset for step 1 needs to be raised later
#   raise_offset_times <- 0
#   optimum_fits_are_found <- FALSE
#
#   # do step 1, check data quality, and get optimum fit
#   # repeat until the upper index is not the last point in the optimum region
#   while (!optimum_fits_are_found) {
#     # step 1, normalize
#     if (!is.null(normalized_data.hint)) {
#       # use given hint
#       normalized_data <- normalized_data.hint
#       normalized_data.hint <- NULL
#     } else {
#       # do the normalization
#       normalized_data <- normalize_data(
#         data = prepared_data,
#         otr.info = otr.info,
#         raise_offset_times = raise_offset_times,
#         denoise.x = denoise.x,
#         denoise.y = denoise.y,
#         vmd.alpha = vmd.alpha,
#         verbose = verbose.all,
#         plot = plot.all,
#         plotFun = plotFun.all
#       )
#     }
#
#     # find optimum fits and determine summary of data & fit quality metrics
#     optimum_fits <- find_linear_regressions.subsampled(
#       normalized_data$data.normalized,
#       1,
#       fit.candidates,
#       optimum.range.size
#     ) %>%
#       # penalize by data quality metrics
#       dplyr::mutate(data.quality.penalty = 1.0) %>%
#       dplyr::mutate(data.quality.penalty = dplyr::case_when(
#         data.quality.check.passed.check.noise.y == FALSE ~ data.quality.penalty * quality_penalty,
#         TRUE ~ data.quality.penalty
#       )) %>%
#       dplyr::mutate(data.quality.penalty = dplyr::case_when(
#         data.quality.check.passed.check.noise.x == FALSE ~ data.quality.penalty * quality_penalty,
#         TRUE ~ data.quality.penalty
#       )) %>%
#       dplyr::mutate(data.quality.penalty = dplyr::case_when(
#         data.quality.check.passed.check.resolution.y == FALSE ~ data.quality.penalty * quality_penalty,
#         TRUE ~ data.quality.penalty
#       )) %>%
#       dplyr::mutate(data.quality.penalty = dplyr::case_when(
#         data.quality.check.passed.check.resolution.x == FALSE ~ data.quality.penalty * quality_penalty,
#         TRUE ~ data.quality.penalty
#       )) %>%
#       dplyr::mutate(offset_raise.weighed = dplyr::case_when(
#         offset_raise_required ~ -data.quality.penalty,
#         TRUE ~ data.quality.penalty
#       ))
#
#     # judge by weighed majority decision if an offset raise is required
#     if (mean(optimum_fits$offset_raise.weighed) >= 0) {
#       optimum_fits_are_found <- TRUE
#     }
#
#     # repeat otherwise
#     if (!optimum_fits_are_found) {
#       raise_offset_times <- raise_offset_times + 1
#     }
#
#     # next offset is out of data range
#     if (raise_offset_times > 17) {
#       stop("Data is unfit for processing: upper index of the optimum fit region is beyond the last point in the truncated test record.")
#     }
#   }
#
#   # discard fits when offset raise would be required or numerical problem in fitting occurred
#   optimum_fits.nrow <- nrow(optimum_fits)
#   optimum_fits <- optimum_fits %>%
#     dplyr::filter(.data$offset_raise_required == FALSE) %>%
#     dplyr::mutate("norm.residual" = dplyr::case_when(
#       .data$norm.residual < 10 * .Machine$double.eps ~ 10 * .Machine$double.eps,
#       TRUE ~ .data$norm.residual
#     ))
#
#   # calculate success rate of quality checks
#   passed.check.data <- nrow(optimum_fits %>%
#     dplyr::filter(.data$data.quality.check.passed.check == TRUE)) / optimum_fits.nrow
#   passed.check.fit <- nrow(optimum_fits %>%
#     dplyr::filter(.data$passed.check == TRUE)) / optimum_fits.nrow
#   passed.check <- nrow(optimum_fits %>%
#     dplyr::filter(.data$data.quality.check.passed.check == TRUE) %>%
#     dplyr::filter(.data$passed.check == TRUE)) / optimum_fits.nrow
#
#   # give a message on success rate of data and fit quality checks
#   if (verbose) {
#     if (denoise.x || denoise.y) {
#       message("  Random sub-sampling information (after de-noising was applied):")
#     } else {
#       message("  Random sub-sampling information:")
#     }
#
#     message(paste0(
#       "      ", optimum.range.size,
#       " points of ", nrow(normalized_data$data.normalized),
#       " points in the normalized range were used."
#     ))
#     message(paste0("      ", round(passed.check.data * 100, 1), " % of sub-sampled normalized ranges passed the data quality checks."))
#     message(paste0("      ", round(passed.check.fit * 100, 1), " % of linear regressions passed the fit quality checks."))
#     message(paste0("      ", round(passed.check * 100, 1), " % of linear regressions passed all quality checks.\n  "))
#   }
#
#   # get numerically stable fits and add weights
#   optimum_fits <- optimum_fits %>%
#     # penalize by fit quality metrics
#     dplyr::mutate("fit.quality.penalty" = 1.0) %>%
#     dplyr::mutate("fit.quality.penalty" = dplyr::case_when(
#       fit.quality.passed.Fit_range == FALSE ~ .data$fit.quality.penalty * quality_penalty,
#       TRUE ~ .data$fit.quality.penalty
#     )) %>%
#     dplyr::mutate("fit.quality.penalty" = dplyr::case_when(
#       .data$fit.quality.passed.fourth_quartile == FALSE ~ .data$fit.quality.penalty * quality_penalty,
#       TRUE ~ .data$fit.quality.penalty
#     )) %>%
#     dplyr::mutate("fit.quality.penalty" = dplyr::case_when(
#       .data$fit.quality.passed.first_quartile == FALSE ~ .data$fit.quality.penalty * quality_penalty,
#       TRUE ~ .data$fit.quality.penalty
#     )) %>%
#     dplyr::mutate("quality.penalty" = .data$fit.quality.penalty * .data$data.quality.penalty) %>%
#     dplyr::mutate("weight" = .data$quality.penalty / .data$norm.residual) %>%
#     dplyr::select(-c(
#       "offset_raise_required",
#       "data.quality.penalty",
#       "fit.quality.penalty",
#       "quality.penalty",
#       "offset_raise.weighed"
#     ))
#
#   # re-do normalization without de-noising
#   plot.otr.denoised <- NULL
#   if (denoise.x || denoise.y) {
#     if(savePlots) {
#       # save the plot of de-noised data
#       plot.otr.denoised <- normalized_data$plots$plot.otr.denoised
#     }
#
#     normalized_data <- normalize_data(
#       data = prepared_data,
#       otr.info = otr.info,
#       raise_offset_times = raise_offset_times,
#       denoise.x = FALSE,
#       denoise.y = FALSE,
#       verbose = verbose.all,
#       plot = plot.all,
#       plotFun = plotFun.all
#     )
#   }
#
#   # shorthands for calculations
#   y.tangent <- normalized_data[["tangency.point"]][["y.tangent"]]
#   x.tangent <- normalized_data[["tangency.point"]][["x.tangent"]]
#   x.shift <- normalized_data[["shift"]][["x.shift"]]
#   y.shift <- normalized_data[["shift"]][["y.shift"]]
#
#   # Un-normalize fits
#   optimum_fits <- optimum_fits %>%
#     dplyr::mutate(y.min = normalized_data$data.normalized[.data$otr.idx.start, "y.normalized"]) %>%
#     dplyr::mutate(y.max = normalized_data$data.normalized[.data$otr.idx.end, "y.normalized"]) %>%
#     dplyr::mutate(x.min = normalized_data$data.normalized[.data$otr.idx.start, "x.normalized"]) %>%
#     dplyr::mutate(x.max = normalized_data$data.normalized[.data$otr.idx.end, "x.normalized"]) %>%
#     dplyr::mutate(y.lowerBound = .data$y.min * y.tangent + y.shift) %>%
#     dplyr::mutate(y.upperBound = .data$y.max * y.tangent + y.shift) %>%
#     dplyr::mutate(x.lowerBound = .data$x.min * x.tangent + x.shift) %>%
#     dplyr::mutate(x.upperBound = .data$x.max * x.tangent + x.shift)
#
#   # use linear regression to find weighed means of results
#   # (otr.idx.start and .end, y and x bounds)
#   # will be used for reporting the results and one final fit
#   weighed.results <- optimum_fits %>%
#     dplyr::select(c(
#       "weight", "otr.idx.start", "otr.idx.end",
#       "y.lowerBound", "y.upperBound",
#       "x.lowerBound", "x.upperBound"
#     )) %>%
#     as.list() %>%
#     {
#       res <- .
#       weight <- .[["weight"]]
#       res$weight <- NULL
#       res %>%
#         purrr::map(carrier::crate(function(value) {
#           list(
#             "value" = value,
#             "weight" = weight
#           )
#         }, weight = weight)) %>%
#         purrr::map(carrier::crate(function(value) {
#           fit <- stats::lm(value ~ 1,
#             data = data.frame(
#               "value" = value$value,
#               "weights" = value$weight
#             ),
#             weights = weights
#           )
#           conf <- stats::confint(fit)
#           list(
#             "value" = fit$coefficients[[1]],
#             "conf.low" = conf[[1]],
#             "conf.high" = conf[[2]]
#           )
#         }))
#     } %>%
#     as.data.frame() %>%
#     {
#       value <- .
#       new.names <- names(value) %>% stringr::str_replace(stringr::coll(".value"), "")
#       value %>% magrittr::set_names(new.names)
#     }
#
#   # check data quality on the finally selected range
#   data_quality_metrics <- check_data_quality(
#     data.normalized = normalized_data$data.normalized,
#     verbose = verbose.all,
#     plot = plot.all,
#     plotFun = plotFun.all
#   )
#
#   # linear regression on (complete) normalized data
#   optimum.otr.idx.start <- round(weighed.results[[1, "otr.idx.start"]], 0)
#   optimum.otr.idx.end <- round(weighed.results[[1, "otr.idx.end"]], 0)
#   optimum.fit <- stats::lm(y.normalized ~ x.normalized,
#     data = normalized_data$data.normalized %>%
#       dplyr::filter(dplyr::between(
#         .data$otr.idx,
#         optimum.otr.idx.start,
#         optimum.otr.idx.end
#       ))
#   )
#
#   optimum_fit <- data.frame(
#     "m" = optimum.fit$coefficients[[2]],
#     "b" = optimum.fit$coefficients[[1]],
#     "otr.idx.start" = optimum.otr.idx.start,
#     "otr.idx.end" = optimum.otr.idx.end
#   )
#
#   # check fit quality
#   fit_quality_metrics <- check_fit_quality(
#     normalized_data$data.normalized,
#     optimum_fit,
#     verbose.all,
#     showPlots.all,
#     savePlots
#   )
#
#   # un-normalize data and summarize
#   assembled_results <- assemble_report(
#     normalized_data = normalized_data,
#     otr.info = otr.info,
#     dataQualityMetrics = data_quality_metrics,
#     fit = optimum_fit,
#     fitQualityMetrics = fit_quality_metrics,
#     verbose = verbose,
#     plot = plot,
#     plotFun.all = plotFun.all,
#     plotFun = plotFun
#   )
#
#   # re-do the final fit to find confidence intervals for slope and intercepts
#   final.fit.confint <- prepared_data[optimum_fit$otr.idx.start:optimum_fit$otr.idx.end, ] %>%
#     {
#       prepared_data.fit.lm <- stats::lm(y.data ~ x.data, data = .)
#     } %>%
#     stats::confint()
#
#   results.labels <- list(
#     "finalSlope" = assembled_results$Slope_Determination_Results$finalSlope.unit,
#     "finalSlope.conf.low" = assembled_results$Slope_Determination_Results$finalSlope.unit,
#     "finalSlope.conf.high" = assembled_results$Slope_Determination_Results$finalSlope.unit,
#     "trueIntercept" = assembled_results$Slope_Determination_Results$trueIntercept.unit,
#     "trueIntercept.conf.low" = assembled_results$Slope_Determination_Results$trueIntercept.unit,
#     "trueIntercept.conf.high" = assembled_results$Slope_Determination_Results$trueIntercept.unit,
#     "y.lowerBound" = assembled_results$Slope_Determination_Results$y.unit,
#     "y.upperBound" = assembled_results$Slope_Determination_Results$y.unit,
#     "x.lowerBound" = assembled_results$Slope_Determination_Results$x.unit,
#     "x.upperBound" = assembled_results$Slope_Determination_Results$x.unit,
#     "y.lowerBound.conf.low" = assembled_results$Slope_Determination_Results$y.unit,
#     "y.upperBound.conf.low" = assembled_results$Slope_Determination_Results$y.unit,
#     "x.lowerBound.conf.low" = assembled_results$Slope_Determination_Results$x.unit,
#     "x.upperBound.conf.low" = assembled_results$Slope_Determination_Results$x.unit,
#     "y.lowerBound.conf.high" = assembled_results$Slope_Determination_Results$y.unit,
#     "y.upperBound.conf.high" = assembled_results$Slope_Determination_Results$y.unit,
#     "x.lowerBound.conf.high" = assembled_results$Slope_Determination_Results$x.unit,
#     "x.upperBound.conf.high" = assembled_results$Slope_Determination_Results$x.unit,
#     "Num.Obs.fit" = "(points in final fit)",
#     "Num.Obs.downsampled" = "(points in down-sampled data range)",
#     "Num.Obs.normalized" = "(points in normalized data range)",
#     "passed.check" = "(all checks)",
#     "passed.check.data" = "(data quality checks)",
#     "passed.check.fit" = "(fit quality checks)"
#   )
#
#   sdar.results <- data.frame(
#     "finalSlope" = assembled_results$Slope_Determination_Results$finalSlope,
#     "finalSlope.conf.low" = final.fit.confint["x.data", 1],
#     "finalSlope.conf.high" = final.fit.confint["x.data", 2],
#     "trueIntercept" = assembled_results$Slope_Determination_Results$trueIntercept,
#     "trueIntercept.conf.low" = final.fit.confint["(Intercept)", 1],
#     "trueIntercept.conf.high" = final.fit.confint["(Intercept)", 2],
#     "y.lowerBound" = assembled_results$Slope_Determination_Results$y.lowerBound,
#     "y.lowerBound.conf.low" = weighed.results[[1, "y.lowerBound.conf.low"]],
#     "y.lowerBound.conf.high" = weighed.results[[1, "y.lowerBound.conf.high"]],
#     "y.upperBound" = assembled_results$Slope_Determination_Results$y.upperBound,
#     "y.upperBound.conf.low" = weighed.results[[1, "y.upperBound.conf.low"]],
#     "y.upperBound.conf.high" = weighed.results[[1, "y.upperBound.conf.high"]],
#     "x.lowerBound" = assembled_results$Slope_Determination_Results$x.lowerBound,
#     "x.lowerBound.conf.low" = weighed.results[[1, "x.lowerBound.conf.low"]],
#     "x.lowerBound.conf.high" = weighed.results[[1, "x.lowerBound.conf.high"]],
#     "x.upperBound" = assembled_results$Slope_Determination_Results$x.upperBound,
#     "x.upperBound.conf.low" = weighed.results[[1, "x.upperBound.conf.low"]],
#     "x.upperBound.conf.high" = weighed.results[[1, "x.upperBound.conf.high"]],
#     "Num.Obs.fit" = nrow(prepared_data[optimum_fit$otr.idx.start:optimum_fit$otr.idx.end, ]),
#     "Num.Obs.downsampled" = optimum.range.size,
#     "Num.Obs.normalized" = assembled_results$Data_Quality_Metrics$Num.Obs.normalized,
#     "passed.check" = assembled_results$Slope_Determination_Results$passed.check,
#     "passed.check.data" = assembled_results$Data_Quality_Metrics$passed.check,
#     "passed.check.fit" = assembled_results$Fit_Quality_Metrics$passed.check,
#     "method" = "SDAR with down-sampling"
#   ) %>%
#     labelled::set_variable_labels(.labels = results.labels)
#
#   # return full results
#   results <- list(
#     "sdar" = sdar.results,
#     "Data_Quality_Metrics" = assembled_results$Data_Quality_Metrics,
#     "Fit_Quality_Metrics" = assembled_results$Fit_Quality_Metrics
#   )
#
#   # append plot
#   if (plotFun || plotFun.all) {
#     results <- results %>% append(list("plots" = list(
#       "final.fit" = assembled_results$plots$plot.fit
#     )))
#   }
#
#   # append more plots
#   if (plotFun.all) {
#     results$plots <- results$plots %>% append(list(
#       "otr" = assembled_results$plots$plot.otr,
#       "normalized" = assembled_results$plots$plot.normalized,
#       "hist.x" = assembled_results$Data_Quality_Metrics$digitalResolution$plots$hist.x,
#       "hist.y" = assembled_results$Data_Quality_Metrics$digitalResolution$plots$hist.y,
#       "residuals" = assembled_results$Fit_Quality_Metrics$plots$plot.residuals
#     ))
#
#     # append plot of de-noised (partial) data, if available
#     if(!is.null(plot.otr.denoised)) {
#       results$plots <- results$plots %>% append(list(
#         "otr.denoised" = plot.otr.denoised
#       ))
#     }
#
#     results$Data_Quality_Metrics$digitalResolution$plots <- NULL
#     results$Fit_Quality_Metrics$plots <- NULL
#   }
#
#   # return results
#   return(results)
# }


#' @title Down-sampling variant of the SDAR-algorithm
#'
#' @description Run a down-sampling modification of the SDAR algorithm as
#'   originally standardized in ASTM E3076-18. As the original version uses
#'   numerous linear regressions (`.lm.fit()` from the stats-package), which can
#'   be painfully slow for test data with high resolution. The lazy variant of
#'   the algorithm will use several random sub-samples of the data to find the
#'   best estimate for the fit-range within the data. Additionally, the test
#'   data will be de-noised using Variational Mode Decomposition in case initial
#'   data quality checks have failed. See the articles on [validation](
#'   https://soylentorange.github.io/sdarr/articles/sdarr_validation.html) and
#'   [robustness against
#'   noise](https://soylentorange.github.io/sdarr/articles/excessive_noise_levels.html)
#'   on the [package-website](https://soylentorange.github.io/sdarr/) for
#'   further information.
#'
#' @note The function can use parallel processing via the
#'   [furrr-package](https://furrr.futureverse.org/). To use this feature, set
#'   up a plan other than the default sequential strategy beforehand. Also, as
#'   random values are drawn, set a [`random seed`][base::RNG] beforehand to get
#'   reproducible results.
#'
#' @references Lucon, E. (2019). Use and validation of the slope determination
#'   by the analysis of residuals (SDAR) algorithm (NIST TN 2050; p. NIST TN
#'   2050). National Institute of Standards and Technology.
#'   https://doi.org/10.6028/NIST.TN.2050
#'
#' @references Standard Practice for Determination of the Slope in the Linear
#'   Region of a Test Record (ASTM E3076-18). (2018).
#'   https://doi.org/10.1520/E3076-18
#'
#' @references Graham, S., & Adler, M. (2011). Determining the Slope and Quality
#'   of Fit for the Linear Part of a Test Record. Journal of Testing and
#'   Evaluation - J TEST EVAL, 39. https://doi.org/10.1520/JTE103038
#'
#' @references Dragomiretskiy, K., & Zosso, D. (2014). Variational Mode
#'   Decomposition. IEEE Transactions on Signal Processing, 62(3), 531–544.
#'   https://doi.org/10.1109/TSP.2013.2288675
#'
#' @inheritParams sdar
#'
#' @param cutoff_probability Cut-off probability for estimating optimum size of
#'   down-sampled data range via logistic regression.
#'
#' @seealso [sdar()] for the standard SDAR-algorithm.
#'
#' @returns A list containing a data.frame with the results of the final fit,
#'   lists with the quality- and fit-metrics, and a list containing the crated
#'   plot-function(s) (if `plotFun = TRUE`).
#'
#' @examples
#' # Synthesize a test record resembling Al 6060 T66
#' # (Values according to Metallic Material Properties
#' # Development and Standardization (MMPDS) Handbook).
#' # Explicitly set names to "strain" and "stress".
#'
#' Al_6060_T66 <- synthesize_test_data(
#'   slope = 68000,
#'   yield.y = 160,
#'   ultimate.y = 215,
#'   ultimate.x = 0.091,
#'   x.name = "strain",
#'   y.name = "stress"
#' )
#'
#'
#' # use sdar.downsample() to analyze the synthetic test record
#' # (using relaxed settings for the noise-free synthetic data)
#' # will print a report and give a plot of the final fit
#' \donttest{
#' result <- sdar.downsample(Al_6060_T66, strain, stress)
#' }
#'
#' @noRd

# @export
sdar.downsample <- function(data, x, y,
                            verbose = TRUE,
                            plot = TRUE,
                            plotFun = FALSE,
                            cutoff_probability = 0.875,
                            ...) {

  # to be furrr-safe, enquote the tidy arguments here
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)

  # take care of dynamic dots
  additional_parameters <- rlang::list2(...)
  Nmin_factor <- try({additional_parameters$Nmin_factor}, silent = TRUE)
  if(!is.numeric(Nmin_factor)) Nmin_factor <- 0.2
  verbose.all <- try({additional_parameters$verbose.all}, silent = TRUE)
  if(!is.logical(verbose.all)) verbose.all <- FALSE
  plot.all <- try({additional_parameters$plot.all}, silent = TRUE)
  if(!is.logical(plot.all)) plot.all <- FALSE
  plotFun.all <- try({additional_parameters$plotFun.all}, silent = TRUE)
  if(!is.logical(plotFun.all)) plotFun.all <- FALSE
  n.candidates <- try({additional_parameters$n.candidates}, silent = TRUE)
  if(!is.numeric(n.candidates)) n.candidates <- 20
  enforce_denoising <- try({additional_parameters$enforce_denoising}, silent = TRUE)
  if(!is.logical(enforce_denoising)) enforce_denoising <- FALSE
  quality_penalty <- try({additional_parameters$quality_penalty}, silent = TRUE)
  if(!is.numeric(quality_penalty)) quality_penalty <- 0.1
  vmd.alpha <- try({additional_parameters$vmd.alpha}, silent = TRUE)
  if(!is.numeric(vmd.alpha)) vmd.alpha <- 1000
  min_K <- try({additional_parameters$min_K}, silent = TRUE)
  if(!is.numeric(min_K)) min_K <- 3

  # save final fit plot, when plotFun.all is set
  plotFun <- plotFun || plotFun.all

  # give messages for report, when verbose.all is set
  verbose <- verbose || verbose.all

  # show final fit plot, when plot.all is set
  plot <- plot || plot.all

  # get units for data
  x.label.unit <- data %>%
    dplyr::select(!!x) %>%
    labelled::var_label() %>%
    {
      .[[1]]
    }
  y.label.unit <- data %>%
    dplyr::select(!!y) %>%
    labelled::var_label() %>%
    {
      .[[1]]
    }

  # get names of data
  x.name <- data %>%
    dplyr::select(!!x) %>%
    names() %>%
    {
      .[[1]]
    }
  y.name <- data %>%
    dplyr::select(!!y) %>%
    names() %>%
    {
      .[[1]]
    }

  # prepare data, add an index and remove NA from data
  prepared_data <- data.frame(
    x.data = data %>% dplyr::select(!!x) %>% unlist(TRUE, FALSE),
    y.data = data %>% dplyr::select(!!y) %>% unlist(TRUE, FALSE),
    otr.idx = seq.int(nrow(data))
  ) %>%
    tidyr::drop_na()

  # assemble information of the original test record
  otr.info <- list(
    "x" = list(
      "name" = x.name,
      "unit" = x.label.unit
    ),
    "y" = list(
      "name" = y.name,
      "unit" = y.label.unit
    )
  )

  # Give a welcome message
  if (verbose) {
    message("Determination of Slope in the Linear Region of a Test Record:")
  }

  if(enforce_denoising) {
    if (verbose) {
      message("  Proceeding with de-noised data...\n")
      denoise.x <- TRUE
      denoise.y <- TRUE

      # use de-noising, as requested
      normalized_data <- normalize_data(
        data = prepared_data,
        otr.info = otr.info,
        raise_offset_times = 0,
        denoise.x = denoise.x,
        denoise.y = denoise.y,
        vmd.alpha = vmd.alpha,
        min_K = min_K,
        verbose = FALSE,
        plot = FALSE,
        plotFun = TRUE
      )
    }
  } else {
    # normalize data for initial quality check
    normalized_data <- normalize_data(
      data = prepared_data,
      otr.info = otr.info,
      raise_offset_times = 0,
      denoise.x = FALSE,
      denoise.y = FALSE,
      verbose = verbose.all,
      plot = plot.all,
      plotFun = plotFun.all
    )

    # get initial data quality metrics
    data_quality_metrics <- check_data_quality.lazy(normalized_data$data.normalized)
    denoise.x <- !(data_quality_metrics$passed.check.noise.x &&
                     data_quality_metrics$passed.check.resolution.x)
    denoise.y <- !(data_quality_metrics$passed.check.noise.y &&
                     data_quality_metrics$passed.check.resolution.y)

    # Check for noise in data quality metrics
    if (denoise.x || denoise.y) {
      # normalize data for fall-back quality check
      # use de-noising, when quality-metrics indicate a problem
      normalized_data <- normalize_data(
        data = prepared_data,
        otr.info = otr.info,
        raise_offset_times = 0,
        denoise.x = denoise.x,
        denoise.y = denoise.y,
        vmd.alpha = vmd.alpha,
        min_K = min_K,
        verbose = FALSE,
        plot = FALSE,
        plotFun = TRUE
      )

      # get data quality metrics of de-noised data
      data_quality_metrics <- check_data_quality.lazy(normalized_data$data.normalized)

      # check results...
      if (data_quality_metrics$passed.check) {
        if (verbose) {
          message("  Data quality checks of original test record failed!")
          message("  Proceeding with de-noised data...\n")
        }
      }
    }

    # Check data quality metrics and use standard-variant in case of failing to
    # pass check (again)
    if (!data_quality_metrics$passed.check) {
      if (verbose) {
        message("  Data quality checks of de-noised test record failed!")
        message("  Examine plots of Standard SDAR-algorithm to determine how to proceed.\n")
      }

      # execute the standard SDAR-algorithm
      result <- sdar_execute(
        prepared_data = prepared_data,
        otr.info = otr.info,
        verbose.all = verbose.all,
        verbose = verbose,
        plot.all = plot.all,
        plot = plot,
        plotFun.all = plotFun.all,
        plotFun = plotFun,
        Nmin_factor = Nmin_factor
      )
      result$sdar <- result$sdar %>% dplyr::mutate("method" = "SDAR")
      return(result)
    }
  }

  # execute the standard SDAR-algorithm
  result <- sdar_execute(
    prepared_data = prepared_data,
    otr.info = otr.info,
    verbose.all = verbose.all,
    verbose = verbose,
    plot.all = plot.all,
    plot = plot,
    plotFun.all = plotFun.all,
    plotFun = plotFun,
    Nmin_factor = Nmin_factor
  )
  result$sdar <- result$sdar %>% dplyr::mutate("method" = "SDAR")
  return(result)

  # # find optimum range for sub-sampling
  # optimum_size <- optimum_size_for_downsampling(normalized_data$data.normalized,
  #   cutoff_probability = cutoff_probability,
  #   showPlots = showPlots.all,
  #   verbose = verbose.all,
  #   savePlots = savePlots
  # )
  #
  # # check for possible errors in determination of optimum size
  # if (optimum_size$cutoff_probability.matched == FALSE) {
  #   if (enforce_subsampling == TRUE) {
  #     # Give an informational message and proceed using the found quasi-optimum
  #     # size
  #     if (verbose) {
  #       message("  Failed to satisfy cutoff_probability. Lowering our expectations...\n")
  #     }
  #   } else {
  #     if (verbose) {
  #       message("  Failed to satisfy cutoff_probability. Standard SDAR-algorithm will be used...\n")
  #     }
  #
  #     # execute the standard SDAR-algorithm
  #     result <- sdar_execute(
  #       prepared_data = prepared_data,
  #       otr.info = otr.info,
  #       verbose.all = verbose.all,
  #       verbose = verbose,
  #       plot.all = plot.all,
  #       plot = plot,
  #       plotFun.all = plotFun.all,
  #       plotFun = plotFun,
  #       Nmin_factor = Nmin_factor
  #     )
  #     result$sdar <- result$sdar %>% dplyr::mutate("method" = "SDAR")
  #     return(result)
  #   }
  # }
  #
  # # check for viability of down-sampling-approach
  # viability <- optimum_size$optimum.range.size < nrow(normalized_data$data.normalized)
  #
  # # check, whether sub-sampling would save time and execute standard or lazy variant
  # if (viability == FALSE) {
  #   if (verbose) {
  #     message("  Standard SDAR-algorithm will be used...\n")
  #   }
  #
  #   # execute the standard SDAR-algorithm
  #   result <- sdar_execute(
  #     prepared_data = prepared_data,
  #     otr.info = otr.info,
  #     verbose.all = verbose.all,
  #     verbose = verbose,
  #     plot.all = plot.all,
  #     plot = plot,
  #     plotFun.all = plotFun.all,
  #     plotFun = plotFun,
  #     Nmin_factor = Nmin_factor
  #   )
  #   result$sdar <- result$sdar %>% dplyr::mutate("method" = "SDAR")
  #   return(result)
  # }
  #
  # # execute the down-sampling - SDAR-algorithm
  # results <- sdar_execute.downsample(
  #   prepared_data = prepared_data,
  #   otr.info = otr.info,
  #   normalized_data.hint = normalized_data,
  #   denoise.x = denoise.x,
  #   denoise.y = denoise.y,
  #   vmd.alpha = vmd.alpha,
  #   fit.candidates = fit.candidates,
  #   optimum.range.size = optimum_size$optimum.range.size,
  #   cutoff_probability = cutoff_probability,
  #   quality_penalty = quality_penalty,
  #   verbose.all = verbose.all,
  #   verbose = verbose,
  #   showPlots.all = showPlots.all,
  #   showPlots.report = showPlots.report,
  #   savePlots = savePlots
  # )
  #
  # # add a de-noising remark
  # if (denoise.x || denoise.y || enforce_denoising) {
  #   results$sdar$method <- "SDAR with down-sampling on de-noised data"
  # }
  #
  # # append/sort plots to the results
  # if (savePlots) {
  #   results$plots <- results$plots %>% append(list(
  #     "glm.optimum_size" = optimum_size$plots$plot.modelpredictions
  #   ))
  # }
  #
  # # return results
  # return(results)
}
