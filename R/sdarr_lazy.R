#' Run the SDAR algorithm with random sub-sampling
#' @noRd
sdar_execute.lazy <- function(prepared_data,
                              otr.info,
                              normalized_data.hint = NULL,
                              n.fit,
                              n.candidates,
                              optimum.range.size,
                              cutoff_probability,
                              quality_penalty,
                              verbose.all = FALSE,
                              verbose = TRUE,
                              plot.all = FALSE,
                              plot = TRUE,
                              plotFun.all = FALSE,
                              plotFun = FALSE,
                              Nmin_factor = 0.2) {
  # Give a (long) welcome message
  if (verbose) {
    message("Random sub-sampling modification of the SDAR-algorithm\n")
  }

  # maybe the offset for step 1 needs to be raised later
  raise_offset_times <- 0
  optimum_fits_are_found <- FALSE

  # do step 1, check data quality, and get optimum fit
  # repeat until the upper index is not the last point in the optimum region
  while (!optimum_fits_are_found) {
    # step 1, normalize
    if (!is.null(normalized_data.hint)) {
      # use given hint
      normalized_data <- normalized_data.hint
      normalized_data.hint <- NULL
    } else {
      # do the normalization
      normalized_data <- normalize_data(
        data = prepared_data,
        otr.info = otr.info,
        raise_offset_times = raise_offset_times,
        verbose = verbose.all,
        plot = plot.all,
        plotFun = plotFun.all
      )
    }

    # find optimum fits and determine summary of data & fit quality metrics
    optimum_fits <- find_linear_regressions.subsampled(
      normalized_data,
      n.fit,
      n.candidates,
      optimum.range.size,
      Nmin_factor = Nmin_factor
    ) %>%
      # penalize by data quality metrics
      dplyr::mutate(data.quality.penalty = 1.0) %>%
      dplyr::mutate(data.quality.penalty = dplyr::case_when(
        data.quality.check.passed.check.noise.y == FALSE ~ data.quality.penalty * quality_penalty,
        TRUE ~ data.quality.penalty
      )) %>%
      dplyr::mutate(data.quality.penalty = dplyr::case_when(
        data.quality.check.passed.check.noise.x == FALSE ~ data.quality.penalty * quality_penalty,
        TRUE ~ data.quality.penalty
      )) %>%
      dplyr::mutate(data.quality.penalty = dplyr::case_when(
        data.quality.check.passed.check.resolution.y == FALSE ~ data.quality.penalty * quality_penalty,
        TRUE ~ data.quality.penalty
      )) %>%
      dplyr::mutate(data.quality.penalty = dplyr::case_when(
        data.quality.check.passed.check.resolution.x == FALSE ~ data.quality.penalty * quality_penalty,
        TRUE ~ data.quality.penalty
      )) %>%
      dplyr::mutate(offset_raise.weighed = dplyr::case_when(
        offset_raise_required ~ -data.quality.penalty,
        TRUE ~ data.quality.penalty
      ))

    # judge by weighed majority decision if an offset raise is required
    if (mean(optimum_fits$offset_raise.weighed) >= 0) {
      optimum_fits_are_found <- TRUE
    }

    # repeat otherwise
    if (!optimum_fits_are_found) {
      raise_offset_times <- raise_offset_times + 1
    }

    # next offset is out of data range
    if (raise_offset_times > 17) {
      stop("Data is unfit for processing: upper index of the optimum fit region is beyond the last point in the truncated test record.")
    }
  }

  # discard fits when offset raise would be required or numerical problem in fitting occurred
  optimum_fits.nrow <- nrow(optimum_fits)
  optimum_fits <- optimum_fits %>%
    dplyr::filter(.data$offset_raise_required == FALSE) %>%
    dplyr::mutate("norm.residual" = dplyr::case_when(
      .data$norm.residual < 10 * .Machine$double.eps ~ 10 * .Machine$double.eps,
      TRUE ~ .data$norm.residual
    ))

  # calculate success rate of quality checks
  passed.check.data <- nrow(optimum_fits %>%
    dplyr::filter(.data$data.quality.check.passed.check == TRUE)) / optimum_fits.nrow
  passed.check.fit <- nrow(optimum_fits %>%
    dplyr::filter(.data$passed.check == TRUE)) / optimum_fits.nrow
  passed.check <- nrow(optimum_fits %>%
    dplyr::filter(.data$data.quality.check.passed.check == TRUE) %>%
    dplyr::filter(.data$passed.check == TRUE)) / optimum_fits.nrow

  # give a message on success rate of data and fit quality checks
  if (verbose) {
    message("  Random sub-sampling information:")

    message(paste0(
      "      ", optimum.range.size,
      " points of ", nrow(normalized_data$data.normalized),
      " points in the normalized range were used."
    ))
    message(paste0("      ", round(passed.check.data * 100, 1), " % of sub-sampled normalized ranges passed the data quality checks."))
    message(paste0("      ", round(passed.check.fit * 100, 1), " % of linear regressions passed the fit quality checks."))
    message(paste0("      ", round(passed.check * 100, 1), " % of linear regressions passed all quality checks.\n  "))
  }

  # get numerically stable fits and add weights
  optimum_fits <- optimum_fits %>%
    # penalize by fit quality metrics
    dplyr::mutate("fit.quality.penalty" = 1.0) %>%
    dplyr::mutate("fit.quality.penalty" = dplyr::case_when(
      fit.quality.passed.Fit_range == FALSE ~ .data$fit.quality.penalty * quality_penalty,
      TRUE ~ .data$fit.quality.penalty
    )) %>%
    dplyr::mutate("fit.quality.penalty" = dplyr::case_when(
      .data$fit.quality.passed.fourth_quartile == FALSE ~ .data$fit.quality.penalty * quality_penalty,
      TRUE ~ .data$fit.quality.penalty
    )) %>%
    dplyr::mutate("fit.quality.penalty" = dplyr::case_when(
      .data$fit.quality.passed.first_quartile == FALSE ~ .data$fit.quality.penalty * quality_penalty,
      TRUE ~ .data$fit.quality.penalty
    )) %>%
    dplyr::mutate("quality.penalty" = .data$fit.quality.penalty * .data$data.quality.penalty) %>%
    dplyr::mutate("weight" = .data$quality.penalty / .data$norm.residual) %>%
    dplyr::select(-c(
      "offset_raise_required",
      "data.quality.penalty",
      "fit.quality.penalty",
      "quality.penalty",
      "offset_raise.weighed"
    ))

  # use linear regression to find weighed means of results
  # (otr.idx.start and .end, y and x bounds)
  # will be used for reporting the results
  weighed.results <- optimum_fits %>%
    dplyr::select(c(
      "weight", "otr.idx.start", "otr.idx.end",
      "y.lowerBound", "y.upperBound",
      "x.lowerBound", "x.upperBound"
    )) %>%
    as.list() %>%
    {
      res <- .
      weight <- .[["weight"]]
      res$weight <- NULL
      res %>%
        purrr::map(carrier::crate(function(value) {
          list(
            "value" = value,
            "weight" = weight
          )
        }, weight = weight)) %>%
        purrr::map(carrier::crate(function(value) {
          fit <- stats::lm(value ~ 1,
            data = data.frame(
              "value" = value$value,
              "weights" = value$weight
            ),
            weights = weights
          )
          conf <- stats::confint(fit)
          list(
            "value" = fit$coefficients[[1]],
            "conf.low" = conf[[1]],
            "conf.high" = conf[[2]]
          )
        }))
    } %>%
    as.data.frame() %>%
    {
      value <- .
      new.names <- names(value) %>% stringr::str_replace(stringr::coll(".value"), "")
      value %>% magrittr::set_names(new.names)
    }

  # check data quality on the finally selected range
  data_quality_metrics <- check_data_quality(
    data.normalized = normalized_data$data.normalized,
    verbose = verbose.all,
    plot = plot.all,
    plotFun = plotFun.all
  )

  # linear regression on (complete) normalized data
  optimum.otr.idx.start <- round(weighed.results[[1, "otr.idx.start"]], 0)
  optimum.otr.idx.end <- round(weighed.results[[1, "otr.idx.end"]], 0)
  fit.input <- normalized_data$data.normalized %>%
    dplyr::filter(dplyr::between(
      .data$otr.idx,
      optimum.otr.idx.start,
      optimum.otr.idx.end
      ))

  # use basic .lm.fit for normalized residual
  fit.result <- .lm.fit(
    cbind(1, fit.input$x.normalized),
    fit.input$y.normalized
  )

  # calculate normalized residual
  mean_y <- mean(fit.input$y.normalized)
  optimum_fit.norm_residual <- sum((fit.result$residuals)^2) /
    sum((fit.input$y.normalized - mean_y)^2)

  optimum_fit <- data.frame(
    "m" = fit.result$coefficients[2],
    "b" = fit.result$coefficients[1],
    "otr.idx.start" = optimum.otr.idx.start,
    "otr.idx.end" = optimum.otr.idx.end
  )

  # use plain lm() for statistics
  optimum_fit.lm <- stats::lm(
    y.data ~ x.data,
    data = prepared_data[optimum_fit$otr.idx.start:optimum_fit$otr.idx.end,])
  optimum_fit.var <- diag(stats::vcov(optimum_fit.lm))
  optimum_fit.detailed <- optimum_fit %>%
    cbind(data.frame(
      "n.samples" = nrow(
        prepared_data[optimum_fit$otr.idx.start:optimum_fit$otr.idx.end,]),
      "trueIntercept" = optimum_fit.lm$coefficients[[1]],
      "finalSlope" = optimum_fit.lm$coefficients[[2]],
      "trueIntercept.var" = optimum_fit.var[[1]],
      "finalSlope.var" = optimum_fit.var[[2]]
    ))

  # remove remains from calculation
  rm(mean_y, fit.result, fit.input)

  # check fit quality
  fit_quality_metrics <- check_fit_quality(
    data.normalized = normalized_data$data.normalized,
    fit = optimum_fit,
    verbose = verbose.all,
    plot = plot.all,
    plotFun = plotFun.all
  )

  # calculate weight for optimum fit and merge fit on complete data and fits
  # from random sub-sampling
  optimum_fits.mixed <- optimum_fit.detailed %>%
    cbind(data.frame("norm.residual" = optimum_fit.norm_residual)) %>%
    dplyr::mutate("norm.residual" = dplyr::case_when(
      .data$norm.residual < 10 * .Machine$double.eps ~ 10 * .Machine$double.eps,
      TRUE ~ .data$norm.residual
    )) %>%
    # penalize by data quality metrics
    dplyr::mutate("data.quality.penalty" = 1.0) %>%
    dplyr::mutate("data.quality.penalty" = dplyr::case_when(
      data_quality_metrics$digitalResolution$x$passed.check == FALSE ~ .data$data.quality.penalty * quality_penalty,
      TRUE ~ .data$data.quality.penalty
    )) %>%
    dplyr::mutate(data.quality.penalty = dplyr::case_when(
      data_quality_metrics$digitalResolution$y$passed.check == FALSE ~ .data$data.quality.penalty * quality_penalty,
      TRUE ~ .data$data.quality.penalty
    )) %>%
    dplyr::mutate(data.quality.penalty = dplyr::case_when(
      data_quality_metrics$Noise$x$passed.check == FALSE ~ .data$data.quality.penalty * quality_penalty,
      TRUE ~ .data$data.quality.penalty
    )) %>%
    dplyr::mutate(data.quality.penalty = dplyr::case_when(
      data_quality_metrics$Noise$y$passed.check == FALSE ~ .data$data.quality.penalty * quality_penalty,
      TRUE ~ .data$data.quality.penalty
    )) %>%
    # penalize by fit quality metrics
    dplyr::mutate("fit.quality.penalty" = 1.0 ) %>%
    dplyr::mutate("fit.quality.penalty" = dplyr::case_when(
      fit_quality_metrics$Curvature$first_quartile$passed.check == FALSE ~ .data$fit.quality.penalty * quality_penalty,
      TRUE ~ .data$fit.quality.penalty
    )) %>%
    dplyr::mutate("fit.quality.penalty" = dplyr::case_when(
      fit_quality_metrics$Curvature$fourth_quartile$passed.check == FALSE ~ .data$fit.quality.penalty * quality_penalty,
      TRUE ~ .data$fit.quality.penalty
    )) %>%
    dplyr::mutate("fit.quality.penalty" = dplyr::case_when(
      fit_quality_metrics$Fit_range$passed.check == FALSE ~ .data$fit.quality.penalty * quality_penalty,
      TRUE ~ .data$fit.quality.penalty
    )) %>%
    dplyr::mutate("quality.penalty" = .data$fit.quality.penalty * .data$data.quality.penalty) %>%
    dplyr::mutate("weight" = .data$quality.penalty / .data$norm.residual) %>%
    dplyr::select(c(
      "finalSlope",
      "trueIntercept",
      "weight",
      "n.samples",
      "trueIntercept.var",
      "finalSlope.var"
    )) %>%
    rbind(optimum_fits %>%
            dplyr::select(c(
              "finalSlope",
              "trueIntercept",
              "weight",
              "n.samples",
              "trueIntercept.var",
              "finalSlope.var"
            )))

  # calculate sd from std.err
  optimum_fits.mixed <- optimum_fits.mixed %>%
    dplyr::mutate(
      "trueIntercept.sd" = sqrt(.data$n.samples * .data$trueIntercept.var)) %>%
    dplyr::mutate(
      "finalSlope.sd" = sqrt(.data$n.samples * .data$finalSlope.var))

  # calculate pooled sd for trueIntercept
  trueIntercept.sd.pooled <- optimum_fits.mixed %>%
    dplyr::select(c("n.samples", "trueIntercept.sd")) %>%
    dplyr::mutate("nvar" = (.data$n.samples - 1) * .data$trueIntercept.sd) %>% {
      data <- .
      k <- nrow(data)
      sqrt(sum(data$nvar) / (sum(data$n.samples) - k))
    }

  # calculate pooled sd for finalSlope
  finalSlope.sd.pooled <- optimum_fits.mixed %>%
    dplyr::select(c("n.samples", "finalSlope.sd")) %>%
    dplyr::mutate("nvar" = (.data$n.samples - 1) * .data$finalSlope.sd) %>% {
      data <- .
      k <- nrow(data)
      sqrt(sum(data$nvar) / (sum(data$n.samples) - k))
    }

  # use linear regression to find weighed means of (mixed) results
  # (finalSlope and trueIntercept)
  # will be used for reporting the results
  weighed.fit.results <- optimum_fits.mixed %>%
    dplyr::select(c("weight", "finalSlope", "trueIntercept", "n.samples")) %>%
    as.list() %>%
    {
      res <- .
      weight <- .[["weight"]]
      res$weight <- NULL
      res %>%
        purrr::map(carrier::crate(function(value) {
          list(
            "value" = value,
            "weight" = weight
          )
        }, weight = weight)) %>%
        purrr::map(carrier::crate(function(value) {
          fit <- stats::lm(value ~ 1,
                           data = data.frame(
                             "value" = value$value,
                             "weights" = value$weight
                           ),
                           weights = weights
          )
          # conf <- stats::confint(fit)
          list(
            "value" = fit$coefficients[[1]]#,
            #"conf.low" = conf[[1]],
            #"conf.high" = conf[[2]]
          )
        }))
    } %>%
    as.data.frame() %>%
    {
      value <- .
      new.names <- c("finalSlope", "trueIntercept", "n.samples")
      value %>% magrittr::set_names(new.names)
    }

  # get CI from pooled sd
  t.val <- stats::qt(p = 0.025,
                     df = round(weighed.fit.results$n.samples, 0 ) - 1)
  trueIntercept.conf.low <- weighed.fit.results$trueIntercept +
    t.val * trueIntercept.sd.pooled
  trueIntercept.conf.high <- weighed.fit.results$trueIntercept -
    t.val * trueIntercept.sd.pooled
  finalSlope.conf.low <- weighed.fit.results$finalSlope +
    t.val * trueIntercept.sd.pooled
  finalSlope.conf.high <- weighed.fit.results$finalSlope -
    t.val * trueIntercept.sd.pooled

  # un-normalize data and summarize
  assembled_results <- assemble_report(
    normalized_data = normalized_data,
    otr.info = otr.info,
    dataQualityMetrics = data_quality_metrics,
    fit = optimum_fit,
    fitQualityMetrics = fit_quality_metrics,
    verbose = verbose,
    plot = plot,
    plotFun.all = plotFun.all,
    plotFun = plotFun
  )

  # get labels for results
  results.labels <- list(
    "finalSlope" = assembled_results$Slope_Determination_Results$finalSlope.unit,
    "finalSlope.conf.low" = assembled_results$Slope_Determination_Results$finalSlope.unit,
    "finalSlope.conf.high" = assembled_results$Slope_Determination_Results$finalSlope.unit,
    "trueIntercept" = assembled_results$Slope_Determination_Results$trueIntercept.unit,
    "trueIntercept.conf.low" = assembled_results$Slope_Determination_Results$trueIntercept.unit,
    "trueIntercept.conf.high" = assembled_results$Slope_Determination_Results$trueIntercept.unit,
    "y.lowerBound" = assembled_results$Slope_Determination_Results$y.unit,
    "y.upperBound" = assembled_results$Slope_Determination_Results$y.unit,
    "x.lowerBound" = assembled_results$Slope_Determination_Results$x.unit,
    "x.upperBound" = assembled_results$Slope_Determination_Results$x.unit,
    "y.lowerBound.conf.low" = assembled_results$Slope_Determination_Results$y.unit,
    "y.upperBound.conf.low" = assembled_results$Slope_Determination_Results$y.unit,
    "x.lowerBound.conf.low" = assembled_results$Slope_Determination_Results$x.unit,
    "x.upperBound.conf.low" = assembled_results$Slope_Determination_Results$x.unit,
    "y.lowerBound.conf.high" = assembled_results$Slope_Determination_Results$y.unit,
    "y.upperBound.conf.high" = assembled_results$Slope_Determination_Results$y.unit,
    "x.lowerBound.conf.high" = assembled_results$Slope_Determination_Results$x.unit,
    "x.upperBound.conf.high" = assembled_results$Slope_Determination_Results$x.unit,
    "Num.Obs.fit" = "(points in final fit)",
    "Num.Obs.subsampled" = "(points in sub-sampled data range)",
    "Num.Obs.normalized" = "(points in normalized data range)",
    "passed.check" = "(all checks)",
    "passed.check.data" = "(data quality checks)",
    "passed.check.fit" = "(fit quality checks)"
  )

  # assemble results
  sdar.results <- data.frame(
    "finalSlope" = weighed.fit.results$finalSlope,
    "finalSlope.conf.low" = finalSlope.conf.low,
    "finalSlope.conf.high" = finalSlope.conf.high,
    "trueIntercept" = weighed.fit.results$trueIntercept,
    "trueIntercept.conf.low" = trueIntercept.conf.low,
    "trueIntercept.conf.high" = trueIntercept.conf.high,
    "y.lowerBound" = assembled_results$Slope_Determination_Results$y.lowerBound,
    "y.lowerBound.conf.low" = weighed.results[[1, "y.lowerBound.conf.low"]],
    "y.lowerBound.conf.high" = weighed.results[[1, "y.lowerBound.conf.high"]],
    "y.upperBound" = assembled_results$Slope_Determination_Results$y.upperBound,
    "y.upperBound.conf.low" = weighed.results[[1, "y.upperBound.conf.low"]],
    "y.upperBound.conf.high" = weighed.results[[1, "y.upperBound.conf.high"]],
    "x.lowerBound" = assembled_results$Slope_Determination_Results$x.lowerBound,
    "x.lowerBound.conf.low" = weighed.results[[1, "x.lowerBound.conf.low"]],
    "x.lowerBound.conf.high" = weighed.results[[1, "x.lowerBound.conf.high"]],
    "x.upperBound" = assembled_results$Slope_Determination_Results$x.upperBound,
    "x.upperBound.conf.low" = weighed.results[[1, "x.upperBound.conf.low"]],
    "x.upperBound.conf.high" = weighed.results[[1, "x.upperBound.conf.high"]],
    "Num.Obs.fit" = nrow(prepared_data[optimum_fit$otr.idx.start:optimum_fit$otr.idx.end, ]),
    "Num.Obs.subsampled" = optimum.range.size,
    "Num.Obs.normalized" = assembled_results$Data_Quality_Metrics$Num.Obs.normalized,
    "passed.check" = assembled_results$Slope_Determination_Results$passed.check,
    "passed.check.data" = assembled_results$Data_Quality_Metrics$passed.check,
    "passed.check.fit" = assembled_results$Fit_Quality_Metrics$passed.check,
    "method" = "SDAR with random sub-sampling"
  ) %>%
    labelled::set_variable_labels(.labels = results.labels)

  # return full results
  results <- list(
    "sdar" = sdar.results,
    "Data_Quality_Metrics" = assembled_results$Data_Quality_Metrics,
    "Fit_Quality_Metrics" = assembled_results$Fit_Quality_Metrics
  )

  # append plot
  if (plotFun || plotFun.all) {
    results <- results %>% append(list("plots" = list(
      "final.fit" = assembled_results$plots$plot.fit
    )))
  }

  # append more plots
  if (plotFun.all) {
    results$plots <- results$plots %>% append(list(
      "otr" = assembled_results$plots$plot.otr,
      "normalized" = assembled_results$plots$plot.normalized,
      "hist.x" = assembled_results$Data_Quality_Metrics$digitalResolution$plots$hist.x,
      "hist.y" = assembled_results$Data_Quality_Metrics$digitalResolution$plots$hist.y,
      "residuals" = assembled_results$Fit_Quality_Metrics$plots$plot.residuals
    ))

    results$Data_Quality_Metrics$digitalResolution$plots <- NULL
    results$Fit_Quality_Metrics$plots <- NULL
  }

  # return results
  return(results)
}

#' @title Random Sub-Sampling Variant of the SDAR-Algorithm
#'
#' @description Run a random sub-sampling modification of the SDAR algorithm as
#'   standardized in "ASTM E3076-18". As the original version uses numerous
#'   linear regressions (`.lm.fit()` from the stats-package), it can be
#'   painfully slow for test data with high resolution. The lazy variant of the
#'   algorithm will use several random sub-samples of the data to find an
#'   estimate for the fit-range within the data and thus can improve processing
#'   speed. See the article [Speed Benchmarking the SDAR-algorithm](https://soylentorange.github.io/sdarr/articles/speed_improvment.html)
#'   for further information.
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
#'
#' @inheritParams sdar
#'
#' @param n.fit Repetitions of drawing a random sub-sample from the data in the
#'   normalized and finding a fitting range to find an estimate for the final
#'   fitting range.
#'
#' @param cutoff_probability Cut-off probability for estimating optimum size of
#'   sub-sampled data range via logistic regression, which is predicting if
#'   sub-sampled data will pass the quality checks.
#'
#' @param ... <[`dynamic-dots`][rlang::dyn-dots]> Pass parameters to downstream
#'   functions: e.g. set `verbose.all` or `plot.all` to `TRUE` to get
#'   additional diagnostic information during processing data. Set
#'   `enforce_subsampling` to `TRUE` to run the random sub-sampling algorithm
#'   even though it might be slower than the standard SDAR-algorithm.
#'
#' @seealso [sdar()] for the standard SDAR-algorithm.
#'
#' @returns A list containing a data.frame with the results of the final fit,
#'   lists with the quality- and fit-metrics, and a list containing the crated
#'   plot-functions.
#'
#' @examples
#' # Synthesize a test record resembling EN AW-6060-T66.
#' # Explicitly set names to "strain" and "stress".
#'
#' Al_6060_T66 <- synthesize_test_data(
#'   slope = 69000,
#'   yield.y = 160,
#'   ultimate.y = 215,
#'   ultimate.x = 0.08,
#'   x.name = "strain",
#'   y.name = "stress",
#'   toe.start.y = 3, toe.end.y = 10,
#'   toe.start.slope = 13600
#' )
#'
#'
#' # use sdar_lazy() to analyze the (noise-free) synthetic test record
#' # will print a report and give a plot of the final fit
#' \donttest{
#' result <- sdar_lazy(Al_6060_T66, strain, stress)
#' }
#'
#' @export
sdar_lazy <- function(data, x, y,
                      verbose = TRUE,
                      plot = TRUE,
                      n.fit = 5,
                      cutoff_probability = 0.5,
                      ...) {
  # to be furrr-safe, enquote the tidy arguments here
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)

  # @param n.candidates Give a number for selecting optimum fit candidates
  #   (ordered decreasingly by normalized residuals) for each of the repetitions
  #   of random sub-sampling.

  # @param quality_penalty Factor to down-weight fits with inferior data- and
  #   fit-quality metrics.

  # @param enforce_subsampling Set to `TRUE`, to use sub-sampling method even
  #   when it is computationally more expensive than the standard SDAR-algorithm.

  # take care of dynamic dots
  additional_parameters <- rlang::list2(...)
  Nmin_factor <- try({additional_parameters$Nmin_factor}, silent = TRUE)
  if (!is.numeric(Nmin_factor)) Nmin_factor <- 0.2

  verbose.all <- try({additional_parameters$verbose.all}, silent = TRUE)
  if (!is.logical(verbose.all)) verbose.all <- FALSE

  plot.all <- try({additional_parameters$plot.all}, silent = TRUE)
  if (!is.logical(plot.all)) plot.all <- FALSE

  n.candidates <- try({additional_parameters$n.candidates}, silent = TRUE)
  if (!is.numeric(n.candidates)) n.candidates <- 5

  enforce_subsampling <- try({additional_parameters$enforce_subsampling},
                             silent = TRUE)
  if (!is.logical(enforce_subsampling)) enforce_subsampling <- FALSE

  quality_penalty <- try({additional_parameters$quality_penalty}, silent = TRUE)
  if (!is.numeric(quality_penalty)) quality_penalty <- 0.1

  # save final fit plot and all other plots
  plotFun <- TRUE
  plotFun.all <- TRUE

  # give messages for report, when verbose.all is set
  verbose <- verbose || verbose.all

  # show final fit plot, when plot.all is set
  plot <- plot || plot.all

  # get units for data
  x.label.unit <- data %>%
    dplyr::select(!!x) %>%
    labelled::var_label() %>% {
      .[[1]]
    }
  y.label.unit <- data %>%
    dplyr::select(!!y) %>%
    labelled::var_label() %>% {
      .[[1]]
    }

  # get names of data
  x.name <- data %>%
    dplyr::select(!!x) %>%
    names() %>% {
      .[[1]]
    }
  y.name <- data %>%
    dplyr::select(!!y) %>%
    names() %>% {
      .[[1]]
    }

  # prepare data, add an index and remove NA from data
  prepared_data <- data.frame(
    x.data = data %>% dplyr::select(!!x) %>% unlist(TRUE, FALSE),
    y.data = data %>% dplyr::select(!!y) %>% unlist(TRUE, FALSE)) %>%
    dplyr::mutate("otr.idx" = as.numeric(rownames(.))) %>%
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

  # normalize data for initial quality check
  normalized_data <- normalize_data(
    data = prepared_data,
    otr.info = otr.info,
    raise_offset_times = 0,
    verbose = verbose.all,
    plot = plot.all,
    plotFun = plotFun.all
  )

  # find optimum range for sub-sampling
  optimum_size <- optimum_size_for_subsampling(
    data.normalized = normalized_data$data.normalized,
    cutoff_probability = cutoff_probability,
    verbose = verbose.all,
    plot = plot.all,
    plotFun = plotFun.all
  )

  # check for viability of sub-sampling-approach
  viability <- subsampling_viability(
    nrow(normalized_data$data.normalized),
    optimum_size$optimum.range.size,
    n.fit = n.fit,
    verbose = FALSE,
    Nmin_factor = Nmin_factor
  )

  # check, whether sub-sampling would save time and execute standard or lazy variant
  if (viability$viable == FALSE) {
    if (verbose) {
      message(paste0(
        "  lazy algorithm requires more fits than standard SDAR-algorithm:  \n    ",
        viability$Nfits.subsampling, " vs. ", viability$Nfits.plain, " fits."
      ))
      if (enforce_subsampling == FALSE) {
        message("  Standard SDAR-algorithm will be used...\n")
      } else {
        message("  Anyways, random sub-sampling will be used...\n")
      }
    }

    # use standard SDAR-algorithm or force sub-sampling method
    if (enforce_subsampling == FALSE) {
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

  # execute the sub-sampling - SDAR-algorithm
  results <- sdar_execute.lazy(
    prepared_data = prepared_data,
    otr.info = otr.info,
    normalized_data.hint = normalized_data,
    n.fit = n.fit,
    n.candidates = n.candidates,
    optimum.range.size = optimum_size$optimum.range.size,
    cutoff_probability = cutoff_probability,
    quality_penalty = quality_penalty,
    verbose.all = verbose.all,
    verbose = verbose,
    plot.all = plot.all,
    plot = plot,
    plotFun.all = plotFun.all,
    plotFun = plotFun,
    Nmin_factor = Nmin_factor
  )

  # append/sort plots to the results
  if (plotFun.all) {
    results$plots <- results$plots %>% append(list(
      "glm.optimum_size" = optimum_size$plots$plot.modelpredictions
    ))
  }

  # return results
  return(results)
}
