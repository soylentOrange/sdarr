#' Run the SDAR algorithm with random sub-sampling
#' @noRd
sdarr_execute.lazy <- function(prepared_data, otr.info,
                               fit.rep, fit.candidates,
                               optimum.range.size, cutoff_probability,
                               quality_penalty,
                               verbose.all = F,
                               verbose.report = T,
                               showPlots.all = F,
                               showPlots.report = T,
                               savePlots = F) {

  # Give a (long) welcome message
  if(verbose.report) {
    message("Random sub-sampling mofification of the SDAR-algorithm\n")
  }

  # maybe the offset for step 1 needs to be raised later
  raise_offset_times <- 0
  optimum_fits_are_found <- FALSE

  # do step 1, check data quality, and get optimum fit
  # repeat until the upper index is not the last point in the optimum region
  while(!optimum_fits_are_found) {
    # do step 1:
    normalized_data <- normalize_data(prepared_data,
                                      otr.info,
                                      raise_offset_times,
                                      verbose.all,
                                      showPlots.all,
                                      savePlots)

    # find optimum fits and determine summary of data & fit quality metrics
    optimum_fits <- find_linear_regressions.subsampled(normalized_data$data.normalized,
                                                       fit.rep,
                                                       fit.candidates,
                                                       optimum.range.size) %>%
      dplyr::mutate(data.quality.penalty = dplyr::case_when(
        data.quality.check.passed.check == F ~ quality_penalty,
        T ~ 1.0)) %>%
      dplyr::mutate(offset_raise.weighed = dplyr::case_when(
        offset_raise_required ~ -data.quality.penalty,
        T ~ data.quality.penalty))

    # judge by weighed majority decision if an offset raise is required
    if(mean(optimum_fits$offset_raise.weighed) >= 0) {
      optimum_fits_are_found <- T
    }

    # repeat otherwise
    if(!optimum_fits_are_found) {
      raise_offset_times <- raise_offset_times + 1
    }

    # next offset is out of data range
    if(raise_offset_times > 17) {
      stop("Data is unfit for processing: upper index of the optimum fit region is beyond the last point in the truncated test record.")
    }
  }

  # discard fits when offset raise would be required or numerical problem in fitting occurred
  optimum_fits.nrow <- nrow(optimum_fits)
  optimum_fits <- optimum_fits %>%
    dplyr::filter(offset_raise_required == F) %>%
    dplyr::filter(norm.residual > 10*.Machine$double.eps)

  # calculate success rate of quality checks
  passed.check.data <- nrow(optimum_fits %>%
                              dplyr::filter(data.quality.check.passed.check == T))/optimum_fits.nrow
  passed.check.fit <- nrow(optimum_fits %>%
                             dplyr::filter(passed.check == T))/optimum_fits.nrow
  passed.check <- nrow(optimum_fits %>%
                         dplyr::filter(data.quality.check.passed.check == T) %>%
                         dplyr::filter(passed.check == T))/optimum_fits.nrow

  # give a message on success rate of data and fit quality checks
  if(verbose.report) {
    message(paste0("  ", round(passed.check.data*100,1), " % of sub-sampled normalized ranges passed the data quality checks."))
    message(paste0("  ", round(passed.check.fit*100,1), " % of linear regressions passed the fit quality checks."))
    message(paste0("  ", round(passed.check*100,1), " % of linear regressions passed all quality checks."))
  }

  # get numerically stable fits and add weights
  optimum_fits <- optimum_fits %>%
    dplyr::mutate(fit.quality.penalty = dplyr::case_when(
      passed.check == F ~ quality_penalty,
      T ~ 1.0)) %>%
    dplyr::mutate(quality.penalty = fit.quality.penalty * data.quality.penalty) %>%
    dplyr::mutate(weight = quality.penalty/norm.residual) %>%
    dplyr::select(-offset_raise_required,
                  -data.quality.penalty,
                  -fit.quality.penalty,
                  -quality.penalty,
                  -offset_raise.weighed)

  # shorthands for calculations
  y.tangent <- normalized_data[["tangency.point"]][["y.tangent"]]
  x.tangent <- normalized_data[["tangency.point"]][["x.tangent"]]
  x.shift <- normalized_data[["shift"]][["x.shift"]]
  y.shift <- normalized_data[["shift"]][["y.shift"]]

  # Un-normalize fits
  optimum_fits <- optimum_fits %>%
    dplyr::mutate(y.min = normalized_data$data.normalized[otr.idx.start, "y.normalized"]) %>%
    dplyr::mutate(y.max = normalized_data$data.normalized[otr.idx.end, "y.normalized"]) %>%
    dplyr::mutate(x.min = normalized_data$data.normalized[otr.idx.start, "x.normalized"]) %>%
    dplyr::mutate(x.max = normalized_data$data.normalized[otr.idx.end, "x.normalized"]) %>%
    dplyr::mutate(y.lowerBound = y.min * y.tangent + y.shift)  %>%
    dplyr::mutate(y.upperBound = y.max * y.tangent + y.shift)  %>%
    dplyr::mutate(x.lowerBound = x.min * x.tangent + x.shift)  %>%
    dplyr::mutate(x.upperBound = x.max * x.tangent + x.shift)

  # use linear regression to find weighed means of results
  # (otr.idx.start and .end, y and x bounds)
  # will be used for reporting the results and one final fit
  weighed.results <- optimum_fits %>%
    dplyr::select(weight, otr.idx.start, otr.idx.end,
                  y.lowerBound, y.upperBound,
                  x.lowerBound, x.upperBound) %>%
    as.list() %>% {
      res <- .
      weight <- .[["weight"]]
      res$weight <- NULL
      res %>% purrr::map(carrier::crate(function(value) {
        list("value" = value,
             "weight" = weight)
      }, weight = weight)) %>%
        purrr::map(carrier::crate(function(value) {
          fit <- stats::lm(value ~ 1,
                           data = data.frame("value" = value$value,
                                             "weights" = value$weight),
                           weights = weights)
          conf <- stats::confint(fit)
          list("value" = fit$coefficients[[1]],
               "conf.low" = conf[[1]],
               "conf.high" = conf[[2]])
        }))
    } %>% as.data.frame() %>% {
      value <- .
      new.names <- names(value) %>% stringr::str_replace(stringr::coll(".value"), "")
      value %>% magrittr::set_names(new.names)
    }

  # check data quality on the finally selected range
  data_quality_metrics <- check_data_quality(normalized_data$data.normalized,
                                             verbose.all,
                                             showPlots.all,
                                             savePlots)


  # linear regression on (complete) normalized data
  optimum.otr.idx.start <- round(weighed.results[[1, "otr.idx.start"]], 0)
  optimum.otr.idx.end <- round(weighed.results[[1, "otr.idx.end"]], 0)
  optimum.fit <- stats::lm(y.normalized ~ x.normalized,
                           data = normalized_data$data.normalized %>%
                             dplyr::filter(dplyr::between(otr.idx,
                                                          optimum.otr.idx.start,
                                                          optimum.otr.idx.end)))

  optimum_fit <- data.frame(
    "m" = optimum.fit$coefficients[[2]],
    "b" = optimum.fit$coefficients[[1]],
    "otr.idx.start" = optimum.otr.idx.start,
    "otr.idx.end" = optimum.otr.idx.end)

  # check fit quality
  fit_quality_metrics <- check_fit_quality(normalized_data$data.normalized,
                                           optimum_fit,
                                           verbose.all,
                                           showPlots.all,
                                           savePlots)

  # un-normalize data and summarize
  assembled_results <- assemble_report(normalized_data,
                                       otr.info,
                                       data_quality_metrics,
                                       optimum_fit,
                                       fit_quality_metrics,
                                       verbose.report,
                                       showPlots.report,
                                       savePlots)

  # re-do the final fit to find confidence intervals for slope and intercepts
  final.fit.confint <- prepared_data[optimum_fit$otr.idx.start:optimum_fit$otr.idx.end,] %>% {
    prepared_data.fit.lm <- stats::lm(y.data ~ x.data, data = .)
  } %>% stats::confint()

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

  sdarr.results <- data.frame(
    "finalSlope" = assembled_results$Slope_Determination_Results$finalSlope,
    "finalSlope.conf.low" = final.fit.confint["x.data", 1],
    "finalSlope.conf.high" = final.fit.confint["x.data", 2],
    "trueIntercept" = assembled_results$Slope_Determination_Results$trueIntercept,
    "trueIntercept.conf.low" = final.fit.confint["(Intercept)", 1],
    "trueIntercept.conf.high" = final.fit.confint["(Intercept)", 2],
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
    "Num.Obs.fit" = nrow(prepared_data[optimum_fit$otr.idx.start:optimum_fit$otr.idx.end,]),
    "Num.Obs.subsampled" = optimum.range.size,
    "Num.Obs.normalized" = assembled_results$Data_Quality_Metrics$Num.Obs.normalized,
    "passed.check" = assembled_results$Slope_Determination_Results$passed.check,
    "passed.check.data" = assembled_results$Data_Quality_Metrics$passed.check,
    "passed.check.fit" = assembled_results$Fit_Quality_Metrics$passed.check,
    "method" = "SDAR with random sub-sampling"
  ) %>%
    labelled::set_variable_labels(.labels = results.labels)

  # return full results
  results <- list("sdar" = sdarr.results,
                  "Data_Quality_Metrics" = assembled_results$Data_Quality_Metrics,
                  "Fit_Quality_Metrics" = assembled_results$Fit_Quality_Metrics)

  # append/sort plots to the results
  if(savePlots) {
    results <- results %>% append(list("plots" = list(
      "final.fit" = assembled_results$plots$plot.fit,
      "otr" = assembled_results$plots$plot.otr,
      "normalized" = assembled_results$plots$plot.normalized,
      "hist.x" = assembled_results$Data_Quality_Metrics$digitalResolution$plots$hist.x,
      "hist.y" = assembled_results$Data_Quality_Metrics$digitalResolution$plots$hist.y,
      "residuals" = assembled_results$Fit_Quality_Metrics$plots$plot.residuals
    )))

    results$Data_Quality_Metrics$digitalResolution$plots <- NULL
    results$Fit_Quality_Metrics$plots <- NULL
  }

  # return results
  return(results)
}


#' @title sdarr.lazy: Random sub-sampling variant of the SDAR-algorithm
#'
#' @description Run a random sub-sampling modification of the SDAR algorithm as
#'   originally standardized in ASTM E3076-18. As the original version uses
#'   numerous linear regressions (.lm.fit form the stats-package), can be
#'   painfully slow for test data with high resolution. The lazy variant of the
#'   algorithm will use several random sub-samples of the data to find the best
#'   estimate for the fit-range within the data.
#'
#' @note The function will use parallel processing via the furrr-package. To use
#'   this feature, set up a plan other than the default sequential strategy
#'   beforehand.
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
#' @param data Data record to analyze. Labels of the data columns will be used
#'   as units.
#'
#' @param x Use tidy selection to specify x within the data.
#'
#' @param y Use tidy selection to specify y within the data.
#'
#' @param fit.rep Repetitions of random sub-sampling and fitting.
#'
#' @param fit.candidates Give a number for selecting optimum fit candidates
#'   (ordered decreasingly by normalized residuals) for each of the repetitions
#'   of random sub-sampling.
#'
#' @param cutoff_probability Cut-off probability for estimating optimum size of
#'   sub-sampled data range via logistic regression.
#'
#' @param quality_penalty Factor to down-weight fits with inferior data- and
#'   fit-quality metrics.
#'
#' @param enforce_subsampling Set to TRUE, to use sub-sampling method even when
#'   it is computationally more expensive than the standard SDAR-algorithm.
#'
#' @param verbose Give informational messages during computation defaults to
#'   "report" to only show a summarizing information set to "all" to also give
#'   messages from the individual steps set to "none" to be quiet. Can be
#'   abbreviated.
#'
#' @param showPlots Show plots during computation defaults to "report" to only
#'   show the plot of the final fit set to "all" to also show plots from the
#'   individual steps set to "none" to be quiet. Can be abbreviated.
#'
#' @param savePlots Give plot functions with the result.
#'
#' @returns A list containing a data-frame with the results of the final fit, a
#'   list with the quality- and fit-metrics, and a list containing the crated
#'   plot-functions (if savePlots was set to TRUE).
#'
#' @export
sdarr.lazy <- function(data, x, y, fit.rep = 5,
                       fit.candidates = 20,
                       cutoff_probability = 0.975,
                       quality_penalty = 0.1,
                       enforce_subsampling = F,
                       verbose = "report",
                       showPlots = "report",
                       savePlots = F) {

  # to be furrr-safe, enquote the tidy arguments here
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)

  # determine verbosity level
  verbose.all <- ifelse(verbose == "a" || verbose == "all", T, F)
  verbose.report <- ifelse(verbose == "r" || verbose == "report", T, F)
  verbose.report <- verbose.report || verbose.all

  # determine plot level
  showPlots.all <- ifelse(showPlots == "a" || showPlots == "all", T, F)
  showPlots.report <- ifelse(showPlots == "r" || showPlots == "report", T, F)
  showPlots.report <- showPlots.report || showPlots.all

  # get units for data
  x.label.unit <- data %>% dplyr::select(!!x) %>% labelled::var_label() %>% {.[[1]]}
  y.label.unit <- data %>% dplyr::select(!!y) %>% labelled::var_label() %>% {.[[1]]}

  # get names of data
  x.name <- data %>% dplyr::select(!!x) %>% names() %>% {.[[1]]}
  y.name <- data %>% dplyr::select(!!y) %>% names() %>% {.[[1]]}

  # prepare data, add an index and remove NA from data
  prepared_data <- data.frame(x.data = data %>% dplyr::select(!!x) %>% unlist(T,F),
                              y.data = data %>% dplyr::select(!!y) %>% unlist(T,F),
                              otr.idx = seq.int(nrow(data))) %>%
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
  if(verbose.report) {
    message("Determination of Slope in the Linear Region of a Test Record:")
  }

  # normalize data for finding the optimum range
  normalized_data <- normalize_data(prepared_data, otr.info, 0, F, F, F)

  # Check data quality metrics
  if(!check_data_quality.lazy(normalized_data$data.normalized)$passed.check) {
    if(verbose.report) {
      message("  Data quality checks of original test record failed!")
      message("  Examine plots of Standard SDAR-algorithm to determine how to proceed.\n")
    }

    # execute the standard SDAR-algorithm
    result <- sdarr_execute(prepared_data, otr.info, verbose.all, verbose.report,
                            showPlots.all, showPlots.report, savePlots)
    result$sdar <- result$sdar %>% dplyr::mutate("method" = "SDAR")
    return(result)
  }

  # find optimum range for subsampling
  optimum_size <- optimum_size_for_subsampling(normalized_data$data.normalized,
                                               cutoff_probability = cutoff_probability,
                                               showPlots = showPlots.all,
                                               verbose = verbose.all,
                                               savePlots = savePlots)

  # check for viability of subsampling-approach
  viability <- subsampling_viability(nrow(normalized_data$data.normalized), optimum_size$optimum.range.size, fit.rep)

  if(viability$viable == F) {
    if(verbose.report) {
      message(paste0("  lazy algorithm requires more fits than standard SDAR-algorithm:  \n    ",
                     viability$Nfits.subsampling, " vs. ", viability$Nfits.plain, " fits."))
      if(enforce_subsampling == F) {
        message("  Standard SDAR-algorithm will be used...\n")
      } else {
        message("  Anyways, random sub-sampling will be used...\n")
      }
    }

    # use standard SDAR-algorithm or force sub-sampling method
    if(enforce_subsampling == F) {
      # execute the standard SDAR-algorithm
      result <- sdarr_execute(prepared_data, otr.info, verbose.all, verbose.report,
                              showPlots.all, showPlots.report, savePlots)
      result$sdar <- result$sdar %>% dplyr::mutate("method" = "SDAR")
      return(result)
    } else {
      # execute the sub-sampling - SDAR-algorithm
      results <- sdarr_execute.lazy(prepared_data, otr.info, fit.rep, fit.candidates,
                                    optimum_size$optimum.range.size,
                                    cutoff_probability, quality_penalty,
                                    verbose.all, verbose.report,
                                    showPlots.all, showPlots.report, savePlots)

      # append/sort plots to the results
      if(savePlots) {
        results$plots <- results$plots %>% append(list(
          "glm.optimum_size" = optimum_size$plots$plot.modelpredictions
        ))
      }

      # return results
      return(results)
    }
  }
}
