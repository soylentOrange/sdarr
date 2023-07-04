#' report results
#'
#' @note
#' Refer to chapter 7 in ASTM E3076-18
#'
#' @noRd
assemble_report <- function(normalized_data,
                            otr.info = NULL,
                            dataQualityMetrics,
                            fit,
                            fitQualityMetrics,
                            verbose = FALSE,
                            showPlots = FALSE,
                            savePlots = FALSE) {

  # print messages
  if(verbose) {
    message("  Data Quality Metric: Digital Resolution")
    message("    x")
    message(paste0("      Relative x-Resolution:   ",
                   dataQualityMetrics[["digitalResolution"]][["x"]][["Relative_resolution"]]))
    message(paste0("      % at this resolution:    ",
                   dataQualityMetrics[["digitalResolution"]][["x"]][["Percent_at_this_resolution"]]))
    message(paste0("      % in zeroth bin:         ",
                   dataQualityMetrics[["digitalResolution"]][["x"]][["Percent_in_zeroth_bin"]]))
    message(paste0("      --> ",
                   ifelse(dataQualityMetrics[["digitalResolution"]][["x"]][["passed.check"]],
                          "pass",
                          "fail")))

    message("    y")
    message(paste0("      Relative y-Resolution:   ",
                   dataQualityMetrics[["digitalResolution"]][["y"]][["Relative_resolution"]]))
    message(paste0("      % at this resolution:    ",
                   dataQualityMetrics[["digitalResolution"]][["y"]][["Percent_at_this_resolution"]]))
    message(paste0("      % in zeroth bin:         ",
                   dataQualityMetrics[["digitalResolution"]][["y"]][["Percent_in_zeroth_bin"]]))
    message(paste0("      --> ",
                   ifelse(dataQualityMetrics[["digitalResolution"]][["y"]][["passed.check"]],
                          "pass",
                          "fail")), "\n")

    message("  Data Quality Metric: Noise")
    message("    x")
    message(paste0("      Relative x-Noise:        ",
                   dataQualityMetrics[["Noise"]][["x"]][["Relative_noise"]]))
    message(paste0("      --> ",
                   ifelse(dataQualityMetrics[["Noise"]][["x"]][["passed.check"]],
                          "pass",
                          "fail")))

    message("    y")
    message(paste0("      Relative y-Noise:        ",
                   dataQualityMetrics[["Noise"]][["y"]][["Relative_noise"]]))
    message(paste0("      --> ",
                   ifelse(dataQualityMetrics[["Noise"]][["y"]][["passed.check"]],
                          "pass",
                          "fail")), "\n")


    message("  Fit Quality Metric: Curvature")
    message("    1st Quartile")
    message(paste0("      Relative Residual Slope: ",
                   fitQualityMetrics[["Curvature"]][["first_quartile"]][["relative_residual_slope"]]))
    message(paste0("      Number of Points:        ",
                   fitQualityMetrics[["Curvature"]][["first_quartile"]][["Number_of_points"]]))
    message(paste0("      --> ",
                   ifelse(fitQualityMetrics[["Curvature"]][["first_quartile"]][["passed.check"]],
                          "pass",
                          "fail")))

    message("    4th Quartile")
    message(paste0("      Relative Residual Slope: ",
                   fitQualityMetrics[["Curvature"]][["fourth_quartile"]][["relative_residual_slope"]]))
    message(paste0("      Number of Points:        ",
                   fitQualityMetrics[["Curvature"]][["fourth_quartile"]][["Number_of_points"]]))
    message(paste0("      --> ",
                   ifelse(fitQualityMetrics[["Curvature"]][["fourth_quartile"]][["passed.check"]],
                          "pass",
                          "fail")), "\n")

    message("  Fit Quality Metric: Fit Range")
    message(paste0("      relative fit range:      ",
                   fitQualityMetrics[["Fit_range"]][["relative_fit_range"]]))
    message(paste0("      --> ",
                   ifelse(fitQualityMetrics[["Fit_range"]][["passed.check"]],
                          "pass",
                          "fail")), "\n")
  }


  # get values (and possibly print output)
  unnormalized_results <- unnormalize_results(normalized_data,
                                              otr.info,
                                              dataQualityMetrics,
                                              fit,
                                              fitQualityMetrics,
                                              verbose = verbose)

  # Let's get plotty...
  if(showPlots || savePlots) {

    # shorthands
    if(!is.null(otr.info)) {
      # paste names and unit for plotting
      x.lab <- otr.info$x$name
      if(!is.null(otr.info$x$unit)) {
        x.lab <- paste0(x.lab, " (", otr.info$x$unit, ")")
      }

      y.lab <- otr.info$y$name
      if(!is.null(otr.info$y$unit)) {
        y.lab <- paste0(y.lab, " (", otr.info$y$unit, ")")
      }
    } else {
      x.lab <- "x"
      y.lab <- "y"
    }

    # un-normalize data
    x.tangent <- normalized_data[["tangency.point"]][["x.tangent"]]
    y.tangent <- normalized_data[["tangency.point"]][["y.tangent"]]
    x.shift <- normalized_data[["shift"]][["x.shift"]]
    y.shift <- normalized_data[["shift"]][["y.shift"]]
    y.lowerBound <- unnormalized_results[["y.lowerBound"]]
    y.upperBound <- unnormalized_results[["y.upperBound"]]
    finalSlope <- unnormalized_results[["finalSlope"]]
    trueIntercept <- unnormalized_results[["trueIntercept"]]
    data.normalized.unnormalized <- normalized_data[["data.normalized"]] %>%
      dplyr::mutate(x.data = x.normalized * x.tangent + x.shift) %>%
      dplyr::mutate(y.data = y.normalized * y.tangent + y.shift) %>%
      dplyr::mutate(y.fit = x.data * finalSlope + trueIntercept) %>%
      dplyr::select(-otr.idx, -x.normalized, -y.normalized)


    # generate function for plotting the Test Record with fit
    plot.fit <- carrier::crate(function(value) {
      # satisfy pipe addiction...
      `%>%` <- magrittr::`%>%`

      # Plot Test Record with fit
      plot(x = unlist(plot.data$x.data, TRUE, FALSE),
           y = unlist(plot.data$y.data, TRUE, FALSE),
           type ="l",
           col = "red",
           lwd = 1.25,
           main = plot.main,
           xlab = plot.xlab,
           ylab = plot.ylab)

      # add the fit
      data.fit <- plot.data %>%
        dplyr::filter(y.fit <= max(y.data)) %>%
        dplyr::filter(y.fit >= min(y.data))
      graphics::lines(data.fit$x.data, data.fit$y.fit,
                      type ="l",
                      lwd = 1.25,
                      col = "darkgrey")

      # indicate fit range
      graphics::abline(h = c(y.lowerBound, y.upperBound),
                       lty= c("dashed", "dashed"))
      graphics::text(min(plot.data$x.data),
                     (y.upperBound + y.lowerBound)/2,
                     srt = 90,
                     adj = c(0.5, 1),
                     "fit range")
    },
    plot.data = data.frame("x.data" = data.normalized.unnormalized$x.data %>% unlist(TRUE, FALSE),
                           "y.data" = data.normalized.unnormalized$y.data %>% unlist(TRUE, FALSE),
                           "y.fit" = data.normalized.unnormalized$y.fit %>% unlist(TRUE, FALSE)),
    plot.x = "x.data",
    plot.y = "y.data",
    plot.y.fit = "y.fit",
    plot.xlab = x.lab,
    plot.ylab = y.lab,
    plot.main = "Final Fit for the Un-Normalized Test Record",
    y.lowerBound = y.lowerBound,
    y.upperBound = y.upperBound
    )


    # Plot Residuals
    if(showPlots) {
      plot.fit()
    }
  }

  # assemble results and report
  results <- list("Data_Quality_Metrics" = dataQualityMetrics,
                  "Fit_Quality_Metrics" = fitQualityMetrics,
                  "Slope_Determination_Results" = unnormalized_results
  )

  # append plot
  if(savePlots) {
    results <- results %>% append(list("plots" = list(
      "plot.fit" = plot.fit,
      "plot.otr" = normalized_data$plots$plot.otr,
      "plot.normalized" = normalized_data$plots$plot.normalized
    )))
  }

  # return results
  return(results)
}
