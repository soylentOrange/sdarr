#' Check viability of sub-sampling approach
#'
#' @description Checks whether to use sub-sampling aproach or not by estimating
#' the number required fits for the standard algorithm vs. the sub-sampling
#' approach.
#'
#' @noRd
subsampling_viability <- function(rangeSize.original,
                                  rangeSize.optimum,
                                  n.fit,
                                  verbose = FALSE,
                                  Nmin_factor = 0.2) {
  Nmin.plain <- max(floor(rangeSize.original * 0.2), 10)
  Nfits.plain <- 0.5 * ((rangeSize.original - Nmin.plain + 1)^2 +
    (rangeSize.original - Nmin.plain + 1))

  Nmin.optimum <- max(floor(rangeSize.optimum * Nmin_factor), 10)
  Nfits.subsampling <- n.fit * 0.5 * ((rangeSize.optimum - Nmin.optimum + 1)^2 +
    (rangeSize.optimum - Nmin.optimum + 1))

  list(
    "viable" = (Nfits.subsampling < Nfits.plain),
    "Nfits.plain" = Nfits.plain,
    "Nfits.subsampling" = Nfits.subsampling
  ) %>% return()
}

#' Find optimum size
#'
#' @description
#' Find the optimum size for sub-sampling.
#' Set a random seed beforehand for reproducible results.
#'
#' @noRd
optimum_size_for_subsampling <- function(data.normalized,
                                         cutoff_probability,
                                         verbose = FALSE,
                                         plot = FALSE,
                                         plotFun = FALSE) {

  # if data is rather short, just return the size of data
  if(nrow(data.normalized) <= 50) {
    return(list(
      "optimum.range.size" = nrow(data.normalized),
      "cutoff_probability.adjusted" = cutoff_probability
    ))
  }

  # find optimum fit-range for sub-sampling
  normalized.ranges <- seq.int(
    50,
    nrow(data.normalized),
    (nrow(data.normalized) - 50) / 99
  ) %>%
    round(0) %>%
    {
      data.frame("idx" = seq.int(length(.)), "range.size" = .)
    } %>%
    tidyr::nest(.by = "idx", .key = "range.size") %>%
    dplyr::mutate("otr.idcs" = purrr::map(.data$range.size, carrier::crate(
      function(value) {
        # satisfy pipe addiction...
        `%>%` <- magrittr::`%>%`

        # create data.frame with indices of original test record
        c(
          otr.idx[1],
          otr.idx[length(otr.idx)],
          sample(otr.idx[2:length(otr.idx) - 1], size = value[1, 1] %>%
            as.numeric() - 2)
        ) %>%
          sort() %>%
          {
            data.frame("otr.idx" = .)
          }
      },
      otr.idx = data.normalized$otr.idx %>%
        unlist(TRUE, FALSE)
    ))) %>%
    dplyr::mutate("data.quality.passed" = furrr::future_map(
      .$otr.idcs, carrier::crate(
        function(value) {
          # satisfy pipe addiction...
          `%>%` <- magrittr::`%>%`

          # select data by index
          data.normalized %>%
            dplyr::filter(.data$otr.idx %in% value$otr.idx) %>%
            check_data_quality.lazy() %>%
            {
              1.0 - ifelse(.$passed.check.resolution.y == TRUE, 0, 0.25) -
                ifelse(.$passed.check.resolution.x == TRUE, 0, 0.25) -
                ifelse(.$passed.check.noise.y == TRUE, 0, 0.25) -
                ifelse(.$passed.check.noise.x == TRUE, 0, 0.25)
            }
        },
        data.normalized = data.normalized,
        check_data_quality.lazy = check_data_quality.lazy
      ),
      .options = furrr::furrr_options(globals = FALSE)
    )) %>%
    dplyr::select(c("idx", "range.size", "data.quality.passed")) %>%
    tidyr::unnest(cols = c("range.size", "data.quality.passed")) %>%
    dplyr::select(c("range.size", "data.quality.passed")) %>%
    data.frame()

  # find logistic regression of data quality passed over sub-sampled range-size
  normalized.ranges.model <- stats::glm(data.quality.passed ~ range.size,
    data = normalized.ranges,
    family = stats::quasibinomial(link = "logit")
  )

  # predict probability of passing data quality check over all...
  # ...possible range-sizes
  normalized.ranges.modelpredictions <- data.frame(
    "range.size" = seq.int(50, nrow(data.normalized), 1)
  ) %>%
    {
      probs <- stats::predict(normalized.ranges.model,
        newdata = .,
        type = "link",
        se.fit = TRUE
      )

      ilink <- stats::family(normalized.ranges.model)$linkinv

      dplyr::bind_cols(.,
        pm = ilink(probs$fit),
        pu = ilink(probs$fit + probs$se.fit * 1.96),
        pl = ilink(probs$fit - probs$se.fit * 1.96)
      )
    }

  # check if we can satisfy cutoff_probability,
  # possibly adjust it
  max.pm <- normalized.ranges.modelpredictions$pm %>% max(na.rm = TRUE)
  cutoff_probability <- max.pm * cutoff_probability

  # find optimum range size by (possibly adjusted) cutoff_probability
  optimum.range.size <- normalized.ranges.modelpredictions %>%
    dplyr::filter(.data$pm >= cutoff_probability) %>%
    {
      .[1, "range.size"] %>%
        unlist(TRUE, FALSE)
    }

  # prepare results
  results <- list(
    "optimum.range.size" = optimum.range.size,
    "cutoff_probability.adjusted" = cutoff_probability
  )

  # print message
  if (verbose) {
    message(paste0(
      "  Optimum size for sub-sampling is ",
      optimum.range.size, " (from ",
      nrow(data.normalized), ") samples.\n"
    ))
  }

  # Let's get plotty...
  if (plot || plotFun) {
    # create plot function
    plot.modelpredictions <- carrier::crate(
      function(value) {
        plot(plot.data$range.size,
          plot.data$pm,
          type = "l",
          ylab = plot.ylab,
          xlab = plot.xlab,
          main = plot.main
        )

        graphics::polygon(
          c(
            rev(plot.data$range.size),
            plot.data$range.size
          ),
          c(
            rev(plot.data$pl),
            plot.data$pu
          ),
          col = "grey90", border = NA
        )

        graphics::lines(plot.data$range.size,
          plot.data$pm,
          lwd = 2
        )
        graphics::lines(plot.data$range.size,
          plot.data$pu,
          lwd = 2, col = "red"
        )
        graphics::lines(plot.data$range.size,
          plot.data$pl,
          lwd = 2, col = "red"
        )

        # # add horizontal line for cutoff
        # graphics::abline(h = cutoff_probability, lty= "dotdash")
      },
      plot.data = data.frame(
        "range.size" = normalized.ranges.modelpredictions$range.size %>%
          unlist(TRUE, FALSE),
        pm = normalized.ranges.modelpredictions$pm %>%
          unlist(TRUE, FALSE),
        pu = normalized.ranges.modelpredictions$pu %>%
          unlist(TRUE, FALSE),
        pl = normalized.ranges.modelpredictions$pl %>%
          unlist(TRUE, FALSE)
      ),
      # cutoff_probability = cutoff_probability,
      plot.x = "range.size",
      plot.y = "pm",
      plot.y.upr = "pu",
      plot.y.lwr = "pl",
      plot.ylab = "Probability",
      plot.xlab = "Range Size",
      plot.main = "Probability of passing Data Quality Check"
    )

    if (plotFun) {
      results <- results %>% append(list(
        "plots" = list(
          "plot.modelpredictions" = plot.modelpredictions
        )
      ))
    }

    # Plot Shifted, Truncated and Normalized Test Record
    if (plot) {
      plot.modelpredictions()
    }
  }

  # return results
  return(results)
}
