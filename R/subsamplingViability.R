#' Check viability of sub-sampling approach
#'
#' @description Checks whether to use sub-sampling aproach or not by estimating
#' the number required fits for the standard algorithm vs. the sub-sampling
#' approach.
#'
#' @noRd
subsampling_viability <- function(rangeSize.original, rangeSize.optimum,
                                  rep = 1, verbose = F) {

  Nmin.plain <- max(floor(rangeSize.original * 0.2), 10)
  Nfits.plain <- 0.5*((rangeSize.original - Nmin.plain + 1)^2 + (rangeSize.original - Nmin.plain + 1))

  Nmin.optimum <- max(floor(rangeSize.optimum * 0.2), 10)
  Nfits.subsampling <- rep * 0.5*((rangeSize.optimum - Nmin.optimum + 1)^2 + (rangeSize.optimum - Nmin.optimum + 1))

  list("viable" = (Nfits.subsampling < Nfits.plain),
       "Nfits.plain" = Nfits.plain,
       "Nfits.subsampling" = Nfits.subsampling
  ) %>% return()
}

#' Find optimum size
#'
#' @description
#' Find the optimum size for sub-sampling.
#' It might be a wise idea to set the random seed beforehand for reproducible results.
#'
#' @noRd
optimum_size_for_subsampling <- function(data.normalized,
                                         cutoff_probability,
                                         verbose = F,
                                         showPlots = F,
                                         savePlots = F) {

  # find optimum fit-range for sub-sampling
  normalized.ranges <- seq.int(50,nrow(data.normalized),(nrow(data.normalized) - 50)/99) %>%
    round(0) %>% {
    data.frame(idx = seq.int(length(.)), range.size = .)
  } %>% tidyr::nest(.by = idx, .key = "range.size") %>%
    dplyr::mutate(otr.idcs = purrr::map(range.size, carrier::crate(function(value) {
      # satisfy pipe addiction...
      `%>%` <- magrittr::`%>%`

      # create data.frame with indices of original test record
      c(otr.idx[1],
        otr.idx[length(otr.idx)],
        sample(otr.idx[2:length(otr.idx)-1], size = value[1,1] %>% as.numeric() - 2)) %>%
        sort() %>% {
          data.frame("otr.idx" = .)
        }
    },
    otr.idx = data.normalized$otr.idx %>%
      unlist(T,F)))) %>%
    dplyr::mutate(data.quality.passed = furrr::future_map(otr.idcs, carrier::crate(function(value) {
      # satisfy pipe addiction...
      `%>%` <- magrittr::`%>%`

      # select data by index
      data.normalized %>%
        dplyr::filter(otr.idx %in% value$otr.idx) %>%
        check_data_quality.lazy() %>% {
          1.0 - ifelse(.$passed.check.resolution.y == T, 0, 0.25) -
            ifelse(.$passed.check.resolution.x == T, 0, 0.25) -
            ifelse(.$passed.check.noise.y == T, 0, 0.25) -
            ifelse(.$passed.check.noise.x == T, 0, 0.25)
        }
    },
    data.normalized = data.normalized,
    check_data_quality.lazy = check_data_quality.lazy),
    .options = furrr::furrr_options(globals = FALSE))) %>%
    dplyr::select(idx, range.size, data.quality.passed) %>%
    tidyr::unnest(c(range.size, data.quality.passed)) %>%
    dplyr::select(range.size, data.quality.passed) %>%
    data.frame()

  # find logistic regression of data quality passed over sub-sampled range-size
  normalized.ranges.model <- stats::glm(data.quality.passed ~ range.size,
                                        data = normalized.ranges,
                                        family = stats::quasibinomial(link='logit'))

  # predict probability of passing data quality check over all possible range-sizes
  normalized.ranges.modelpredictions <- data.frame(range.size = seq.int(50,nrow(data.normalized),1)) %>%  {
    probs <- stats::predict(normalized.ranges.model,
                            newdata = .,
                            type = 'link',
                            se.fit = T)

    ilink <- stats::family(normalized.ranges.model)$linkinv

    dplyr::bind_cols(.,
                     pm = ilink(probs$fit),
                     pu = ilink(probs$fit + probs$se.fit * 1.96),
                     pl = ilink(probs$fit - probs$se.fit * 1.96)
    )
  }

  # check if we can satisfy cutoff_probability,
  # possibly adjust it
  max.pm <- normalized.ranges.modelpredictions$pm %>% max(na.rm = T)
  cutoff_probability.matched <- TRUE
  if(max.pm < cutoff_probability) {
    warning("Failed to determine optimum size by given cutoff_probability.", call. = F)
    cutoff_probability <- max.pm
    cutoff_probability.matched <- FALSE
  }

  # find optimum range size by (possibly adjusted) cutoff_probability
  optimum.range.size <- normalized.ranges.modelpredictions %>%
    dplyr::filter(pm >= cutoff_probability) %>%
    {.[1, "range.size"] %>%
        unlist(T,F)}

  # prepare results
  results <- list("optimum.range.size" = optimum.range.size,
                  "cutoff_probability.matched" = cutoff_probability.matched,
                  "max.pm" = max.pm)

  # Let's get plotty...
  if(showPlots || savePlots) {

    # create plot function
    plot.modelpredictions <- carrier::crate(function(value) {
      plot(plot.data$range.size,
           plot.data$pm,
           type ="l",
           ylab = plot.ylab,
           xlab = plot.xlab,
           main = plot.main)

      graphics::polygon(c(rev(plot.data$range.size),
                          plot.data$range.size),
                        c(rev(plot.data$pl),
                          plot.data$pu),
                        col = "grey90", border = NA)

      graphics::lines(plot.data$range.size,
                      plot.data$pm,
                      lwd = 2)
      graphics::lines(plot.data$range.size,
                      plot.data$pu,
                      lwd = 2, col = "red")
      graphics::lines(plot.data$range.size,
                      plot.data$pl,
                      lwd = 2, col = "red")

      # # add horizontal line for cutoff
      # graphics::abline(h = cutoff_probability, lty= "dotdash")

    }, plot.data = data.frame("range.size" = normalized.ranges.modelpredictions$range.size %>% unlist(T,F),
                              pm = normalized.ranges.modelpredictions$pm %>% unlist(T,F),
                              pu = normalized.ranges.modelpredictions$pu %>% unlist(T,F),
                              pl = normalized.ranges.modelpredictions$pl %>% unlist(T,F)),
    #cutoff_probability = cutoff_probability,
    plot.x = "range.size",
    plot.y = "pm",
    plot.y.upr = "pu",
    plot.y.lwr = "pl",
    plot.ylab = "Probability",
    plot.xlab = "Range Size",
    plot.main = "Probability of passing Data Quality Check")

    if(savePlots) {
      results <- results %>% append(list(
        "plots" = list(
          "plot.modelpredictions" = plot.modelpredictions)))
    }

    # Plot Shifted, Truncated and Normalized Test Record
    if(showPlots) {
      plot.modelpredictions()
    }
  }

  # return results
  return(results)
}
