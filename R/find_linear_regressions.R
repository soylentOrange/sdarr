#' 2nd step of SDAR-algorithm: successive linear regressions
#'
#' @note
#' Refer to chapter 6.6 in ASTM E3076-18
#'
#' @description
#' the actual implementation
#'
#' @returns a list of the best candidates
#'
#' @noRd
execute_find_linear_regressions <- function(data.normalized,
                                            tangency.point = NULL,
                                            shift = NULL,
                                            fit.candidates,
                                            Nmin_factor = 0.2) {
  # should at least return 1 fit
  if (fit.candidates < 1) {
    fit.candidates <- 1
  }

  # find minimum size of window for linear regression
  samples <- data.normalized %>% nrow()
  Nmin <- max(floor(samples * Nmin_factor), 10)

  # find number of required fits
  Nfits <- 0.5 * ((samples - Nmin + 1)^2 + (samples - Nmin + 1))

  # how many rounds do we have to consider when setting up the search-ranges?
  Nrounds <- samples - Nmin + 1
  round.idx <- seq(1, Nrounds, 1)
  round.count <- seq(Nrounds, 1, -1)

  # prepare data frame for start and end indices
  # and do the fitting
  # and find fit with lowest normalized residual
  optimum.fits <- data.frame(
    round = round.idx,
    count = round.count,
    round.start = round.idx
  ) %>%
    tidyr::nest(.by = "round", .key = "limits") %>%
    dplyr::mutate("limits" = purrr::map(.data$limits, carrier::crate(\(limit) {
      # satisfy pipe addiction...
      `%>%` <- magrittr::`%>%`

      # get start and count as numbers without much ado...
      limit.start <- limit[1, "round.start"] %>% unlist(TRUE, TRUE)
      limit.count <- limit[1, "count"] %>% unlist(TRUE, TRUE)

      # "calculate" starts indices
      starts <- rep(limit.start, limit.count)

      # calculate end indices
      ends <- seq(from = limit.start + Nmin - 1, by = 1, along.with = starts)

      # store in data.frame
      data.frame(
        start = starts,
        end = ends
      )
    }, Nmin = Nmin))) %>%
    dplyr::select(c("round", "limits")) %>%
    tidyr::unnest(cols = "limits") %>%
    dplyr::select(c("start", "end")) %>%
    dplyr::mutate("fit.idx" = seq.int(nrow(.))) %>%
    tidyr::nest(.by = "fit.idx", .key = "bounds") %>%
    dplyr::mutate("fit" = furrr::future_map(.data$bounds,
      carrier::crate(\(bound) {
        # satisfy pipe addiction...
        `%>%` <- magrittr::`%>%`

        # get indices for fit
        fit.range <- bound %>% unlist(TRUE, FALSE)

        # grab data to analyze
        fit.input <- data.normalized %>%
          dplyr::filter(dplyr::between(
            as.numeric(rownames(.)),
            fit.range[1], fit.range[2]
          ))

        # use basic .lm.fit for increasing speed
        fit.result <- .lm.fit(
          cbind(1, fit.input$x.normalized),
          fit.input$y.normalized
        )

        mean_y <- mean(fit.input$y.normalized)
        norm_residual <- sum((fit.result$residuals)^2) /
          sum((fit.input$y.normalized - mean_y)^2)

        data.frame(
          norm.residual = norm_residual,
          m = fit.result$coefficients[2],
          b = fit.result$coefficients[1],
          otr.idx.start = fit.input$otr.idx[1],
          otr.idx.end = fit.input$otr.idx[length(fit.input$otr.idx)]
        )
      }, data.normalized = data.normalized, .lm.fit = .lm.fit),
      .options = furrr::furrr_options(globals = FALSE)
    )) %>%
    tidyr::unnest(cols = "fit") %>%
    dplyr::select(-c("fit.idx", "bounds")) %>%
    dplyr::arrange(.data$norm.residual) %>%
    utils::head(fit.candidates) %>%
    dplyr::mutate("offset_raise_required" = dplyr::case_when(
      otr.idx.end == samples ~ TRUE,
      TRUE ~ FALSE
    ))

  # un-normalize data and find statistics for the fits
  if(!is.null(tangency.point) && !is.null(shift)) {

    data.unnormalized <- execute_unnormalize.data(
      data.normalized, tangency.point, shift)

    optimum.fits.detailed <- optimum.fits %>%
      dplyr::mutate("rownum" = rownames(.)) %>%
      tidyr::nest(.by = "rownum", .key = "fit") %>%
      dplyr::mutate("fit.details" = furrr::future_map(
        .data$fit, carrier::crate(function(value) {
          # satisfy pipe addiction...
          `%>%` <- magrittr::`%>%`

          # grab data to re-analyze
          fit.input <- data.unnormalized %>%
            dplyr::filter(dplyr::between(
              .data$otr.idx,
              value$otr.idx.start, value$otr.idx.end
            ))

          # use plain lm() for statistics
          fit.result.lm <- stats::lm(
            y.unnormalized ~ x.unnormalized, data = fit.input)
          fit.result.var <- diag(stats::vcov(fit.result.lm))

          # un-normalize fit
          fit.unnormalized <- execute_unnormalize(
            data.normalized, tangency.point, shift, value
          )

          data.frame(
            "n.samples" = nrow(fit.input),
            "trueIntercept" = fit.result.lm$coefficients[[1]],
            "finalSlope" = fit.result.lm$coefficients[[2]],
            "y.lowerBound" = fit.unnormalized[["y.lowerBound"]],
            "y.upperBound" = fit.unnormalized[["y.upperBound"]],
            "x.lowerBound" = fit.unnormalized[["x.lowerBound"]],
            "x.upperBound" = fit.unnormalized[["x.upperBound"]],
            "trueIntercept.var" = fit.result.var[[1]],
            "finalSlope.var" = fit.result.var[[2]]
          )
        }, data.unnormalized = data.unnormalized,
        data.normalized = data.normalized,
        tangency.point = tangency.point,
        shift = shift,
        execute_unnormalize = execute_unnormalize),
        .options = furrr::furrr_options(globals = FALSE))) %>%
      tidyr::unnest(cols = c("fit", "fit.details"))

    # return detailed results
    return(optimum.fits.detailed %>%
             dplyr::select(c("m", "b",
                             "otr.idx.start", "otr.idx.end",
                             "norm.residual", "offset_raise_required",
                             "n.samples", "trueIntercept",
                             "finalSlope",
                             "y.lowerBound",
                             "y.upperBound",
                             "x.lowerBound",
                             "x.upperBound",
                             "trueIntercept.var",
                             "finalSlope.var")))
  } else {
    # return simple results
    return(optimum.fits %>%
             dplyr::select(c("m", "b",
                             "otr.idx.start", "otr.idx.end",
                             "norm.residual", "offset_raise_required")))
  }
}


#' 2nd step of SDAR-algorithm: successive linear regressions
#'
#' @note
#' Refer to chapter 6.6 in ASTM E3076-18
#'
#' @description
#' A wrapper for the actual implementation
#'
#' @noRd
find_linear_regressions <- function(normalized_data,
                                    verbose = FALSE,
                                    Nmin_factor = 0.2) {
  # find optimum fit
  optimum.fit <- execute_find_linear_regressions(
    data.normalized = normalized_data$data.normalized,
    fit.candidates = 1,
    Nmin_factor = Nmin_factor
  )

  # Assemble results
  results <- optimum.fit %>%
    utils::head(1) %>%
    as.list()

  # print messages
  if (verbose) {
    # find minimum size of window for linear regression
    samples <- normalized_data$data.normalized %>% nrow()
    Nmin <- max(floor(samples * Nmin_factor), 10)

    # find number of required fits
    Nfits <- 0.5 * ((samples - Nmin + 1)^2 + (samples - Nmin + 1))

    message("  linear regressions")
    message(paste0(
      "      total points in the search-range:            ",
      samples
    ))
    message(paste0(
      "      --> minimum number of points for linear-fit: ",
      Nmin
    ))
    message(paste0(
      "      --> total number of linear fits:             ",
      Nfits
    ))
    message(paste0(
      "      --> optimum fit found between indices:      ",
      results$otr.idx.start %>% unlist(TRUE, FALSE), " and ",
      results$otr.idx.end %>% unlist(TRUE, FALSE)
    ))
    message(paste0("      m = ", results$m %>% unlist(TRUE, FALSE)))
    message(paste0("      b = ", results$b %>% unlist(TRUE, FALSE)))
    message(paste0(
      "      --> Upper index of the optimum region is ",
      ifelse(results$offset_raise_required %>% unlist(TRUE, FALSE),
        "", "not "
      ),
      "the last point in the truncated test record\n"
    ))
  }

  # return results
  return(results)
}

#' 2nd step of SDAR-algorithm: successive linear regressions
#'
#' @note
#' Refer to chapter 6.6 in ASTM E3076-18
#'
#' @description
#' A wrapper for the actual implementation in a lazy variant
#'
#' @noRd
find_linear_regressions.subsampled <- function(normalized_data,
                                               fit.rep,
                                               fit.candidates,
                                               range.size,
                                               Nmin_factor) {
  # do the fitting (several times) on sub-sampled data
  best.fits <- seq.int(fit.rep) %>%
    purrr::map(carrier::crate(
      function(value) {
        # satisfy pipe addiction...
        `%>%` <- magrittr::`%>%`

        # create data.frame with indices of original test record
        c(
          otr.idx[1],
          otr.idx[length(otr.idx)],
          sample(otr.idx[2:length(otr.idx) - 1], size = size - 2)
        ) %>%
          sort() %>%
          {
            data.frame("mc.idx" = value, "otr.idx" = .)
          }
      },
      size = range.size,
      otr.idx = normalized_data$data.normalized$otr.idx %>%
        unlist(TRUE, FALSE)
    )) %>%
    purrr::list_rbind() %>%
    tidyr::nest(.by = "mc.idx", .key = "otr.idcs") %>%
    dplyr::mutate("data.quality.check" = furrr::future_map(.data$otr.idcs,
      carrier::crate(
        function(value) {
          # satisfy pipe addiction...
          `%>%` <- magrittr::`%>%`

          # select data by index
          normalized_data$data.normalized %>%
            dplyr::filter(.data$otr.idx %in% value$otr.idx) %>%
            check_data_quality.lazy() %>%
            as.data.frame()
        },
        normalized_data = normalized_data,
        check_data_quality.lazy = check_data_quality.lazy
      ),
      .options = furrr::furrr_options(globals = FALSE)
    )) %>%
    dplyr::mutate("fit" = purrr::map(.data$otr.idcs, carrier::crate(
      function(value) {
        # satisfy pipe addiction...
        `%>%` <- magrittr::`%>%`

        # select data by index
        data.normalized.fits <- data.normalized %>%
          dplyr::filter(.data$otr.idx %in% value$otr.idx)

        # do linear regressions
        fits <- data.normalized.fits %>%
          execute_find_linear_regressions(
            tangency.point = tangency.point,
            shift = shift,
            fit.candidates = fit.candidates,
            Nmin_factor = Nmin_factor
          ) %>%
          split(rownames(.)) %>%
          purrr::map(carrier::crate(
            function(value) {
              # satisfy pipe addiction...
              `%>%` <- magrittr::`%>%`

              fit_quality_metrics <- check_fit_quality.lazy(
                data.normalized.fits, value
              ) %>% as.data.frame()
              value <- value %>%
                dplyr::bind_cols(fit_quality_metrics)
            },
            data.normalized.fits = data.normalized.fits,
            check_fit_quality.lazy = check_fit_quality.lazy
          )) %>%
          purrr::list_rbind()
      },
      data.normalized = normalized_data$data.normalized,
      execute_find_linear_regressions = execute_find_linear_regressions,
      fit.candidates = fit.candidates,
      check_fit_quality.lazy = check_fit_quality.lazy,
      Nmin_factor = Nmin_factor,
      tangency.point = normalized_data$tangency.point,
      shift = normalized_data$shift
    ))) %>%
    dplyr::select(c("mc.idx", "data.quality.check", "fit")) %>%
    tidyr::unnest("fit") %>%
    tidyr::unnest("data.quality.check", names_sep = ".")

  # return best fits
  return(best.fits)
}
