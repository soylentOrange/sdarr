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

  # Assemble results
  data.frame(
    "m" = optimum.fits$m %>% unlist(TRUE, FALSE),
    "b" = optimum.fits$b %>% unlist(TRUE, FALSE),
    "otr.idx.start" = optimum.fits$otr.idx.start %>% unlist(TRUE, FALSE),
    "otr.idx.end" = optimum.fits$otr.idx.end %>% unlist(TRUE, FALSE),
    "norm.residual" = optimum.fits$norm.residual %>% unlist(TRUE, FALSE),
    "offset_raise_required" = optimum.fits$offset_raise_required %>%
      unlist(TRUE, FALSE)
  ) %>%
    return()
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
find_linear_regressions <- function(data.normalized,
                                    verbose = FALSE,
                                    Nmin_factor = 0.2) {
  # find optimum fit
  optimum.fit <- execute_find_linear_regressions(
    data.normalized = data.normalized,
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
    samples <- data.normalized %>% nrow()
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
find_linear_regressions.subsampled <- function(data.normalized,
                                               fit.rep = 10,
                                               fit.candidates = 10,
                                               range.size,
                                               Nmin_factor = 0.2) {
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
      otr.idx = data.normalized$otr.idx %>%
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
          data.normalized %>%
            dplyr::filter(.data$otr.idx %in% value$otr.idx) %>%
            check_data_quality.lazy() %>%
            as.data.frame()
        },
        data.normalized = data.normalized,
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
      data.normalized = data.normalized,
      execute_find_linear_regressions = execute_find_linear_regressions,
      fit.candidates = fit.candidates,
      check_fit_quality.lazy = check_fit_quality.lazy,
      Nmin_factor = Nmin_factor
    ))) %>%
    dplyr::select(c("mc.idx", "data.quality.check", "fit")) %>%
    tidyr::unnest("fit") %>%
    tidyr::unnest("data.quality.check", names_sep = ".")

  # return best fits
  return(best.fits)
}
