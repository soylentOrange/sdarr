#' 1st step of SDAR-algorithm: offset, truncation, and normalization
#' actual implementation of normalization
#' @noRd
normalize_data.execute <- function(data, otr.info = NULL,
                                   raise_offset_times = 0,
                                   return_data.normalized = FALSE) {
  # find shift for test data
  shift <- data %>%
    utils::head(1) %>%
    magrittr::set_names(c("x.shift", "y.shift"))

  # Shift data
  data.offset <- data %>%
    dplyr::mutate("x.data" = .data$x.data - shift$x.shift) %>%
    dplyr::mutate("y.data" = .data$y.data - shift$y.shift)

  # find offset
  y.max <- max(data.offset$y.data, na.rm = TRUE)

  y.offset <- data.offset %>%
    dplyr::filter(.data$y.data > 0.05 * y.max) %>%
    utils::head(1) %>%
    .$y.data + (0.15 + 0.05 * raise_offset_times) * y.max

  x.offset <- data.offset %>%
    dplyr::filter(.data$y.data > 0.05 * y.max) %>%
    utils::head(1) %>%
    .$x.data

  # compensate offset test data
  data.trunc <- data.offset %>%
    dplyr::filter(.data$y.data > y.offset) %>%
    dplyr::mutate(tangent_slope = (.data$y.data - y.offset) / (.data$x.data - x.offset))

  # try to find tangency point
  tangent_slope.peak <- pracma::findpeaks(data.trunc$tangent_slope,
    nups = 2,
    ndowns = 2,
    npeaks = 1,
    sortstr = TRUE
  )

  # retry with relaxed settings, when no tangency point was found
  if (is.null(tangent_slope.peak)) {
    tangent_slope.peak <- pracma::findpeaks(data.trunc$tangent_slope,
      nups = 1,
      ndowns = 1,
      npeaks = 1,
      sortstr = TRUE
    )
  }

  # give an error, when no tangency point was found
  if (is.null(tangent_slope.peak)) {
    stop("No Tangency point identified.")
  }

  # truncate data by tangency point
  x.tangent <- data.trunc[tangent_slope.peak[1, 2], ]$x.data
  y.tangent <- data.trunc[tangent_slope.peak[1, 2], ]$y.data
  otr.idx.tangency <- data.trunc[tangent_slope.peak[1, 2], ]$otr.idx

  # assemble results
  results <- list(
    "tangency.point" = list(
      "x.tangent" = x.tangent,
      "y.tangent" = y.tangent,
      "otr.idx.tangent" = otr.idx.tangency
    ),
    "shift" = list(
      "x.shift" = shift$x.shift,
      "y.shift" = shift$y.shift
    )
  )

  if (return_data.normalized) {
    # normalize data
    data.normalized <- data.offset %>%
      dplyr::filter(.data$y.data <= y.tangent) %>%
      dplyr::filter(.data$x.data <= x.tangent) %>%
      dplyr::mutate(x.normalized = .data$x.data / x.tangent) %>%
      dplyr::mutate(y.normalized = .data$y.data / y.tangent) %>%
      dplyr::select(-c("x.data", "y.data"))

    # return results with data
    results %>%
      append(list("data.normalized" = data.normalized)) %>%
      return()
  } else {
    # return bare results
    return(results)
  }
}

#' 1st step of SDAR-algorithm: offset, truncation, and normalization
#' @noRd
normalize_data <- function(data,
                           otr.info = NULL,
                           raise_offset_times = 0,
                           denoise.x = FALSE,
                           denoise.y = FALSE,
                           alpha = 2000,
                           K = 100,
                           verbose = FALSE,
                           plot = FALSE,
                           plotFun = FALSE) {
  # just normalize and return normalized data (with info and (possibly) plots)
  if (denoise.x == FALSE && denoise.y == FALSE) {
    normalized_data <- normalize_data.execute(
      data = data,
      otr.info = otr.info,
      raise_offset_times = raise_offset_times,
      return_data.normalized = TRUE
    )
  } else {
    # do a first normalization on the data to get info on the range
    normalized_data <- normalize_data.execute(
      data = data,
      otr.info = otr.info,
      raise_offset_times = raise_offset_times,
      return_data.normalized = FALSE
    )

    # denoise x and y (in a limited range of the original data, saving time...)
    data_denoised <- data %>%
      dplyr::filter(.data$otr.idx <= normalized_data$tangency.point$otr.idx.tangent * 2)


    # de-noise x-data
    if (denoise.x == TRUE) {
      data_denoised$x.data <- denoise_vector(
        data_denoised$x.data %>%
          as.numeric(),
        alpha = alpha,
        K = K,
        verbose = verbose
      )
    }

    # de-noise y-data
    if (denoise.y == TRUE) {
      data_denoised$y.data <- denoise_vector(
        data_denoised$y.data %>%
          as.numeric(),
        alpha = alpha,
        K = K,
        verbose = verbose
      )
    }

    # do the normalization on the de-noised data
    normalized_data <- normalize_data.execute(
      data = data_denoised,
      otr.info = otr.info,
      raise_offset_times = raise_offset_times,
      return_data.normalized = TRUE
    )
  }

  # print messages
  if (verbose) {
    if (denoise.x || denoise.y) {
      message("Normalize de-noised Test Data\n")
    } else {
      message("Normalize Test Data\n")
    }

    message("  shift")
    message(paste0(
      "      shift x:                 ",
      normalized_data$shift$x.shift
    ))
    message(paste0(
      "      shift y:                 ",
      normalized_data$shift$y.shift
    ))
    message("  tangency point")
    message(paste0(
      "      tangency point x:        ",
      normalized_data$tangency.point$x.tangent
    ))
    message(paste0(
      "      tangency point y:        ",
      normalized_data$tangency.point$y.tangent, "\n"
    ))
  }

  # Let's get plotty...
  if (plot || plotFun) {
    if (!is.null(otr.info)) {
      # paste names and unit for plotting
      x.lab <- otr.info$x$name
      if (!is.null(otr.info$x$unit)) {
        x.lab <- paste0(x.lab, " (", otr.info$x$unit, ")")
      }

      y.lab <- otr.info$y$name
      if (!is.null(otr.info$y$unit)) {
        y.lab <- paste0(y.lab, " (", otr.info$y$unit, ")")
      }
    } else {
      x.lab <- "x"
      y.lab <- "y"
    }


    # create function for plot of Original Test Record
    plot.otr <- carrier::crate(
      function(value) {
        plot(
          x = unlist(plot.data$x.data, TRUE, FALSE),
          y = unlist(plot.data$y.data, TRUE, FALSE),
          type = "l",
          col = "red",
          lwd = 1.25,
          main = plot.main,
          xlab = plot.xlab,
          ylab = plot.ylab
        )
      },
      plot.data = data.frame(
        "x.data" = data$x.data %>% unlist(TRUE, FALSE),
        "y.data" = data$y.data %>% unlist(TRUE, FALSE)
      ),
      plot.x = "x.data",
      plot.y = "y.data",
      plot.xlab = x.lab,
      plot.ylab = y.lab,
      plot.main = "Original Test Record"
    )

    # Plot Original Test Record
    if (plot) {
      plot.otr()
    }

    # create function for plot of de-noised Original Test Record
    if (denoise.y || denoise.x) {
      plot.otr.denoised <- carrier::crate(
        function(value) {
          plot(
            x = unlist(plot.data$x.data.denoised, TRUE, FALSE),
            y = unlist(plot.data$y.data.denoised, TRUE, FALSE),
            type = "l",
            col = "red",
            lwd = 1.25,
            main = plot.main,
            xlab = plot.xlab,
            ylab = plot.ylab
          )
        },
        plot.data = data.frame(
          "x.data.denoised" = data_denoised$x.data %>% unlist(TRUE, FALSE),
          "y.data.denoised" = data_denoised$y.data %>% unlist(TRUE, FALSE)
        ),
        plot.x = "x.data.denoised",
        plot.y = "y.data.denoised",
        plot.xlab = x.lab,
        plot.ylab = y.lab,
        plot.main = "De-noised Original Test Record (limited range)"
      )

      # Plot Original Test Record
      if (plot) {
        plot.otr.denoised()
      }
    }

    # create function for plot of Shifted, Truncated and Normalized Test Record
    plot.normalized <- carrier::crate(
      function(value) {
        plot(
          x = unlist(plot.data$x.normalized, TRUE, FALSE),
          y = unlist(plot.data$y.normalized, TRUE, FALSE),
          type = "l",
          col = "red",
          lwd = 1.25,
          main = plot.main,
          xlab = plot.xlab,
          ylab = plot.ylab
        )
      },
      plot.data = data.frame(
        "x.normalized" = normalized_data$data.normalized$x.normalized %>% unlist(TRUE, FALSE),
        "y.normalized" = normalized_data$data.normalized$y.normalized %>% unlist(TRUE, FALSE)
      ),
      plot.x = "x.normalized",
      plot.y = "y.normalized",
      plot.xlab = "normalized x",
      plot.ylab = "normalized y",
      plot.main = "Shifted, Truncated and Normalized Test Record"
    )

    # Plot Shifted, Truncated and Normalized Test Record
    if (plot) {
      plot.normalized()
    }
  }

  # assemble results
  results <- normalized_data

  if (plotFun) {
    results <- results %>% append(list(
      "plots" = list(
        "plot.otr" = plot.otr,
        "plot.normalized" = plot.normalized
      )
    ))

    # append plot of de-noised (partial) data, if available
    if (denoise.y || denoise.x) {
      results$plots <- results$plots %>% append(list(
        "plot.otr.denoised" = plot.otr.denoised
      ))
    }
  }

  # return results
  return(results)
}
