#' 1st step of SDAR-algorithm: offset, truncation, and normalization
#' @noRd
normalize_data <- function(data, otr.info = NULL,
                           raise_offset_times = 0,
                           verbose = F,
                           showPlots = F,
                           savePlots = F) {

  # find shift for test data
  shift <- data %>% utils::head(1) %>% magrittr::set_names(c("x.shift", "y.shift"))

  # Shift data
  data.offset <- data %>%
    dplyr::mutate(x.data = x.data - shift$x.shift) %>%
    dplyr::mutate(y.data = y.data - shift$y.shift)

  # find offset
  y.max <- max(data.offset$y.data)

  y.offset <- data.offset %>%
    dplyr::filter(y.data > 0.05 * y.max) %>%
    utils::head(1) %>%
    .$y.data + (0.15 + 0.05 * raise_offset_times) * y.max

  x.offset <- data.offset %>%
    dplyr::filter(y.data > 0.05 * y.max) %>%
    utils::head(1) %>%
    .$x.data

  # compensate offset test data
  data.trunc <- data.offset %>%
    dplyr::filter(y.data > y.offset) %>%
    dplyr::mutate(tangent_slope = (y.data - y.offset) / (x.data - x.offset))

  # try to find tangency point
  tangent_slope.peak <- pracma::findpeaks(data.trunc$tangent_slope,
                                          nups = 2,
                                          ndowns = 2,
                                          npeaks = 1,
                                          sortstr = T)

  # retry with relaxed settings, when no tangency point was found
  if(is.null(tangent_slope.peak)) {
    tangent_slope.peak <- pracma::findpeaks(data.trunc$tangent_slope,
                                            nups = 1,
                                            ndowns = 1,
                                            npeaks = 1,
                                            sortstr = T)
  }

  # give an error, when no tangency point was found
  if(is.null(tangent_slope.peak)) {
    stop("No Tangency point identified.")
  }

  # truncate data by tangency point
  x.tangent <- data.trunc[tangent_slope.peak[1,2],]$x.data
  y.tangent <- data.trunc[tangent_slope.peak[1,2],]$y.data
  otr.idx.tangency <- data.trunc[tangent_slope.peak[1,2],]$otr.idx

  # print messages
  if(verbose) {
    message("Normalize Test Data according to NIST Technical Note 205 / ASTM E3076-18\n")
    message("  shift")
    message(paste0("      shift x:                 ", shift$x.shift))
    message(paste0("      shift y:                 ", shift$y.shift))
    message("  tangency point")
    message(paste0("      tangency point x:        ", x.tangent))
    message(paste0("      tangency point y:        ", y.tangent, "\n"))
  }

  # normalize data
  data.normalized <- data.offset %>%
    dplyr::filter(y.data <= y.tangent) %>%
    dplyr::filter(x.data <= x.tangent) %>%
    dplyr::mutate(x.normalized = x.data/x.tangent) %>%
    dplyr::mutate(y.normalized = y.data/y.tangent) %>%
    dplyr::select(-x.data, -y.data)

  # Let's get plotty...
  if(showPlots || savePlots) {

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


    # create function for plot of Original Test Record
    plot.otr <- carrier::crate(function(value) {
      plot(x = unlist(plot.data$x.data, T, F),
           y = unlist(plot.data$y.data, T, F),
           type ="l",
           col = "red",
           lwd = 1.25,
           main = plot.main,
           xlab = plot.xlab,
           ylab = plot.ylab)
    },
    plot.data = data.frame("x.data" = data$x.data %>% unlist(T,F),
                           "y.data" = data$y.data %>% unlist(T,F)),
    plot.x = "x.data",
    plot.y = "y.data",
    plot.xlab = x.lab,
    plot.ylab = y.lab,
    plot.main = "Original Test Record")

    # Plot Original Test Record
    if(showPlots) {
      plot.otr()
    }

    # create function for plot of Shifted, Truncated and Normalized Test Record
    plot.normalized <- carrier::crate(function(value) {
      plot(x = unlist(plot.data$x.normalized, T, F),
           y = unlist(plot.data$y.normalized, T, F),
           type ="l",
           col = "red",
           lwd = 1.25,
           main = plot.main,
           xlab = plot.xlab,
           ylab = plot.ylab)
    },
    plot.data = data.frame("x.normalized" = data.normalized$x.normalized %>% unlist(T, F),
                           "y.normalized" = data.normalized$y.normalized %>% unlist(T, F)),
    plot.x = "x.normalized",
    plot.y = "y.normalized",
    plot.xlab = "normalized x",
    plot.ylab = "normalized y",
    plot.main = "Shifted, Truncated and Normalized Test Record")

    # Plot Shifted, Truncated and Normalized Test Record
    if(showPlots) {
      plot.normalized()
    }
  }

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

  if(savePlots) {
    results <- results %>% append(list(
      "plots" = list(
        "plot.otr" = plot.otr,
        "plot.normalized" = plot.normalized
      )
    ))
  }

  # return results
  results %>% append(list("data.normalized" = data.normalized))
}
