#' Universal ggplot export function from github.com/jbengler/bro with minor fix for current ggplot versions
#'
#' This function takes a ggplot (for single page) or list of ggplots (for multi page) and writes them to file.
#' In case the input has absolute dimensions (e.g. as a result of `egg::set_panel_size()` or
#' `patchwork::plot_layout()`), width and height of the output are adjusted to fit the content.
#'
#' @param gg ggplot or list of ggplots
#' @param filename Filename for export
#' @param device Device to use. Can either be a device function (e.g. `png()`), or one of "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf" (windows only).
#' @param path Path to save plot to (combined with filename).
#' @param scale Multiplicative scaling factor.
#' @param width,height,units Plot size in `units` ("mm", "cm", or "in").
#'   If not supplied, uses the size of current graphics device.
#'   In case the input has absolute dimensions after applying
#'   `patchwork::plot_layout()`, width and height of the output are adjusted to fit the content.
#' @param dpi Plot resolution. Also accepts a string input: "retina" (320), "print" (300), or "screen" (72). Applies only to raster output types.
#' @param limitsize When TRUE (the default), ggsave will not save images larger than 50x50 inches, to prevent the common error of specifying dimensions in pixels.
#' @param ... Other arguments passed on to the graphics device function, as specified by device.
#' @param return_input If `TRUE` the input ggplot or plotlist is returned after saving.
#' This enables the use of `bro_ggsave_paged()` within `dplyr` pipes.
#' @param burst_to_multiple_files If `TRUE` multiple plots are saved to separate files instead of a single multi-page file.
#' @export
bro_ggsave_paged <- function(gg = ggplot2::last_plot(), filename, device = NULL, path = NULL, scale = 1,
                             width = NA, height = NA, units = "mm", dpi = 300, limitsize = TRUE,
                             return_input = FALSE, burst_to_multiple_files = FALSE, ...) {
  if (class(gg)[1] %in% c("patchwork", "gg", "ggplot")) gg <- list(gg)
  dimensions <- bro::bro_get_ggsize(gg, units)

  width_defined_by <- dplyr::case_when(
    is.na(width) && is.na(dimensions[["width"]]) ~ "was not defined - default device used",
    !is.na(width) ~ "was provided as parameter 'width' to bro_ggsave_paged()",
    TRUE ~ "was inferred from plot width"
  )
  height_defined_by <- dplyr::case_when(
    is.na(height) && is.na(dimensions[["height"]]) ~ "was not defined - default device used",
    !is.na(height) ~ "was provided as parameter 'height' to bro_ggsave_paged()",
    TRUE ~ "was inferred from plot height"
  )

  if (is.na(width)) width <- dimensions[["width"]]
  if (is.na(height)) height <- dimensions[["height"]]

  message("Device width ", width_defined_by)
  message("Device height ", height_defined_by)
  if (!is.na(width) && !is.na(height)) message("Saving ", round(width), " x ", round(height), " mm image")

  if (burst_to_multiple_files) {
    filenames <- bro::burst_filename(filename, length(gg))
    purrr::map2(gg, filenames, function(x, y) {
      ggplot2::ggsave(
        plot = x, filename = y, device = device, path = path, scale = scale,
        width = width, height = height, units = units, dpi = dpi, limitsize = limitsize, ...
      )
    })
    if (return_input) return(gg)

  } else {

    # the rest of the code is adapted from ggplot2::ggsave()
    dpi <- ggplot2:::parse_dpi(dpi)

    # Build full filename before device selection (device guessing uses file extension)
    if (!is.null(path)) {
      filename <- file.path(path, filename)
    }

    # ggplot2 >= 4 uses validate_device(); older versions used plot_dev()
    dev <- if (exists("validate_device", envir = asNamespace("ggplot2"), inherits = FALSE)) {
      ggplot2:::validate_device(device, filename, dpi = dpi)
    } else {
      ggplot2:::plot_dev(device, filename, dpi = dpi)
    }

    # ggplot2 >= 4 plot_dim() has a dpi argument (relevant for units = "px")
    dim <- if ("dpi" %in% names(formals(ggplot2:::plot_dim))) {
      ggplot2:::plot_dim(c(width, height), scale = scale, units = units, limitsize = limitsize, dpi = dpi)
    } else {
      ggplot2:::plot_dim(c(width, height), scale = scale, units = units, limitsize = limitsize)
    }

    old_dev <- grDevices::dev.cur()
    dev(filename = filename, width = dim[1], height = dim[2], ...)
    on.exit(utils::capture.output({
      grDevices::dev.off()
      if (old_dev > 1) grDevices::dev.set(old_dev)
    }))

    purrr::map(gg, ~grid::grid.draw(.x))
    invisible()
    if (return_input) return(gg)
  }
}

