#' Costum ggplot2 theme
#' @export
theme_mk <- theme_bw() + theme(
  text = element_text(size = 7),
  axis.text.x = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  legend.text = element_text(size = 7),
  strip.text = element_text(size = 7),
  plot.title = element_text(hjust = 0.5, size = 7),
  strip.background = element_blank(),
  line = element_line(size = 0.3),
  legend.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)
