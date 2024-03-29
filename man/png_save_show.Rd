% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/png_save_show.R
\name{png_save_show}
\alias{png_save_show}
\title{Export ggplot as png and include in knitted document}
\usage{
png_save_show(plot, file, show_plot = TRUE, ...)
}
\arguments{
\item{plot}{A ggplot object}

\item{file}{The output path. Filename has to end in png.}

\item{show_plot}{Whether to include the plot in the knitted document.}

\item{...}{Other parameters to be passed to bro_ggsave_paged}
}
\description{
Wrapper for bro::bro_ggsave_paged with unit set to "mm" and automatic inclusion of the png image in the knitted document
}
\examples{
p <- ggplot(mtcars, aes(mpg, cyl)) +
  geom_point()
  png_save_show(p, "~/mtcars.png", width = 89, height = 89)
}
