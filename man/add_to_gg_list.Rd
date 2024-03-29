% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_to_gg_list.R
\name{add_to_gg_list}
\alias{add_to_gg_list}
\title{Add additional ggplot layers to a list of existing ggplot objects}
\usage{
add_to_gg_list(gg_list, arguments_add, subset = NULL)
}
\arguments{
\item{gg_list}{A list of ggplot objects}

\item{arguments_add}{ggplot layers to add provided as a comma-separated list. Do not use + signs between layers.}

\item{subset}{Define a subset of the list to add the layers to, eg. 1:3, or c(1,5,7)}
}
\value{

}
\description{
Add additional ggplot layers to a list of existing ggplot objects
}
