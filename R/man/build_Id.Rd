% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indicators.R
\name{build_Id}
\alias{build_Id}
\title{Build identity matrix}
\usage{
build_Id(adj_mat)
}
\arguments{
\item{adj_mat}{The original adjacency matrix.}
}
\value{
An identity matrix \code{Id} in sparse form.
}
\description{
Build the sparse identity matrix \code{Id} from a graph adjacency matrix.
}
\keyword{internal}
