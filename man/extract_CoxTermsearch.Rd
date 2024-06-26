% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoxTermsearch.R
\name{extract_CoxTermsearch}
\alias{extract_CoxTermsearch}
\title{This is a function to extract log HR based on the resuls from CoxTermsearch}
\usage{
extract_CoxTermsearch(fit, lag, latency, knots)
}
\arguments{
\item{fit}{A model from CoxKnotsearch}

\item{lag}{A value of lag time}

\item{latency}{prespecified latency}

\item{knots}{A vector of prespecified knots}
}
\value{
The log HR estimated at that lag time
}
\description{
This is a function to extract log HR based on the resuls from CoxTermsearch
}
\examples{
time_start <- "age_start"
time_end <- "age_end"
status <- "failure"
exposure <- paste0("lag", 0:19)
knots <- 0:19
latency <- 20
fit <- CoxTermsearch(sim_data, time_start, time_end, status, exposure, knots = knots, latency,
adjusted_variable = adjusted_variable, adjusted_model = adjusted_model)
#' extract_CoxTermsearch(fit, 10, latency = latency, knot = knots)
}
