% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_base_functions.R
\name{Nglycan_hypergeometric}
\alias{Nglycan_hypergeometric}
\title{Nglycan_hypergeometric()}
\usage{
Nglycan_hypergeometric(query, universe)
}
\arguments{
\item{query}{should be a character vector containing N-Glycans names, the N-glycans in this list have to be present in the universe list}

\item{universe}{should be a character vector containing N-Glycans names}
}
\value{
a data.frame containing enrichment results.
}
\description{
Allows to perform an hypergeometric test for N-glycan lists (query vs universe), by default the universe is the list of all observed N-glycans in Metaspace
}
\author{
Geremy Clair
}
