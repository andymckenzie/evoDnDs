% Generated by roxygen2 (4.1.0.9000): do not edit by hand
% Please edit documentation in R/get_dn_ds.R
\name{get_dn_ds}
\alias{get_dn_ds}
\title{Dn/Ds Call Function}
\usage{
get_dn_ds(species_name, list_of_genes)
}
\arguments{
\item{species_name}{Ensembl/biomaRt namespace consistent species name.}

\item{list_of_genes}{Character vector of hgnc symbol names to search.}
}
\description{
This function queries biomaRt and ensembl to find the dn/ds values for a character vector of gene names.
}
\examples{
get_dn_ds
}
\keyword{dn/ds}
\keyword{evolution}

