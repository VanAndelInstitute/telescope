% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotGene.R
\name{STPlotGene}
\alias{STPlotGene}
\title{STPlotGene}
\usage{
STPlotGene(x, genes, size = 1, normalize = FALSE)
}
\arguments{
\item{x}{An ST data object such as returned from STLoad}

\item{genes}{gene identifier(s). Identifier type is identified
automatically.}

\item{size}{Size of data in plot}

\item{normalize}{Normalize expression to 0:1?}
}
\description{
Lookup gene data and create 2D heatmap of expression
}
