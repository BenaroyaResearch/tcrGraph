% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tcrGraphUtils.R
\name{tcrGraphToGraphNEL}
\alias{tcrGraphToGraphNEL}
\title{Convert a tcrGraph obect to a graphNEL object}
\usage{
tcrGraphToGraphNEL(tcrGraph)
}
\arguments{
\item{tcrGraph}{the tcrGraph object to convert}
}
\value{
an object of class 'graphNEL'
}
\description{
From the graphNEL class description: "This is a class of graphs that are 
represented in terms of nodes and an edge list. This is a suitable 
representation for a graph with a large number of nodes and relatively few edges."
}
\details{
Converting a tcrGraph to a graphNEL object enables the use of utility and 
statistics functions from the \code{graph} package.

Note that tcrGraph does not rely on the \code{graph} package, because it is unavailable in CRAN.
However, the graph package may be useful for performing operations on TCR network graphs, and so this
function is included for portability.
}
\seealso{
\code{\link[graph]{graphNEL}}
}
\author{
Mario G Rosasco, \email{mrosasco@benaroyaresearch.org}
}
