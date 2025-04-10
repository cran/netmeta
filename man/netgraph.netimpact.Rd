% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netgraph.netimpact.R
\name{netgraph.netimpact}
\alias{netgraph.netimpact}
\title{Network graph for objects of class netimpact}
\usage{
\method{netgraph}{netimpact}(
  x,
  col.ignore = "red",
  number.of.studies = TRUE,
  main,
  sub,
  multiarm = gs("multiarm"),
  col.multiarm = NULL,
  alpha.transparency = 0.5,
  col.ignore.multiarm = "transparent",
  col = "black",
  plastic = gs("plastic"),
  ...
)
}
\arguments{
\item{x}{An object of class \code{netimpact}.}

\item{col.ignore}{A character string indicating color for
comparisons removed from network, either \code{"transparent"} or
any color defined in \code{\link[grDevices]{colours}}.}

\item{number.of.studies}{A logical indicating whether number of
studies should be added to network graph.}

\item{main}{Main title.}

\item{sub}{Subtitle.}

\item{multiarm}{A logical indicating whether multi-arm studies
should be marked in plot.}

\item{col.multiarm}{Either a function from R package
\bold{colorspace} or \bold{grDevice} to define colors for
multi-arm studies or a character vector with colors to highlight
multi-arm studies.}

\item{alpha.transparency}{The alpha transparency of colors used to
highlight multi-arm studies (0 means transparent and 1 means
opaque).}

\item{col.ignore.multiarm}{A character string indicating color to
mark multi-arm studies removed from network, either
\code{"transparent"} or any color defined in
\code{\link[grDevices]{colours}}.}

\item{col}{A single color (or vector of colors) for lines
connecting treatments (edges) if argument \code{plastic =
FALSE}. Length of the vector must be equal to the number of
edges.}

\item{plastic}{A logical indicating whether the appearance of the
comparisons should be in '3D look'.}

\item{\dots}{Additional arguments passed on to
\code{\link{netgraph.netmeta}}.}
}
\description{
This function generates a graph of the evidence network.
}
\examples{
# Examples: example(netimpact)

}
\seealso{
\code{\link{netimpact}}, \code{\link{netgraph.netmeta}},
  \code{\link[metadat]{dat.franchini2012}}
}
\author{
Guido Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de},
  Gerta Rücker \email{gerta.ruecker@uniklinik-freiburg.de}
}
\keyword{hplot}
