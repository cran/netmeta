% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.netconnection.R
\name{summary.netconnection}
\alias{summary.netconnection}
\alias{print.summary.netconnection}
\title{Summary method for objects of class netconnection}
\usage{
\method{summary}{netconnection}(object, ...)

\method{print}{summary.netconnection}(x, ...)
}
\arguments{
\item{object}{An object of class \code{netconnection}.}

\item{\dots}{Additional arguments (passed on to
\code{\link{print.netconnection}}.}

\item{x}{An object of class \code{summary.netconnection}.}
}
\description{
Summary method for objects of class \code{netconnection} to print
list of studies in subnetworks.
}
\examples{
# Artificial example with two subnetworks
#
t1 <- c("G", "B", "B", "D", "A", "F")
t2 <- c("B", "C", "E", "E", "H", "A")
#
nc2 <- netconnection(t1, t2)
print(nc2, details = TRUE)

summary(nc2)

}
\seealso{
\code{\link{netconnection}}
}
\author{
Guido Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
}
