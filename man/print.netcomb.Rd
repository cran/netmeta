% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.netcomb.R
\name{print.netcomb}
\alias{print.netcomb}
\title{Print method for objects of class netcomb}
\usage{
\method{print}{netcomb}(
  x,
  common = x$common,
  random = x$random,
  overall.hetstat = x$overall.hetstat,
  backtransf = x$backtransf,
  nchar.comps = x$nchar.comps,
  digits = gs("digits"),
  digits.stat = gs("digits.stat"),
  digits.pval = gs("digits.pval"),
  digits.pval.Q = max(gs("digits.pval.Q"), 2),
  digits.Q = gs("digits.Q"),
  digits.tau2 = gs("digits.tau2"),
  digits.tau = gs("digits.tau"),
  digits.I2 = gs("digits.I2"),
  big.mark = gs("big.mark"),
  scientific.pval = gs("scientific.pval"),
  zero.pval = gs("zero.pval"),
  JAMA.pval = gs("JAMA.pval"),
  print.tau2 = gs("print.tau2"),
  print.tau = gs("print.tau"),
  print.Q = gs("print.Q"),
  print.I2 = gs("print.I2"),
  print.I2.ci = gs("print.I2.ci"),
  text.tau2 = gs("text.tau2"),
  text.tau = gs("text.tau"),
  text.I2 = gs("text.I2"),
  details.methods = gs("details"),
  legend = gs("legend"),
  warn.deprecated = gs("warn.deprecated"),
  ...
)
}
\arguments{
\item{x}{An object of class \code{netcomb} or
\code{summary.netcomb}.}

\item{common}{A logical indicating whether results for the common
effects model should be printed.}

\item{random}{A logical indicating whether results for the random
effects model should be printed.}

\item{overall.hetstat}{A logical indicating whether to print heterogeneity
measures.}

\item{backtransf}{A logical indicating whether results should be
back transformed in printouts and forest plots. If
\code{backtransf = TRUE}, results for \code{sm = "OR"} are
presented as odds ratios rather than log odds ratios, for
example.}

\item{nchar.comps}{A numeric defining the minimum number of
characters used to create unique names for components.}

\item{digits}{Minimal number of significant digits, see
\code{print.default}.}

\item{digits.stat}{Minimal number of significant digits for z- or
t-value, see \code{print.default}.}

\item{digits.pval}{Minimal number of significant digits for p-value
of overall treatment effect, see \code{print.default}.}

\item{digits.pval.Q}{Minimal number of significant digits for
p-value of heterogeneity tests, see \code{print.default}.}

\item{digits.Q}{Minimal number of significant digits for
heterogeneity statistics, see \code{print.default}.}

\item{digits.tau2}{Minimal number of significant digits for
between-study variance, see \code{print.default}.}

\item{digits.tau}{Minimal number of significant digits for square
root of between-study variance, see \code{print.default}.}

\item{digits.I2}{Minimal number of significant digits for I-squared
statistic, see \code{print.default}.}

\item{big.mark}{A character used as thousands separator.}

\item{scientific.pval}{A logical specifying whether p-values should
be printed in scientific notation, e.g., 1.2345e-01 instead of
0.12345.}

\item{zero.pval}{A logical specifying whether p-values should be
printed with a leading zero.}

\item{JAMA.pval}{A logical specifying whether p-values for test of
overall effect should be printed according to JAMA reporting
standards.}

\item{print.tau2}{A logical specifying whether between-study
variance \eqn{\tau^2} should be printed.}

\item{print.tau}{A logical specifying whether \eqn{\tau}, the
square root of the between-study variance \eqn{\tau^2}, should be
printed.}

\item{print.Q}{A logical value indicating whether to print the
results of the test of heterogeneity.}

\item{print.I2}{A logical specifying whether heterogeneity
statistic I\eqn{^2} should be printed.}

\item{print.I2.ci}{A logical specifying whether confidence interval for
heterogeneity statistic I\eqn{^2} should be printed.}

\item{text.tau2}{Text printed to identify between-study variance
\eqn{\tau^2}.}

\item{text.tau}{Text printed to identify \eqn{\tau}, the square
root of the between-study variance \eqn{\tau^2}.}

\item{text.I2}{Text printed to identify heterogeneity statistic
I\eqn{^2}.}

\item{details.methods}{A logical specifying whether details on statistical
methods should be printed.}

\item{legend}{A logical indicating whether a legend should be
printed.}

\item{warn.deprecated}{A logical indicating whether warnings should
be printed if deprecated arguments are used.}

\item{\dots}{Additional arguments (to catch deprecated arguments).}
}
\description{
Print method for objects of class \code{netcomb}.
}
\examples{
# Examples: example(netcomb)

}
\seealso{
\code{\link{netcomb}}, \code{\link{discomb}}
}
\author{
Guido Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
}
\keyword{print}
