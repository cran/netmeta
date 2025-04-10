% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forest.netcomparison.R
\name{forest.netcomparison}
\alias{forest.netcomparison}
\alias{plot.netcomparison}
\title{Forest plot for complex interventions in component network
meta-analysis}
\usage{
\method{forest}{netcomparison}(
  x,
  pooled = ifelse(x$random, "random", "common"),
  leftcols = c("studlab", "treat1", "treat2"),
  leftlabs = c("Comparison", "Trt 1", "Trt 2"),
  rightcols = c("effect", "ci", "statistic", "pval"),
  rightlabs = c(NA, NA, "z", "p-value"),
  nchar.comps = x$nchar.trts,
  digits = gs("digits.forest"),
  digits.stat = gs("digits.stat"),
  digits.pval = gs("digits.pval"),
  smlab = NULL,
  backtransf = x$backtransf,
  lab.NA = gs("lab.NA"),
  equal.size = gs("equal.size"),
  ...
)

\method{plot}{netcomparison}(x, ...)
}
\arguments{
\item{x}{An object of class \code{netcomparison}.}

\item{pooled}{A character string indicating whether results for the
common (\code{"common"}) or random effects model (\code{"random"})
should be plotted. Can be abbreviated.}

\item{leftcols}{A character vector specifying (additional) columns
to be plotted on the left side of the forest plot or a logical
value (see \code{\link[meta]{forest.meta}} help page for details).}

\item{leftlabs}{A character vector specifying labels for
(additional) columns on left side of the forest plot (see
\code{\link[meta]{forest.meta}} help page for details).}

\item{rightcols}{A character vector specifying (additional) columns
to be plotted on the right side of the forest plot or a logical
value (see \code{\link[meta]{forest.meta}} help page for details).}

\item{rightlabs}{A character vector specifying labels for
(additional) columns on right side of the forest plot (see
\code{\link[meta]{forest.meta}} help page for details).}

\item{nchar.comps}{A numeric defining the minimum number of
characters used to create unique names for components.}

\item{digits}{Minimal number of significant digits for treatment
effects and confidence intervals, see \code{print.default}.}

\item{digits.stat}{Minimal number of significant digits for tests
of overall effect, see \code{print.default}.}

\item{digits.pval}{Minimal number of significant digits for p-value
of overall effects, see \code{print.default}.}

\item{smlab}{A label printed at top of figure. By default, text
indicating either common or random effects model is printed.}

\item{backtransf}{A logical indicating whether results should be
back transformed in forest plots. If \code{backtransf = TRUE},
results for \code{sm = "OR"} are presented as odds ratios rather
than log odds ratios, for example.}

\item{lab.NA}{A character string to label missing values.}

\item{equal.size}{A logical indicating whether all squares should
be of equal size. Otherwise, the square size is proportional to
the precision of estimates.}

\item{\dots}{Additional arguments for \code{\link[meta]{forest.meta}}
function.}
}
\description{
Draws a forest plot in the active graphics window (using grid
graphics system).
}
\details{
A forest plot, also called confidence interval plot, is drawn in
the active graphics window. For more information see help page of
\code{\link[meta]{forest.meta}} function.
}
\examples{
# Examples: example(netcomparison)

}
\seealso{
\code{\link{netcomparison}}, \code{\link{netcomb}},
  \code{\link{discomb}}, \code{\link[meta]{forest.meta}}
}
\author{
Guido Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
}
\keyword{hplot}
