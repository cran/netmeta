% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forest.netsplit.R
\name{forest.netsplit}
\alias{forest.netsplit}
\alias{plot.netsplit}
\title{Forest plot for direct and indirect evidence}
\usage{
\method{forest}{netsplit}(
  x,
  pooled = ifelse(x$x$random, "random", "common"),
  show = "both",
  subgroup = "comparison",
  overall = TRUE,
  direct = TRUE,
  indirect = TRUE,
  prediction = x$prediction,
  only.reference = FALSE,
  sortvar = NULL,
  subset = NULL,
  text.overall = "Network estimate",
  text.direct = "Direct estimate",
  text.indirect = "Indirect estimate",
  text.predict = "Prediction interval",
  type.overall,
  type.direct,
  type.indirect,
  col.square = "gray",
  col.square.lines = col.square,
  col.inside = "white",
  col.diamond = "gray",
  col.diamond.lines = "black",
  col.predict = "red",
  col.predict.lines = "black",
  equal.size = TRUE,
  leftcols,
  leftlabs,
  rightcols = c("effect", "ci"),
  rightlabs = NULL,
  digits = gs("digits.forest"),
  digits.prop = max(gs("digits.pval") - 2, 2),
  backtransf = x$backtransf,
  lab.NA = "",
  smlab,
  ...
)

\method{plot}{netsplit}(x, ...)
}
\arguments{
\item{x}{An object of class \code{netsplit}.}

\item{pooled}{A character string indicating whether results for the
common (\code{"common"}) or random effects model (\code{"random"})
should be plotted. Can be abbreviated.}

\item{show}{A character string indicating which comparisons should
be printed (see Details).}

\item{subgroup}{A character string indicating which layout should
be used in forest plot: subgroups by comparisons
(\code{"comparison"}) or subgroups by estimates
(\code{"estimate"}). Can be abbreviated.}

\item{overall}{A logical indicating whether network meta-analysis
estimates should be printed.}

\item{direct}{A logical indicating whether direct estimates should
be printed.}

\item{indirect}{A logical indicating whether indirect estimates
should be printed.}

\item{prediction}{A logical indicating whether prediction intervals
should be printed.}

\item{only.reference}{A logical indicating whether only comparisons
with the reference group should be printed.}

\item{sortvar}{An optional vector used to sort comparisons (must be
of same length as the total number of comparisons).}

\item{subset}{An optional logical vector specifying a subset of
comparisons to print (must be of same length as the total number of
comparisons) .}

\item{text.overall}{A character string used in the plot to label
the network estimates.}

\item{text.direct}{A character string used in the plot to label the
direct estimates.}

\item{text.indirect}{A character string used in the plot to label
the indirect estimates.}

\item{text.predict}{A character string used in the plot to label
the prediction interval.}

\item{type.overall}{A character string specifying how to plot
treatment effects and confidence intervals for the overall
network evidence.}

\item{type.direct}{A character string specifying how to plot
treatment effects and confidence intervals for the direct
evidence.}

\item{type.indirect}{A character string specifying how to plot
treatment effects and confidence intervals for the indirect
evidence.}

\item{col.square}{The colour for squares.}

\item{col.square.lines}{The colour for the outer lines of squares.}

\item{col.inside}{The colour for results and confidence limits if
confidence limits are completely within squares squares.}

\item{col.diamond}{The colour of diamonds.}

\item{col.diamond.lines}{The colour of the outer lines of diamonds.}

\item{col.predict}{Background colour of prediction intervals.}

\item{col.predict.lines}{Colour of outer lines of prediction
intervals.}

\item{equal.size}{A logical indicating whether all squares should
be of equal size. Otherwise, the square size is proportional to
the precision of estimates.}

\item{leftcols}{A character vector specifying columns to be plotted
on the left side of the forest plot (see Details).}

\item{leftlabs}{A character vector specifying labels for columns on
left side of the forest plot.}

\item{rightcols}{A character vector specifying columns to be
plotted on the right side of the forest plot (see Details).}

\item{rightlabs}{A character vector specifying labels for columns
on right side of the forest plot.}

\item{digits}{Minimal number of significant digits for treatment
effects and confidence intervals, see \code{print.default}.}

\item{digits.prop}{Minimal number of significant digits for the
direct evidence proportion.}

\item{backtransf}{A logical indicating whether results should be
back transformed in forest plots. If \code{backtransf = TRUE},
results for \code{sm = "OR"} are presented as odds ratios rather
than log odds ratios, for example.}

\item{lab.NA}{A character string to label missing values.}

\item{smlab}{A label printed at top of figure. By default, text
indicating either common or random effects model is printed.}

\item{\dots}{Additional arguments for \code{\link{forest.meta}}
function.}
}
\description{
Forest plot to show direct and indirect evidence in network
meta-analysis.  Furthermore, estimates from network meta-analysis
as well as prediction intervals can be printed.
}
\details{
A forest plot, also called confidence interval plot, is drawn in
the active graphics window.

The arguments \code{leftcols} and \code{rightcols} can be used to
specify columns which are plotted on the left and right side of the
forest plot, respectively. If argument \code{rightcols} is
\code{FALSE}, no columns will be plotted on the right side.

If direct estimates are included in the forest plot (\code{direct =
TRUE}, default), the following columns will be printed on the left
side of the forest plot: the comparisons (column \code{"studlab"}
in \code{\link{forest.meta}}), number of pairwise comparisons
(\code{"k"}), direct evidence proportion (\code{"k"}), and
I\eqn{^2} from pairwise comparison (\code{"I2"}).

If direct estimates are not included in the forest plot
(\code{direct = FALSE}), only the comparisons (\code{"studlab"})
are printed on the left side of the forest plot.

For more information see help page of \code{\link{forest.meta}}
function.

Argument \code{show} determines which comparisons are printed:
\tabular{ll}{
\dQuote{all} \tab All comparisons \cr
\dQuote{both} \tab Only comparisons contributing both direct and
  indirect evidence \cr
\dQuote{with.direct} \tab Comparisons providing direct evidence \cr
\dQuote{direct.only} \tab Comparisons providing only direct
  evidence \cr
\dQuote{indirect.only} \tab Comparisons providing only indirect
  evidence
}
}
\examples{
data(Senn2013)

# Only consider first five studies (to reduce runtime of example)
#
studies <- unique(Senn2013$studlab)
Senn2013.5 <- subset(Senn2013, studlab \%in\% studies[1:5])

net1 <- netmeta(TE, seTE, treat1.long, treat2.long,
  studlab, data = Senn2013.5, common = FALSE, reference = "plac")
#
ns1 <- netsplit(net1)

# Forest plot showing comparisons contributing both direct and
# indirect evidence
#
forest(ns1, fontsize = 6, spacing = 0.5, addrow.subgroups = FALSE)

\dontrun{
# Forest plot showing comparisons contributing direct evidence
#
forest(ns1, fontsize = 6, spacing = 0.5, addrow.subgroups = FALSE,
  show = "with.direct")


# Forest plot only showing network estimates compared to reference
# group and prediction intervals
#
forest(ns1, fontsize = 8, spacing = 0.75, show = "all",
  only.reference = TRUE, prediction = TRUE,
  direct = FALSE, indirect = FALSE)
}

}
\seealso{
\code{\link{forest.meta}}
}
\author{
Guido Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
}
\keyword{hplot}
