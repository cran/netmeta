% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rankogram.netmeta.R
\name{rankogram.netmeta}
\alias{rankogram.netmeta}
\title{Calculate rankogram}
\usage{
\method{rankogram}{netmeta}(
  x,
  nsim = gs("nsim"),
  common = x$common,
  random = x$random,
  small.values = x$small.values,
  cumulative.rankprob = FALSE,
  keep.samples = FALSE,
  nchar.trts = x$nchar.trts,
  warn.deprecated = gs("warn.deprecated"),
  ...
)
}
\arguments{
\item{x}{An object of class \code{\link{netmeta}}.}

\item{nsim}{Number of simulations.}

\item{common}{A logical indicating to compute ranking probabilities
and SUCRAs for the common effects model.}

\item{random}{A logical indicating to compute ranking probabilities
and SUCRAs for the random effects model.}

\item{small.values}{A character string specifying whether small
treatment effects indicate a beneficial (\code{"desirable"}) or
harmful (\code{"undesirable"}) effect, can be abbreviated.}

\item{cumulative.rankprob}{A logical indicating whether cumulative
ranking probabilities should be printed.}

\item{keep.samples}{A logical indicating whether to keep the generated
samples.}

\item{nchar.trts}{A numeric defining the minimum number of
characters used to create unique treatment names.}

\item{warn.deprecated}{A logical indicating whether warnings should
be printed if deprecated arguments are used.}

\item{\dots}{Additional arguments (ignored).}
}
\value{
An object of class \code{rankogram} with corresponding \code{print}
and \code{plot} function. The object is a list containing the
following components:
\item{ranking.matrix.common}{Numeric matrix giving the probability
  of each treatment being at each possible rank for the common
  effects model.}
\item{ranking.common}{SUCRA values for the common effects model.}
\item{ranking.matrix.random}{Numeric matrix giving the probability
  of each treatment being at each possible rank for the random
  effects model.}
\item{ranking.random}{SUCRA values for the random effects model.}
\item{cumrank.matrix.common}{Numeric matrix giving the cumulative
  ranking probability of each treatment for the
  common effects model.}
\item{cumrank.matrix.random}{Numeric matrix giving the cumulative
  ranking probability of each treatment for the random effects
  model.}
\item{nsim, common, random}{As defined above},
\item{small.values, x}{As defined above},
}
\description{
This function calculates the probabilities of each treatment being
at each possible rank and the SUCRAs (Surface Under the Cumulative
RAnking curve) in frequentist network meta-analysis.
}
\details{
We derive a matrix showing the probability of each treatment being
at each possible rank. To this aim, we use resampling from a
multivariate normal distribution with estimated network effects as
means and corresponding estimated variance covariance matrix. We
then summarise them using the ranking metric SUCRAs (Surface Under
the Cumulative RAnking curve).
}
\examples{
pw1 <- pairwise(treatment, event = r, n = N, studlab = author,
  data = dat.woods2010, sm = "OR")
net1 <- netmeta(pw1, small.values = "desirable")

set.seed(1909) # get reproducible results
ran1 <- rankogram(net1, nsim = 100)
ran1
print(ran1, cumulative.rankprob = TRUE)

# Also print mean ranks
summary(ran1)

plot(ran1)

}
\references{
Salanti G, Ades AE, Ioannidis JP (2011):
Graphical methods and numerical summaries for presenting results
from multiple-treatment meta-analysis: an overview and tutorial.
\emph{Journal of Clinical Epidemiology},
\bold{64}, 163--71
}
\seealso{
\code{\link{netmeta}}, \code{\link{netrank}},
  \code{\link{plot.rankogram}}, \code{\link[metadat]{dat.woods2010}}
}
\author{
Theodoros Papakonstantinou \email{dev@tpapak.com}, Guido
  Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
}
