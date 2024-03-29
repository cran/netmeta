% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hasse.R
\name{hasse}
\alias{hasse}
\title{Hasse diagram}
\usage{
hasse(x, pooled = ifelse(x$random, "random", "common"), newpage = TRUE)
}
\arguments{
\item{x}{An object of class \code{netposet} (mandatory).}

\item{pooled}{A character string indicating whether Hasse diagram
show be drawn for common (\code{"common"}) or random effects model
(\code{"random"}). Can be abbreviated.}

\item{newpage}{A logical value indicating whether a new figure
should be printed in an existing graphics window. Otherwise, the
Hasse diagram is added to the existing figure.}
}
\description{
This function generates a Hasse diagram for a partial order of
treatment ranks in a network meta-analysis.
}
\details{
Generate a Hasse diagram (Carlsen & Bruggemann, 2014) for a partial
order of treatment ranks in a network meta-analysis (Rücker &
Schwarzer, 2017).

This R function is a wrapper function for R function
\code{\link[hasseDiagram]{hasse}} in R package \bold{hasseDiagram}
(Krzysztof Ciomek, \url{https://github.com/kciomek/hasseDiagram}),
i.e., function \code{hasse} can only be used if R package
\bold{hasseDiagram} is installed.
}
\examples{
\dontrun{
# Use depression dataset
#
data(Linde2015)

# Define order of treatments
#
trts <- c("TCA", "SSRI", "SNRI", "NRI",
  "Low-dose SARI", "NaSSa", "rMAO-A", "Hypericum", "Placebo")

# Outcome labels
#
outcomes <- c("Early response", "Early remission")

# (1) Early response
#
p1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
  event = list(resp1, resp2, resp3),
  n = list(n1, n2, n3),
  studlab = id, data = Linde2015, sm = "OR")
#
net1 <- netmeta(p1, common = FALSE,
  seq = trts, ref = "Placebo", small.values = "undesirable")

# (2) Early remission
#
p2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
  event = list(remi1, remi2, remi3),
  n = list(n1, n2, n3),
  studlab = id, data = Linde2015, sm = "OR")
#
net2 <- netmeta(p2, common = FALSE,
  seq = trts, ref = "Placebo", small.values = "undesirable")

# Partial order of treatment rankings
#
po <- netposet(netrank(net1), netrank(net2), outcomes = outcomes)

# Hasse diagram
#
hasse(po)
}

}
\references{
Carlsen L, Bruggemann R (2014):
Partial order methodology: a valuable tool in chemometrics.
\emph{Journal of Chemometrics},
\bold{28}, 226--34

Rücker G, Schwarzer G (2017):
Resolve conflicting rankings of outcomes in network meta-analysis:
Partial ordering of treatments.
\emph{Research Synthesis Methods},
\bold{8}, 526--36
}
\seealso{
\code{\link{netmeta}}, \code{\link{netposet}},
  \code{\link{netrank}}, \code{\link{plot.netrank}}
}
\author{
Gerta Rücker \email{gerta.ruecker@uniklinik-freiburg.de}, Guido
  Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
}
\keyword{plot}
