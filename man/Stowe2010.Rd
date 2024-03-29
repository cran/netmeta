% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Stowe2010.R
\docType{data}
\name{Stowe2010}
\alias{Stowe2010}
\title{Network meta-analysis of adjuvant treatments to levodopa therapy
for Parkinson's disease}
\format{
A data frame with the following columns:
\tabular{rl}{
\bold{\emph{study}}\tab study label \cr
\bold{\emph{id}}\tab study id \cr
\bold{\emph{t1}}\tab treatment 1 \cr
\bold{\emph{y1}}\tab treatment effect arm 1 \cr
\bold{\emph{sd1}}\tab Standard deviation arm 1 \cr
\bold{\emph{n1}}\tab Sample size arm 1 \cr
\bold{\emph{t2}}\tab treatment 2 \cr
\bold{\emph{y2}}\tab treatment effect arm 2 \cr
\bold{\emph{sd2}}\tab Standard deviation arm 2 \cr
\bold{\emph{n2}}\tab Sample size arm 2 \cr
\bold{\emph{t3}}\tab treatment 3 \cr
\bold{\emph{y3}}\tab treatment effect arm 3 \cr
\bold{\emph{sd3}}\tab Standard deviation arm 3 \cr
\bold{\emph{n3}}\tab Sample size arm 3
}
}
\source{
Stowe R, Ives N, Clarke CE, Deane K, Hilten V, Wheatley K, et
al. (2010):
Evaluation of the efficacy and safety of adjuvant treatment to
levodopa therapy in Parkinson's disease patients with motor
complications.
\emph{Cochrane Database of Systematic Reviews}
}
\description{
This data set contains data from a Cochrane review assessing
efficacy and safety of three drug classes as adjuvant treatment to
levodopa therapy in patients with Parkinson’s disease and motor
complications (Stowe et al., 2010). The authors conducted three
pairwise meta-analyses comparing dopamine agonists,
catechol-O-methyl transferase inhibitors (COMTIs), and monoamine
oxidase type B inhibitors (MAOBIs), respectively, with placebo.

The primary outcome was the mean reduction of the time spent in a
relatively immobile ’off’ phase (mean off-time), calculated in
hours per day. Relative treatment effects were expressed as mean
difference. Data on this outcome were available for 5,331 patients
from 28 studies comparing an active treatment with placebo and one
three-arm study comparing two active treatments with placebo.
}
\examples{
data(Stowe2010)

# Transform data from arm-based format to contrast-based format
#
p1 <- pairwise(list(t1, t2, t3),
  n = list(n1, n2, n3),
  mean = list(y1, y2, y3), sd = list(sd1, sd2, sd3),
  data = Stowe2010, studlab = study)
p1

# Conduct network meta-analysis
net1 <- netmeta(p1, ref = "plac")
net1
}
\seealso{
\code{\link{pairwise}}, \code{\link{metacont}},
  \code{\link{netmeta}}, \code{\link{netgraph.netmeta}}
}
\keyword{datasets}
