% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/discomb.R
\name{discomb}
\alias{discomb}
\title{Additive network meta-analysis for combinations of treatments
(disconnected networks)}
\usage{
discomb(
  TE,
  seTE,
  treat1,
  treat2,
  studlab,
  data = NULL,
  subset = NULL,
  inactive = NULL,
  sep.comps = gs("sep.comps"),
  C.matrix,
  sm,
  level = gs("level"),
  level.ma = gs("level.ma"),
  common = gs("common"),
  random = gs("random") | !is.null(tau.preset),
  reference.group,
  baseline.reference = gs("baseline.reference"),
  seq = gs("sep"),
  tau.preset = NULL,
  tol.multiarm = gs("tol.multiarm"),
  tol.multiarm.se = gs("tol.multiarm.se"),
  details.chkmultiarm = gs("details.chkmultiarm"),
  details.chkident = FALSE,
  sep.trts = gs("sep.trts"),
  nchar.comps = gs("nchar.comps"),
  sep.ia = gs("sep.ia"),
  func.inverse = invmat,
  n1 = NULL,
  n2 = NULL,
  event1 = NULL,
  event2 = NULL,
  incr = NULL,
  overall.hetstat = gs("overall.hetstat"),
  backtransf = gs("backtransf"),
  na.unident = gs("na.unident"),
  title = gs("title"),
  keepdata = gs("keepdata"),
  warn = gs("warn"),
  warn.deprecated = gs("warn.deprecated"),
  nchar.trts = nchar.comps,
  ...
)
}
\arguments{
\item{TE}{Estimate of treatment effect, i.e. difference between
first and second treatment (e.g. log odds ratio, mean difference,
or log hazard ratio). Or an R object created with
\code{\link[meta]{pairwise}}.}

\item{seTE}{Standard error of treatment estimate.}

\item{treat1}{Label/Number for first treatment.}

\item{treat2}{Label/Number for second treatment.}

\item{studlab}{An optional - but important! - vector with study
labels (see \code{\link{netmeta}}).}

\item{data}{An optional data frame containing the study
information.}

\item{subset}{An optional vector specifying a subset of studies to
be used.}

\item{inactive}{A character string defining the inactive treatment
component (see Details).}

\item{sep.comps}{A single character to define separator between
treatment components.}

\item{C.matrix}{C matrix (see Details).}

\item{sm}{A character string indicating underlying summary measure,
e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"},
\code{"HR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.}

\item{level}{The level used to calculate confidence intervals for
individual comparisons.}

\item{level.ma}{The level used to calculate confidence intervals
for network estimates.}

\item{common}{A logical indicating whether a common effects /
common effects network meta-analysis should be conducted.}

\item{random}{A logical indicating whether a random effects network
meta-analysis should be conducted.}

\item{reference.group}{Reference treatment (first treatment is used
if argument is missing).}

\item{baseline.reference}{A logical indicating whether results
should be expressed as comparisons of other treatments versus the
reference treatment (default) or vice versa. This argument is
only considered if \code{reference.group} has been specified.}

\item{seq}{A character or numerical vector specifying the sequence
of treatments in printouts.}

\item{tau.preset}{An optional value for the square-root of the
between-study variance \eqn{\tau^2}.}

\item{tol.multiarm}{A numeric for the tolerance for consistency of
treatment estimates in multi-arm studies which are consistent by
design.}

\item{tol.multiarm.se}{A numeric for the tolerance for consistency
of standard errors in multi-arm studies which are consistent by
design. This check is not conducted if the argument is
\code{NULL}.}

\item{details.chkmultiarm}{A logical indicating whether treatment
estimates and / or variances of multi-arm studies with
inconsistent results or negative multi-arm variances should be
printed.}

\item{details.chkident}{A logical indicating whether details on
unidentifiable components should be printed.}

\item{sep.trts}{A character used in comparison names as separator
between treatment labels.}

\item{nchar.comps}{A numeric defining the minimum number of
characters used to create unique names for components (see
Details).}

\item{sep.ia}{A single character to define separator for interactions.}

\item{func.inverse}{R function used to calculate the pseudoinverse
of the Laplacian matrix L (see \code{\link{netmeta}}).}

\item{n1}{Number of observations in first treatment group.}

\item{n2}{Number of observations in second treatment group.}

\item{event1}{Number of events in first treatment group.}

\item{event2}{Number of events in second treatment group.}

\item{incr}{Numerical value added to cell frequencies.}

\item{overall.hetstat}{A logical indicating whether to print heterogeneity
measures.}

\item{backtransf}{A logical indicating whether results should be
back transformed in printouts and forest plots. If
\code{backtransf = TRUE}, results for \code{sm = "OR"} are
presented as odds ratios rather than log odds ratios, for
example.}

\item{na.unident}{A logical indicating whether unidentifiable
components and combinations should be set to missing values.}

\item{title}{Title of meta-analysis / systematic review.}

\item{keepdata}{A logical indicating whether original data(set)
should be kept in netmeta object.}

\item{warn}{A logical indicating whether warnings should be printed
(e.g., if studies are excluded from meta-analysis due to zero
standard errors).}

\item{warn.deprecated}{A logical indicating whether warnings should
be printed if deprecated arguments are used.}

\item{nchar.trts}{Deprecated argument (replaced by
\code{nchar.comps}).}

\item{\dots}{Additional arguments (to catch deprecated arguments).}
}
\value{
An object of classes \code{discomb} and \code{netcomb} with
corresponding \code{print}, \code{summary}, and \code{forest}
functions. The object is a list containing the following
components:
\item{studlab}{Study labels.}
\item{treat1}{Label/Number for first treatment.}
\item{treat2}{Label/Number for second treatment.}
\item{TE}{Estimate of treatment effect, i.e. difference between
  first and second treatment.}
\item{seTE}{Standard error of treatment estimate.}
\item{seTE.adj.common, seTE.adj.random}{Standard error of treatment
  estimate, adjusted for multi-arm studies.}
\item{event1}{Number of events in first treatment group.}
\item{event2}{Number of events in second treatment group.}
\item{n1}{Number of observations in first treatment group.}
\item{n2}{Number of observations in second treatment group.}
\item{k}{Total number of studies.}
\item{m}{Total number of pairwise comparisons.}
\item{n}{Total number of treatments.}
\item{d}{Total number of designs (corresponding to the unique set
  of treatments compared within studies).}
\item{c}{Total number of components.}
\item{trts}{Treatments included in network meta-analysis.}
\item{comps}{Unique list of components present in the network.}
\item{TE.cnma.common, TE.cnma.random}{A vector of length \emph{m} of
  consistent treatment effects estimated by the additive (common and
  random effects) model.}
\item{seTE.cnma.common, seTE.cnma.random}{A vector of length
  \emph{m} with standard errors estimated by the additive (common
  and random effects) model.}
\item{lower.cnma.common, lower.cnma.random}{A vector of length
  \emph{m} of lower confidence interval limits for consistent
  treatment effects estimated by the additive (common and random
  effects) model.}
\item{upper.cnma.common, upper.cnma.random}{A vector of length
  \emph{m} of upper confidence interval limits for consistent
  treatment effects estimated by the additive (common and random
  effects) model.}
\item{statistic.cnma.common, statistic.cnma.random}{A vector of
  length \emph{m} of z-values for the test of an overall effect
  estimated by the additive (common and random effects) model.}
\item{pval.cnma.common, pval.cnma.random}{A vector of length
  \emph{m} of p-values for the test of an overall effect estimated
  by the additive (common and random effects) model.}
\item{TE.common, TE.random}{\emph{n}x\emph{n} matrix with overall
  treatment effects estimated by the additive (common and random
  effects) model.}
\item{seTE.common, seTE.random}{\emph{n}x\emph{n} matrix with
  standard errors estimated by the additive (common and random
  effects) model.}
\item{lower.common, upper.common, lower.random,
  upper.random}{\emph{n}x\emph{n} matrices with lower and upper
  confidence interval limits estimated by the additive (common and
  random effects) model.}
\item{statistic.common, pval.common, statistic.random,
  pval.random}{\emph{n}x\emph{n} matrices with z-values and
  p-values for test of overall effect estimated by the additive
  (common and random effects) model.}
\item{Comp.common, Comp.random}{A vector of component effects (common
  and random effects model).}
\item{seComp.common, seComp.random}{A vector with corresponding
  standard errors (common and random effects model).}
\item{lower.Comp.common, lower.Comp.random}{A vector with lower
  confidence limits for components (common and random effects
  model).}
\item{upper.Comp.common, upper.Comp.random}{A vector with upper
  confidence limits for components (common and random effects
  model).}
\item{statistic.Comp.common, statistic.Comp.random}{A vector with
  z-values for the overall effect of components (common and random
  effects model).}
\item{pval.Comp.common, pval.Comp.random}{A vector with p-values for
  the overall effect of components (common and random effects
  model).}
\item{Comb.common, Comb.random}{A vector of combination effects (common
  and random effects model).}
\item{seComb.common, seComb.random}{A vector with corresponding
  standard errors (common and random effects model).}
\item{lower.Comb.common, lower.Comb.random}{A vector with lower
  confidence limits for combinations (common and random effects
  model).}
\item{upper.Comb.common, upper.Comb.random}{A vector with upper
  confidence limits for combinations (common and random effects
  model).}
\item{statistic.Comb.common, statistic.Comb.random}{A vector with
  z-values for the overall effect of combinations (common and random
  effects model).}
\item{pval.Comb.common, pval.Comb.random}{A vector with p-values for
  the overall effect of combinations (common and random effects
  model).}
\item{Q.additive}{Overall heterogeneity / inconsistency statistic
  (additive model).}
\item{df.Q.additive}{Degrees of freedom for test of heterogeneity /
  inconsistency (additive model).}
\item{pval.Q.additive}{P-value for test of heterogeneity /
  inconsistency (additive model).}
\item{tau}{Square-root of between-study variance (additive model).}
\item{I2}{I-squared (additive model).}
\item{Q.standard}{Overall heterogeneity / inconsistency statistic
  (standard model).}
\item{df.Q.standard}{Degrees of freedom for test of heterogeneity /
  inconsistency (standard model).}
\item{pval.Q.standard}{P-value for test of heterogeneity /
  inconsistency (standard model).}
\item{Q.diff}{Test statistic for difference in goodness of fit
  between standard and additive model.}
\item{df.Q.diff}{Degrees of freedom for difference in goodness of
  fit between standard and additive model.}
\item{pval.Q.diff}{P-value for difference in goodness of fit
  between standard and additive model.}
\item{X.matrix}{Design matrix (\emph{m}x\emph{n}).}
\item{B.matrix}{Edge-vertex incidence matrix (\emph{m}x\emph{n}).}
\item{C.matrix}{As defined above.}
\item{sm}{Summary measure.}
\item{level.ma}{Level for confidence intervals.}
\item{common, random, tau.preset}{As defined above.}
\item{sep.trts}{A character used in comparison names as separator
  between treatment labels.}
\item{nchar.comps}{A numeric defining the minimum number of
  characters used to create unique component names.}
\item{inactive, sep.comps}{As defined above.}
\item{backtransf}{A logical indicating whether results should be
  back transformed in printouts and forest plots.}
\item{title}{Title of meta-analysis / systematic review.}
\item{x}{As defined above.}
\item{call}{Function call.}
\item{version}{Version of R package netmeta used to create
  object.}
}
\description{
Some treatments in a network meta-analysis may be combinations of
other treatments or have common components. The influence of
individual components can be evaluated in an additive network
meta-analysis model assuming that the effect of treatment
combinations is the sum of the effects of its components. This
function implements this additive model in a frequentist way and is
particularly intended for disconnected networks.
}
\details{
Treatments in network meta-analysis (NMA) can be complex
interventions. Some treatments may be combinations of others or
have common components. The standard analysis provided by
\code{\link{netmeta}} is a NMA where all existing (single or
combined) treatments are considered as different nodes in the
network. Exploiting the fact that some treatments are combinations
of common components, an additive component network meta-analysis
(CNMA) model can be used to evaluate the influence of individual
components. This model assumes that the effect of a treatment
combination is the sum of the effects of its components which
implies that common components cancel out in comparisons.

This R function can be used for disconnected networks. Use
\code{\link{netmeta}} and \code{\link{netcomb}} for connected
networks.

The additive CNMA model has been implemented using Bayesian methods
(Mills et al., 2012; Welton et al., 2013). This function implements
the additive model in a frequentist way (Rücker et al., 2020).

The underlying multivariate model is given by

\deqn{\bold{\delta} = \bold{B} \bold{\theta}, \bold{\theta} =
\bold{C} \bold{\beta}}

with
\describe{
\item{\eqn{\bold{\delta}}}{vector of true treatment effects
  (differences) from individual studies,}
\item{\eqn{\bold{B}}}{design matrix describing the structure of the
  network,}
\item{\eqn{\bold{\theta}}}{parameter vector that represents the
  existing combined treatments,}
\item{\eqn{\bold{C}}}{matrix describing how the treatments are
  composed,}
\item{\eqn{\bold{\beta}}}{parameter vector representing the
  treatment components.}
}
All parameters are estimated using weighted least squares
regression.

Argument \code{inactive} can be used to specify a single component
that does not have any therapeutic value. Accordingly, it is
assumed that the treatment effect of the combination of this
component with an additional treatment component is equal to the
treatment effect of the additional component alone.

Argument \code{sep.comps} can be used to specify the separator
between individual components. By default, the matrix \strong{C} is
calculated internally from treatment names. However, it is possible
to specify a different matrix using argument \code{C.matrix}.

By default, component names are not abbreviated in
printouts. However, in order to get more concise printouts,
argument \code{nchar.comps} can be used to define the minimum
number of characters for abbreviated component names (see
\code{\link{abbreviate}}, argument \code{minlength}). R function
\code{\link{treats}} is utilised internally to create abbreviated
component names.
}
\note{
This function calculates effects for individual components and
complex interventions present in the network.

R function \code{\link{netcomplex}} can be used to calculate the
effect for arbitrary complex interventions in a component network
meta-analysis. Furthermore, R function \code{\link{netcomparison}}
can be used to calculate the effect for comparisons of two
arbitrary complex intervention in a component network
meta-analysis.
}
\examples{
# Artificial dataset
#
t1 <- c("A + B", "A + C", "A"    , "A"    , "D", "D", "E")
t2 <- c("C"    , "B"    , "B + C", "A + D", "E", "F", "F")
#
mean    <- c(4.1, 2.05, 0, 0, 0.1, 0.1, 0.05)
se.mean <- rep(0.1, 7)
#
study <- paste("study", c(1:4, 5, 5, 5))
#
dat <- data.frame(mean, se.mean, t1, t2, study,
  stringsAsFactors = FALSE)
#
trts <- c("A", "A + B", "A + C", "A + D",
  "B", "B + C", "C", "D", "E", "F")
#
comps <- LETTERS[1:6]

# Use netconnection() to display network information
#
netconnection(t1, t2, study)

dc1 <- discomb(mean, se.mean, t1, t2, study, seq = trts)
dc1

\donttest{
forest(dc1, ref = "F")

# Define C matrix manually (which will produce the same results)
#
C <- rbind(c(1, 0, 0, 0, 0, 0),  # A
  c(1, 1, 0, 0, 0, 0),  # A + B
  c(1, 0, 1, 0, 0, 0),  # A + C
  c(1, 0, 0, 1, 0, 0),  # A + D
  c(0, 1, 0, 0, 0, 0),  # B
  c(0, 1, 1, 0, 0, 0),  # B + C
  c(0, 0, 1, 0, 0, 0),  # C
  c(0, 0, 0, 1, 0, 0),  # D
  c(0, 0, 0, 0, 1, 0),  # E
  c(0, 0, 0, 0, 0, 1))  # F
#                  
colnames(C) <- comps
rownames(C) <- trts
#
dc2 <- discomb(mean, se.mean, t1, t2, study, seq = trts,
  C.matrix = C)
#
# Compare C matrices
#
all.equal(dc1$C.matrix, dc2$C.matrix)
}

}
\references{
König J, Krahn U, Binder H (2013):
Visualizing the flow of evidence in network meta-analysis and
characterizing mixed treatment comparisons.
\emph{Statistics in Medicine},
\bold{32}, 5414--29

Mills EJ, Thorlund K, Ioannidis JP (2012):
Calculating additive treatment effects from multiple randomized
trials provides useful estimates of combination therapies.
\emph{Journal of Clinical Epidemiology},
\bold{65}, 1282--8

Rücker G, Petropoulou M, Schwarzer G (2020):
Network meta-analysis of multicomponent interventions.
\emph{Biometrical Journal},
\bold{62}, 808--21

Welton NJ, Caldwell DM, Adamopoulos E, Vedhara K (2009):
Mixed treatment comparison meta-analysis of complex interventions:
psychological interventions in coronary heart disease.
\emph{American Journal of Epidemiology},
\bold{169}: 1158--65
}
\seealso{
\code{\link{netcomb}}, \code{\link{forest.netcomb}},
  \code{\link{summary.netcomb}}, \code{\link{netmeta}},
  \code{\link{netconnection}}, \code{\link{netcomplex}},
  \code{\link{netcomparison}}
}
\author{
Gerta Rücker \email{gerta.ruecker@uniklinik-freiburg.de}, Guido
  Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
}
