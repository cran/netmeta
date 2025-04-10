% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netcomb.R
\name{netcomb}
\alias{netcomb}
\title{Additive network meta-analysis for combinations of treatments}
\usage{
netcomb(
  x,
  inactive = NULL,
  sep.comps = gs("sep.comps"),
  C.matrix,
  common = x$common,
  random = x$random | !is.null(tau.preset),
  tau.preset = NULL,
  details.chkident = FALSE,
  nchar.comps = x$nchar.trts,
  sep.ia = gs("sep.ia"),
  func.inverse = invmat,
  overall.hetstat = x$overall.hetstat,
  backtransf = x$backtransf,
  na.unident = gs("na.unident"),
  warn.deprecated = gs("warn.deprecated"),
  ...
)
}
\arguments{
\item{x}{An object of class \code{netmeta}.}

\item{inactive}{A character string defining the inactive treatment
component (see Details).}

\item{sep.comps}{A single character to define separator between
treatment components.}

\item{C.matrix}{C matrix (see Details).}

\item{common}{A logical indicating whether a common effects network
meta-analysis should be conducted.}

\item{random}{A logical indicating whether a random effects network
meta-analysis should be conducted.}

\item{tau.preset}{An optional value for the square-root of the
between-study variance \eqn{\tau^2}.}

\item{details.chkident}{A logical indicating whether details on
unidentifiable components should be printed.}

\item{nchar.comps}{A numeric defining the minimum number of
characters used to create unique names for components (see
Details).}

\item{sep.ia}{A single character to define separator for interactions.}

\item{func.inverse}{R function used to calculate the pseudoinverse
of the Laplacian matrix L (see \code{\link{netmeta}}).}

\item{overall.hetstat}{A logical indicating whether to print heterogeneity
measures.}

\item{backtransf}{A logical indicating whether results should be
back transformed in printouts and forest plots.}

\item{na.unident}{A logical indicating whether unidentifiable
components and combinations should be set to missing values.}

\item{warn.deprecated}{A logical indicating whether warnings should
be printed if deprecated arguments are used.}

\item{\dots}{Additional arguments (to catch deprecated arguments).}
}
\value{
An object of class \code{netcomb} with corresponding \code{print},
\code{summary}, and \code{forest} functions. The object is a list
containing the following components:
\item{studlab}{Study labels.}
\item{treat1}{Label/Number for first treatment.}
\item{treat2}{Label/Number for second treatment.}
\item{TE}{Estimate of treatment effect, i.e. difference between
  first and second treatment.}
\item{seTE}{Standard error of treatment estimate.}
\item{seTE.adj.common, seTE.adj.random}{Standard error of treatment
  estimate, adjusted for multi-arm studies.}
\item{design}{Design of study providing pairwise comparison.}
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
\item{k.trts}{Number of studies evaluating a treatment.}
\item{n.trts}{Number of observations receiving a treatment (if
  arguments \code{n1} and \code{n2} are provided).}
\item{events.trts}{Number of events observed for a treatment (if
  arguments \code{event1} and \code{event2} are provided).}
\item{studies}{Study labels coerced into a factor with its levels
  sorted alphabetically.}
\item{narms}{Number of arms for each study.}
\item{designs}{Unique list of designs present in the network. A
  design corresponds to the set of treatments compared within a
  study.}
\item{comps}{Unique list of components present in the network.}
\item{TE.nma.common, TE.nma.random}{A vector of length \emph{m} of
  consistent treatment effects estimated by network meta-analysis
  (nma) (common and random effects model).}
\item{seTE.nma.common, seTE.nma.random}{A vector of length \emph{m}
  of effective standard errors estimated by network meta-analysis
  (common and random effects model).}
\item{lower.nma.common, lower.nma.random}{A vector of length
  \emph{m} of lower confidence interval limits for consistent
  treatment effects estimated by network meta-analysis (common
  and random effects model).}
\item{upper.nma.common, upper.nma.random}{A vector of length
  \emph{m} of upper confidence interval limits for the consistent
  treatment effects estimated by network meta-analysis (common
  and random effects model).}
\item{statistic.nma.common, statistic.nma.random}{A vector of length \emph{m}
  of z-values for test of treatment effect for individual
  comparisons (common and random effects model).}
\item{pval.nma.common, pval.nma.random}{A vector of length \emph{m}
  of p-values for test of treatment effect for individual
  comparisons (common and random effects model).}
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
\item{statistic.cnma.common, statistic.cnma.random}{A vector of length
  \emph{m} of z-values for the test of an overall effect estimated
  by the additive (common and random effects) model.}
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
\item{statistic.Comp.common, statistic.Comp.random}{A vector with z-values for
  the overall effect of components (common and random effects
  model).}
\item{pval.Comp.common, pval.Comp.random}{A vector with p-values for
  the overall effect of components (common and random effects
  model).}
\item{Comb.common, Comb.random}{A vector of combination effects
  (common and random effects model).}
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
\item{I2, lower.I2, upper.I2}{I-squared, lower and upper confidence
  limits.}
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
\item{A.matrix}{Adjacency matrix (\emph{n}x\emph{n}).}
\item{B.matrix}{Edge-vertex incidence matrix (\emph{m}x\emph{n}).}
\item{C.matrix}{As defined above.}
\item{sm}{Summary measure.}
\item{level.ma}{Level for confidence intervals.}
\item{common, random, tau.preset}{As defined above.}
\item{sep.trts}{A character used in comparison names as separator
  between treatment labels.}
\item{nchar.comps}{A numeric defining the minimum number of
  characters used to create unique component names.}
\item{inactive, sep.comps, sep.ia}{As defined above.}
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
function implements this additive model in a frequentist way.
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

The additive CNMA model has been implemented using Bayesian methods
(Mills et al., 2012; Welton et al., 2013). This function implements
the additive model in a frequentist way (Rücker et al., 2020).

The underlying multivariate model is given by

\deqn{\bold{\delta} = \bold{B} \bold{\theta}, \bold{\theta} =
\bold{C} \bold{\beta}}

with \describe{
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
data(Linde2016)

# Only consider studies including Face-to-face PST (to reduce
# runtime of example)
#
face <- subset(Linde2016, id \%in\% c(16, 24, 49, 118))

# Conduct random effects network meta-analysis
#
net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
  data = face, ref = "placebo", sm = "OR", common = FALSE)
net1

# Additive model for treatment components (with placebo as inactive
# treatment)
#
nc1 <- netcomb(net1, inactive = "placebo")
nc1

# Print results with two digits and test statistics with three digits
print(nc1, digits = 2, digits.stat = 3)

# Print detailed results
print(summary(nc1), digits = 2, digits.stat = 3)

# All available combinations in CNMA
combinations(nc1) # only "Face-to-face PST + SSRI"

# Create C matrix with all available interactions, i.e., one interaction
# for each combination
colnames(createC(nc1))

# Run interaction CNMA model with all available interactions
# (same result as standard NMA)
netcomb(net1, C.matrix = createC(nc1))

\donttest{
forest(net1, xlim = c(0.2, 50))
forest(nc1, xlim = c(0.2, 50))

# Specify, order of treatments
#
trts <- c("TCA", "SSRI", "SNRI", "NRI", "Low-dose SARI", "NaSSa",
  "rMAO-A", "Ind drug", "Hypericum", "Face-to-face CBT",
  "Face-to-face PST", "Face-to-face interpsy", "Face-to-face psychodyn",
  "Other face-to-face", "Remote CBT", "Self-help CBT", "No contact CBT",
  "Face-to-face CBT + SSRI", "Face-to-face interpsy + SSRI",
  "Face-to-face PST + SSRI", "UC", "Placebo")
#
# Note, three treatments are actually combinations of 'SSRI' with
# other components:
# "Face-to-face CBT + SSRI",
# "Face-to-face interpsy + SSRI",
# "Face-to-face PST + SSRI"

# Conduct random effects network meta-analysis
#
net2 <- netmeta(lnOR, selnOR, treat1, treat2, id,
  data = Linde2016, ref = "placebo",
  seq = trts, sm = "OR", common = FALSE)
net2
forest(net2, xlim = c(0.2, 50))

# Additive model for treatment components (with placebo as inactive
# treatment)
#
nc2 <- netcomb(net2, inactive = "placebo")
nc2
forest(nc2, xlim = c(0.2, 50))
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
\code{\link{discomb}}, \code{\link{netmeta}},
  \code{\link{forest.netcomb}}, \code{\link{print.netcomb}},
  \code{\link{netcomplex}}, \code{\link{netcomparison}}
}
\author{
Gerta Rücker \email{gerta.ruecker@uniklinik-freiburg.de}, Guido
  Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
}
