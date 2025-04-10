% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Baker2009.R
\docType{data}
\name{Baker2009}
\alias{Baker2009}
\title{Network meta-analysis of pharmacologic treatments for chronic
obstructive pulmonary disease}
\format{
A data frame with the following columns:
\tabular{rl}{
\bold{\emph{study}}\tab study label \cr
\bold{\emph{year}}\tab year of publication \cr
\bold{\emph{id}}\tab study ID \cr
\bold{\emph{treatment}}\tab treatment \cr
\bold{\emph{exac}}\tab one or more episodes of COPD exacerbation
  \cr
\bold{\emph{total}}\tab number of individuals in treatment arm
}
}
\source{
Baker WL, Baker EL, Coleman CI (2009):
Pharmacologic Treatments for Chronic Obstructive Pulmonary Disease:
A Mixed-Treatment Comparison Meta-analysis.
\emph{Pharmacotherapy: The Journal of Human Pharmacology and Drug
Therapy},
\bold{29}, 891--905
}
\description{
This dataset comes from a systematic review of randomized
controlled trials on pharmacologic treatments for chronic
obstructive pulmonary disease (COPD) (Baker et al., 2009).

The primary outcome, occurrence of one or more episodes of COPD
exacerbation, is binary (yes / no). For this outcome, five drug
treatments (fluticasone, budesonide, salmeterol, formoterol,
tiotropium) and two combinations (fluticasone + salmeterol,
budesonide + formoterol) were compared to placebo. The authors
considered the two combinations as separate treatments instead of
evaluating the individual components.
}
\note{
The dataset Baker2009 is identical to dataset
\code{\link[metadat]{dat.baker2009}} in R package \bold{metadat}.
}
\examples{
head(dat.baker2009)

# Example using pairwise(), netmeta() and netcomb():
# example(dat.baker2009, run.dontrun = TRUE)
}
\seealso{
\code{\link[metadat]{dat.baker2009}},
  \code{\link[meta]{pairwise}}, \code{\link[meta]{metabin}},
  \code{\link{netmeta}}, \code{\link{netcomb}}
}
\keyword{datasets}
