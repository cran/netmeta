#' Calculate rankogram
#'
#' @description
#' This function calculates the probabilities of each treatment being
#' at each possible rank and the SUCRAs (Surface Under the Cumulative
#' RAnking curve) in frequentist network meta-analysis.
#'
#' @param x An object of class \code{\link{netmeta}}.
#' @param nsim Number of simulations.
#' @param common A logical indicating to compute ranking probabilities
#'   and SUCRAs for the common effects model.
#' @param random A logical indicating to compute ranking probabilities
#'   and SUCRAs for the random effects model.
#' @param small.values A character string specifying whether small
#'   treatment effects indicate a beneficial (\code{"desirable"}) or
#'   harmful (\code{"undesirable"}) effect, can be abbreviated.
#' @param cumulative.rankprob A logical indicating whether cumulative
#'   ranking probabilities should be printed.
#' @param keep.samples A logical indicating whether to keep the generated
#'   samples.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' We derive a matrix showing the probability of each treatment being
#' at each possible rank. To this aim, we use resampling from a
#' multivariate normal distribution with estimated network effects as
#' means and corresponding estimated variance covariance matrix. We
#' then summarise them using the ranking metric SUCRAs (Surface Under
#' the Cumulative RAnking curve).
#'
#' @return
#' An object of class \code{rankogram} with corresponding \code{print}
#' and \code{plot} function. The object is a list containing the
#' following components:
#' \item{ranking.matrix.common}{Numeric matrix giving the probability
#'   of each treatment being at each possible rank for the common
#'   effects model.}
#' \item{ranking.common}{SUCRA values for the common effects model.}
#' \item{ranking.matrix.random}{Numeric matrix giving the probability
#'   of each treatment being at each possible rank for the random
#'   effects model.}
#' \item{ranking.random}{SUCRA values for the random effects model.}
#' \item{cumrank.matrix.common}{Numeric matrix giving the cumulative
#'   ranking probability of each treatment for the
#'   common effects model.}
#' \item{cumrank.matrix.random}{Numeric matrix giving the cumulative
#'   ranking probability of each treatment for the random effects
#'   model.}
#' \item{nsim, common, random}{As defined above},
#' \item{small.values, x}{As defined above},
#'
#' @author Theodoros Papakonstantinou \email{dev@@tpapak.com}, Guido
#'   Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @seealso \code{\link{netmeta}}, \code{\link{netrank}},
#'   \code{\link{plot.rankogram}}, \code{\link[metadat]{dat.woods2010}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP (2011):
#' Graphical methods and numerical summaries for presenting results
#' from multiple-treatment meta-analysis: an overview and tutorial.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{64}, 163--71
#'
#' @examples
#' pw1 <- pairwise(treatment, event = r, n = N, studlab = author,
#'   data = dat.woods2010, sm = "OR")
#' net1 <- netmeta(pw1, small.values = "desirable")
#'
#' set.seed(1909) # get reproducible results
#' ran1 <- rankogram(net1, nsim = 100)
#' ran1
#' print(ran1, cumulative.rankprob = TRUE)
#'
#' # Also print mean ranks
#' summary(ran1)
#'
#' plot(ran1)
#'
#' @rdname rankogram.netmeta
#' @method rankogram netmeta
#' @export

rankogram.netmeta <- function(x, nsim = gs("nsim"),
                              common = x$common, random = x$random,
                              small.values = x$small.values,
                              cumulative.rankprob = FALSE,
                              keep.samples = FALSE,
                              nchar.trts = x$nchar.trts,
                              warn.deprecated = gs("warn.deprecated"),
                              ...) {
  
  #
  #
  # (1) Check for netmeta object and upgrade object
  #
  #
  
  chkclass(x, "netmeta")
  chksuitable(x, "Rankograms", "netmetabin")
  x <- updateversion(x)
  
  
  #
  #
  # (2) Check other arguments
  #
  #
  
  chknumeric(nsim, min = 1, length = 1)
  #
  small.values <- setsv(small.values)
  #
  chklogical(cumulative.rankprob)
  chklogical(keep.samples)
  #
  if (is.null(nchar.trts))
    nchar.trts <- 666
  chknumeric(nchar.trts, length = 1)
  #
  # Check for deprecated arguments in '...'
  #
  args  <- list(...)
  chklogical(warn.deprecated)
  #
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "comb.fixed",
                      warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  #
  random <-
    deprecated(random, missing(random), args, "comb.random", warn.deprecated)
  chklogical(random)
  #
  # Additional checks for crossnma and multinma objects (and return results)
  #
  if (inherits(x, c("netmeta.crossnma", "netmeta.multinma"))) {
    if (!x$keep.samples)
      stop("Input for argument 'x' is a netmeta.", x$method,
           " object without samples; recreate ", x$method,
           " object with argument 'keep.samples = TRUE'.",
           call. = FALSE)
    #
    if (common != x$common)
      warning("Argument 'common = ", x$common, "' as netmeta.", x$method,
              " object is based on ", if (x$common) "common" else "random",
              " effects model.",
              call. = FALSE)
    #
    if (random != x$random)
      warning("Argument 'random = ", x$random, "' as netmeta.", x$method,
              " object is based on ", if (x$random) "random" else "common",
              " effects model.",
              call. = FALSE)
    #
    #
    if (!missing(nsim) & x$keep.samples)
      warning("Argument 'nsim' ignored for netmeta.", x$method, " object.",
              call. = FALSE)
    #
    common <- x$common
    random <- x$random
    #
    return(rankogram(x$samples$d,
                     pooled = if (common) "common" else "random",
                     small.values = small.values,
                     cumulative.rankprob = cumulative.rankprob,
                     keep.samples = keep.samples,
                     nchar.trts = nchar.trts))
  }
  
  
  #
  #
  # (3) Resampling to calculate ranking probabilities and SUCRAs
  #
  #
  
  sucras.common  <- ranking.matrix.common  <- cumrank.matrix.common <-
    meanranks.common <- medianranks.common <- NULL
  #
  sucras.random <- ranking.matrix.random <- cumrank.matrix.random <-
    meanranks.random <- medianranks.random <- NULL
  #
  if (common) {
    res.c <- ranksampling(x, nsim, "common", small.values, keep.samples)
    #
    sucras.common <- res.c$sucras
    ranking.matrix.common <- res.c$rankogram
    cumrank.matrix.common <- res.c$cumrank
    #
    meanranks.common <- res.c$meanranks
    medianranks.common <- res.c$medianranks
    #
    samples.common <- res.c$samples
  }
  #
  if (random) {
    res.r <- ranksampling(x, nsim, "random", small.values, keep.samples)
    #
    sucras.random <- res.r$sucras
    ranking.matrix.random <- res.r$rankogram
    cumrank.matrix.random <- res.r$cumrank
    #
    meanranks.random <- res.r$meanranks
    medianranks.random <- res.r$medianranks
    #
    samples.random <- res.r$samples
  }
  
  
  #
  #
  # (4) Create rankogram object
  #
  #
  
  res <- list(ranking.common = sucras.common,
              ranking.matrix.common = ranking.matrix.common,
              cumrank.matrix.common = cumrank.matrix.common,
              #
              meanranks.common = meanranks.common,
              medianranks.common = medianranks.common,
              #
              samples.common =
                if (common & keep.samples) samples.common else NULL,
              #
              ranking.random = sucras.random,
              ranking.matrix.random = ranking.matrix.random,
              cumrank.matrix.random = cumrank.matrix.random,
              #
              meanranks.random = meanranks.random,
              medianranks.random = medianranks.random,
              #
              samples.random =
                if (random & keep.samples) samples.random else NULL,
              #
              nsim = nsim,
              #
              common = common,
              random = random,
              small.values = small.values,
              cumulative.rankprob = cumulative.rankprob,
              #
              nchar.trts = nchar.trts,
              x = x,
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  #
  # Backward compatibility
  #
  res$fixed <- common
  #
  res$ranking.fixed <- sucras.common
  res$ranking.matrix.fixed <- ranking.matrix.common
  res$cumrank.matrix.fixed <- cumrank.matrix.common
  #
  class(res) <- "rankogram"
  
  res
}
