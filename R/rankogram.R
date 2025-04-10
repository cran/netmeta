#' Generic function for rankograms
#' 
#' @description
#' Generic function to calculate the probabilities of each treatment being
#' at each possible rank and the SUCRAs (Surface Under the Cumulative
#' RAnking curve) in network meta-analysis.
#' 
#' @param x An R object.
#' @param \dots Additional arguments.
#' 
#' @details
#' For more details, look at the following functions to generate
#' rankograms:
#' \itemize{
#' \item \code{\link{rankogram.netmeta}}
#' \item \code{\link{rankogram.default}}
#' }
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
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
#' # Examples:
#' # example(rankogram.netmeta)
#' # example(rankogram.default)
#'
#' @rdname rankogram
#' @export rankogram

rankogram <- function(x, ...)
  UseMethod("rankogram")
