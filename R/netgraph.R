#' Generic function for network graphs
#' 
#' @description
#' Generic function for network graphs
#' 
#' @param x An R object.
#' @param \dots Additional arguments.
#' 
#' @details
#' 
#' For more details, look at the following functions to generate
#' network graphs:
#' \itemize{
#' \item \code{\link{netgraph.netmeta}}
#' \item \code{\link{netgraph.netimpact}}
#' \item \code{\link{netgraph.netconnection}}
#' \item \code{\link{netgraph.netcomb}}
#' \item \code{\link{netgraph.discomb}}
#' }
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de
#' }
#' 
#' @seealso \code{\link[metadat]{dat.woods2010}}
#' 
#' @keywords hplot
#'
#' @examples
#' data(smokingcessation)
#' 
#' # Transform data from arm-based format to contrast-based format
#' #
#' pw1 <- pairwise(list(treat1, treat2, treat3),
#'   event = list(event1, event2, event3), n = list(n1, n2, n3),
#'   data = smokingcessation, sm = "OR")
#' 
#' # Conduct random effects network meta-analysis
#' #
#' net1 <- netmeta(pw1, common = FALSE)
#' 
#' # Network graph with default settings
#' #
#' netgraph(net1)
#' 
#' \donttest{
#' data(Senn2013)
#' 
#' # Network graphs with default settings
#' #
#' net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD", reference = "plac")
#' netgraph(net2)
#' #
#' pw3 <- pairwise(treatment, event = r, n = N,
#'   studlab = author, data = dat.woods2010, sm = "OR")
#' net3 <- netmeta(pw3)
#' netgraph(net3)
#' 
#' # Network graph with number of participants for each treatment arm
#' #
#' netgraph(net3, labels = paste0(trts, " (n=", n.trts, ")"))
#' }
#' 
#' @rdname netgraph
#' @export netgraph

netgraph <- function(x, ...)
  UseMethod("netgraph")
