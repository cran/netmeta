#' Contribution matrix in network meta-analysis
#' 
#' @description
#' This function generates the contribution of direct comparisons to
#' every network treatment comparison. The output is a matrix where
#' rows represent network treatment effects and columns represent the
#' contribution of direct treatment effects.
#' 
#' @aliases netcontrib print.netcontrib
#' 
#' @param x An object of class \code{netmeta} or \code{netcontrib}.
#' @param method A character string indicating which method is to
#'   calculate the contribution matrix. Either \code{"randomwalk"} or
#'   \code{"shortestpath"}, can be abbreviated.
#' @param hatmatrix.F1000 A logical indicating whether hat matrix
#'   given in F1000 article should be used for \code{method =
#'   "shortestpath"}.
#' @param fixed A logical indicating whether a contribution matrix
#'   should be printed for the fixed effect / common effect network
#'   meta-analysis.
#' @param random A logical indicating whether a contribution matrix
#'   should be printed for the random effects network meta-analysis.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names (see Details).
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments.
#' 
#' @details
#' In network meta-analysis (NMA), it is important to assess the
#' influence of limitations or other characteristics of individual
#' studies on the estimates obtained from the network. To this end,
#' the contribution matrix shows how much each direct treatment effect
#' contributes to each treatment effect estimate from network
#' meta-analysis.
#' 
#' We use ideas from graph theory to derive the proportion that is
#' contributed by each direct treatment effect. We start with the
#' 'projection' matrix in a two-step network meta-analysis model,
#' called the H matrix, which is analogous to the hat matrix in a
#' linear regression model. H entries are translated to proportion
#' contributions based on the observation that the rows of H can be
#' interpreted as flow networks.  A stream is defined as the
#' composition of a path and its associated flow (Papakonstantinou et
#' al., 2018).
#'
#' To account for multi-arm trials, we use the H matrix from a
#' two-step (aggregate) version of the graph theoretical NMA model
#' (Davies et al., 2021). This H matrix can be obtained from
#' \code{\link{hatmatrix}} with argument \code{method = "davies"}.
#' 
#' Two methods are implemented to estimate the streams and as a
#' result, the proportion contributions:
#' 
#' (1) If argument \code{method = "randomwalk"}, an analytical
#' random-walk (RW) approach is used (Davies et al., 2021). Here, the
#' "full" version of the aggregate H matrix (\code{\link{hatmatrix}}
#' with arguments \code{method = "davies"} and \code{type = "full"})
#' is used to define RW transition matrices. For each pair of
#' treatments (ij) in the network, the elements in the corresponding
#' row of H-full define a transition matrix from node i to node j. We
#' use the \bold{igraph} package to find every (directed) path from
#' node i to node j. The flow through each path is then equal to the
#' probability that a walker takes that path. This is simply the
#' product of the transition probabilities associated with each edge
#' along the path.
#' 
#' (2) If argument \code{method = "shortestpath"}, an iterative
#' algorithm is used (Papakonstantinou et al., 2018). Broadly
#' speaking, each iteration of the algorithm consists of the following
#' steps: (i) A path in the evidence flow network is selected. (ii)
#' The minimum flow through the edges making up the path is
#' identified. This is assigned as the flow associated with the
#' path. (iii) The flow of the path is subtracted from the values of
#' flow in the edges that make up that path. This means that the edge
#' corresponding to the minimum flow in that path is removed from the
#' graph. (iv) A new path is then selected from the remaining
#' graph. The process repeats until all the evidence flow in the edges
#' has been assigned to a path.
#' 
#' In the original F1000 paper (Papakonstantinou et al., 2018), the
#' hat matrix used did not account for correlations due to multi-arm
#' trials. For reproducibility the result of this version can be
#' obtained by specifying \code{hatmatrix.F1000 = TRUE} for
#' \code{method = "shortestpath"}. For other purposes, this method is
#' not recommended.
#' 
#' Once the streams have been identified (either by method (1) or
#' (2)), the proportion contribution of each direct comparison is
#' equal to the sum over the flow of evidence in each path containing
#' that edge divided by the number of edges that make up that path.
#'
#' By default, treatment names are not abbreviated in
#' printouts. However, in order to get more concise printouts,
#' argument \code{nchar.trts} can be used to define the minimum number
#' of characters for abbreviated treatment names (see
#' \code{\link{abbreviate}}, argument \code{minlength}). R function
#' \code{\link{treats}} is utilised internally to create abbreviated
#' treatment names.
#' 
#' @return
#' An object of class \code{netcontrib} with corresponding
#' \code{print} function. The object is a list containing the
#' following components:
#' \item{fixed}{Numeric matrix of percentage contributions of direct
#'   comparisons for each network comparison for the fixed effects
#'   model.}
#' \item{random}{Numeric matrix of percentage contributions of direct
#'   comparisons for each network comparison for the random effects
#'   model.}
#' \item{x}{As defined above.}
#' with the contribution matrices for fixed and random NMA. Each
#' matrix has the percentage contributions of each direct comparison
#' as columns for each network comparison, direct or indirect as rows.
#' 
#' @author Theodoros Papakonstantinou \email{dev@@tpapak.com}, Annabel
#'   Davies \email{annabel.davies@@manchester.ac.uk}
#' 
#' @seealso \code{\link{netmeta}}
#' 
#' @references
#' Davies AL, Papakonstantinou T, Nikolakopoulou A, Rücker G, Galla T
#' (2021):
#' Network meta-analysis and random walks.
#' Available from: http://arxiv.org/abs/2107.02886
#' 
#' Papakonstantinou, T., Nikolakopoulou, A., Rücker, G., Chaimani, A.,
#' Schwarzer, G., Egger, M., Salanti, G. (2018):
#' Estimating the contribution of studies in network meta-analysis:
#' paths, flows and streams.
#' \emph{F1000Research}
#' 
#' @keywords contribution
#' 
#' @examples
#' # Use the Woods dataset
#' #
#' data("Woods2010")
#' p1 <- pairwise(treatment, event = r, n = N,
#'                studlab = author, data = Woods2010, sm = "OR")
#' 
#' net1 <- netmeta(p1)
#' cm <- netcontrib(net1)
#' cm
#'
#' netcontrib(net1, method = "r")
#' 
#' @rdname netcontrib
#' @export netcontrib


netcontrib <- function(x,
                       method = "shortestpath",
                       hatmatrix.F1000 = FALSE,
                       fixed = x$fixed,
                       random = x$random,
                       nchar.trts = x$nchar.trts,
                       warn.deprecated = gs("warn.deprecated"),
                       ...) {
  
  ##
  ##
  ## (1) Check for netmeta object and upgrade object
  ##
  ##
  chkclass(x, "netmeta")
  x <- updateversion(x)
  ##
  is.installed.package("igraph")
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  method <- setchar(method, c("randomwalk", "shortestpath"))
  chklogical(hatmatrix.F1000)
  if (method == "randomwalk" & hatmatrix.F1000) {
    warning("Argument 'hatmatrix.F1000' ignored for random walk method.",
            call. = FALSE)
    hatmatrix.F1000 <- FALSE
  }
  chknumeric(nchar.trts, min = 1, length = 1)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  fixed <- deprecated(fixed, missing(fixed), args, "comb.fixed",
                      warn.deprecated)
  chklogical(fixed)
  ##
  random <- deprecated(random, missing(random), args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  
  
  ##
  ##
  ## (3) Create netcontrib object
  ##
  ##
  x$fixed <- fixed
  x$random <- random
  ##
  res <- list(fixed =
                contribution.matrix(x, method, "fixed", hatmatrix.F1000),
              random =
                contribution.matrix(x, method, "random", hatmatrix.F1000),
              method = method,
              hatmatrix.F1000 = hatmatrix.F1000,
              nchar.trts = nchar.trts,
              x = x,
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- "netcontrib"
  ##
  res
}





#' @rdname netcontrib
#' 
#' @method print netcontrib
#' 
#' @export


print.netcontrib <- function(x,
                             fixed = x$x$fixed,
                             random = x$x$random,
                             digits = 4,
                             nchar.trts = x$nchar.trts,
                             legend = TRUE,
                             warn.deprecated = gs("warn.deprecated"),
                             ...) {
  
  ##
  ##
  ## (1) Check for netcontrib object and upgrade object
  ##
  ##
  chkclass(x, "netcontrib")
  x <- updateversion(x)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  chknumeric(nchar.trts, length = 1)
  chklogical(legend)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  fixed <- deprecated(fixed, missing(fixed), args, "comb.fixed",
                      warn.deprecated)
  chklogical(fixed)
  ##
  random <- deprecated(random, missing(random), args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  
  
  ##
  ##
  ## (3) Print results for contribution matrix
  ##
  ##
  matitle(x$x)
  ##
  cat(paste0("Contribution matrix (",
             if (is.null(x$method) | x$method == "shortestpath")
               "Papakonstantinou et al., 2018, F1000Research"
             else
               "Davies et al., 2021",
             ")"))
  if ((is.null(x$method) | x$method == "shortestpath") & x$hatmatrix.F1000)
    cat(paste(",\nhat matrix does not take correlation of",
              "multi-arm studies into account"))
  ##
  cat("\n\n")
  
  
  ##
  trts <- x$x$trts
  trts.abbr <- treats(trts, nchar.trts)
  ##
  if (fixed) {
    rownames(x$fixed) <- comps(x$fixed, trts, x$x$sep.trts, nchar.trts)
    colnames(x$fixed) <- comps(x$fixed, trts, x$x$sep.trts, nchar.trts,
                               row = FALSE)
    ##
    cat("Fixed effects model:\n\n")
    prmatrix(round(x$fixed, digits))
    if (random)
      cat("\n")
  }
  if (random) {
    rownames(x$random) <- comps(x$random, trts, x$x$sep.trts, nchar.trts)
    colnames(x$random) <- comps(x$random, trts, x$x$sep.trts, nchar.trts,
                                row = FALSE)
    ##
    cat("Random effects model:\n\n")
    prmatrix(round(x$random, digits))
  }
  ##
  ## Add legend
  ##
  if (legend & (fixed | random)) {
    diff.trts <- trts != trts.abbr
    if (any(diff.trts)) {
      tmat <- data.frame(trts, trts.abbr)
      names(tmat) <- c("Abbreviation", "Treatment name")
      tmat <- tmat[diff.trts, ]
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      cat("\nLegend:\n")
      prmatrix(tmat, quote = FALSE, right = TRUE,
               rowlab = rep("", length(trts.abbr)))
    }
  }
  ##
  invisible(NULL)
}