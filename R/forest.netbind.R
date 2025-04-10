#' Forest plot showing results of two or more network meta-analyses
#' 
#' @description
#' Forest plot to show network estimates of two or more network
#' meta-analyses.
#'
#' @aliases forest.netbind plot.netbind
#' 
#' @param x An object of class \code{netbind}.
#' @param pooled A character string indicating whether results for the
#'   common (\code{"common"}) or random effects model (\code{"random"})
#'   should be plotted. Can be abbreviated.
#' @param equal.size A logical indicating whether all squares should
#'   be of equal size. Otherwise, the square size is proportional to
#'   the precision of estimates.
#' @param leftcols A character vector specifying columns to be plotted
#'   on the left side of the forest plot (see Details).
#' @param leftlabs A character vector specifying labels for columns on
#'   left side of the forest plot.
#' @param rightcols A character vector specifying columns to be
#'   plotted on the right side of the forest plot (see Details).
#' @param rightlabs A character vector specifying labels for columns
#'   on right side of the forest plot.
#' @param subset.treatments A character vector specifying treatments
#'   to show in forest plot as comparators to the reference.
#' @param digits Minimal number of significant digits for treatment
#'   effects and confidence intervals, see \code{print.default}.
#' @param digits.prop Minimal number of significant digits for the
#'   direct evidence proportion.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in forest plots. If \code{backtransf = TRUE},
#'   results for \code{sm = "OR"} are presented as odds ratios rather
#'   than log odds ratios, for example.
#' @param lab.NA A character string to label missing values.
#' @param smlab A label printed at top of figure. By default, text
#'   indicating either common or random effects model is printed.
#' @param \dots Additional arguments for \code{\link[meta]{forest.meta}}
#'   function.
#' 
#' @details
#' A forest plot, also called confidence interval plot, is drawn in
#' the active graphics window.
#' 
#' The arguments \code{leftcols} and \code{rightcols} can be used to
#' specify columns which are plotted on the left and right side of the
#' forest plot, respectively. If argument \code{rightcols} is
#' \code{FALSE}, no columns will be plotted on the right side.
#' 
#' For more information see help page of \code{\link[meta]{forest.meta}}
#' function.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netbind}}, \code{\link{netcomb}},
#'   \code{\link[meta]{forest.meta}}
#' 
#' @keywords hplot
#' 
#' @examples
#' # Examples: example(netbind)
#' 
#' @method forest netbind
#' @export

forest.netbind <- function(x,
                           pooled = ifelse(x$x$random, "random", "common"),
                           ##
                           equal.size = gs("equal.size"),
                           ##
                           leftcols = "studlab",
                           leftlabs = "Treatment",
                           rightcols = c("effect", "ci"),
                           rightlabs = NULL,
                           ##
                           subset.treatments,
                           ##
                           digits = gs("digits.forest"),
                           digits.prop = max(gs("digits.pval") - 2, 2),
                           ##
                           backtransf = x$backtransf,
                           lab.NA = gs("lab.NA"),
                           smlab,
                           ...) {
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  chkclass(x, "netbind")
  x <- updateversion(x)
  ##
  pooled <- setchar(pooled, c("common", "random", "fixed"))
  pooled[pooled == "fixed"] <- "common"
  ##
  chklogical(equal.size)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.prop, min = 0, length = 1)
  chklogical(backtransf)
  ##
  chkchar(lab.NA)
  
  
  ##
  ##
  ## (2) Extract results for common and random effects model
  ##
  ##
  if (pooled == "common") {
    if (!missing(subset.treatments)) {
      subset.treatments <- setchar(subset.treatments, unique(x$common$treat))
      sel <- x$common$treat %in% subset.treatments
    }
    else
      sel <- x$common$treat != x$reference.group
    ##
    m <-
      suppressWarnings(metagen(x$common$TE, x$common$seTE,
                               studlab = x$common$name,
                               sm = x$sm,
                               common = FALSE, random = FALSE,
                               subgroup = x$common$treat,
                               print.subgroup.name = FALSE,
                               subset = sel,
                               method.tau = "DL", method.tau.ci = ""))
    ##
    m$studlab <- x$common$name[sel]
    m$TE <- x$common$TE[sel]
    m$seTE <- x$common$seTE[sel]
    m$lower <- x$common$lower[sel]
    m$upper <- x$common$upper[sel]
    m$statistic <- x$common$statistic[sel]
    m$pval <- x$common$pval[sel]
    m$zval <- x$common$statistic[sel]
    ##
    m$col.study <- x$common$col.study[sel]
    m$col.square <- x$common$col.square[sel]
    m$col.square.lines <- x$common$col.square.lines[sel]
    m$col.inside <- x$common$col.inside[sel]
    ##
    text.pooled <- "Common Effects Model"
    #
    if (!is.null(x$common$method))
      x$method <- unique(x$common$method[sel])
  }
  else {
    if (!missing(subset.treatments)) {
      subset.treatments <- setchar(subset.treatments, unique(x$random$treat))
      sel <- x$random$treat %in% subset.treatments
    }
    else
      sel <- x$random$treat != x$reference.group
    ##
    m <-
      suppressWarnings(metagen(x$random$TE, x$random$seTE,
                               studlab = x$random$name,
                               sm = x$sm,
                               common = FALSE, random = FALSE,
                               subgroup = x$random$treat,
                               print.subgroup.name = FALSE,
                               subset = sel,
                               method.tau = "DL", method.tau.ci = ""))
    ##
    m$studlab <- x$random$name[sel]
    m$TE <- x$random$TE[sel]
    m$seTE <- x$random$seTE[sel]
    m$lower <- x$random$lower[sel]
    m$upper <- x$random$upper[sel]
    m$statistic <- x$random$statistic[sel]
    m$pval <- x$random$pval[sel]
    m$zval <- x$random$statistic[sel]
    ##
    m$col.study <- x$random$col.study[sel]
    m$col.square <- x$random$col.square[sel]
    m$col.square.lines <- x$random$col.square.lines[sel]
    m$col.inside <- x$random$col.inside[sel]
    ##
    text.pooled <- "Random Effects Model"
    #
    if (!is.null(x$random$method))
      x$method <- unique(x$random$method[sel])
  }
  ##
  if (missing(smlab))
    if (x$baseline.reference)
      smlab <- paste0("Comparison: other vs '",
                      x$reference.group, "'\n(",
                      text.pooled,
                      ")")
    else
      smlab <- paste0("Comparison: '",
                      x$reference.group, "' vs other\n(",
                      text.pooled,
                      ")")
  #
  m$.text.details.methods <-
    textmeth(x, pooled == pooled, FALSE, FALSE,
             "", "", gs("digits.tau2"), gs("digits.tau"),
             FALSE, gs("text.I2"),
             big.mark = gs("big.mark"), forest = TRUE)
  
  
  ##
  ##
  ## (3) Forest plot
  ##
  ##
  forest(m,
         digits = digits,
         ##
         overall = FALSE, common = FALSE, random = FALSE,
         hetstat = FALSE, test.subgroup = FALSE,
         ##
         subgroup.hetstat = FALSE,
         prediction.subgroup = FALSE,
         calcwidth.subgroup = TRUE,
         ##
         leftcols = leftcols,
         leftlabs = leftlabs,
         rightcols = rightcols,
         rightlabs = rightlabs,
         ##
         lab.NA = lab.NA,
         smlab = smlab,
         backtransf = backtransf,
         ##
         col.study = m$col.study,
         col.square = m$col.square,
         col.square.lines = m$col.square.lines,
         col.inside = m$col.inside,
         col.inside.common = "black",
         col.inside.random = "black",
         ##
         weight.study = if (equal.size) "same" else pooled,
         ...)
  
  ret <- m
  ##
  ret$leftcols <- leftcols
  ret$rightcols <- rightcols
  ret$leftlabs <- leftlabs
  ret$rightlabs <- rightlabs
  ##
  invisible(ret)
}


#' @rdname forest.netbind
#' @method plot netbind
#' @export

plot.netbind <- function(x, ...)
  forest(x, ...)
