#' Forest plot for direct and indirect evidence
#' 
#' @description
#' Forest plot to show direct and indirect evidence in network
#' meta-analysis.  Furthermore, estimates from network meta-analysis
#' as well as prediction intervals can be printed.
#'
#' @aliases forest.netsplit plot.netsplit
#' 
#' @param x An object of class \code{netsplit}.
#' @param pooled A character string indicating whether results for the
#'   common (\code{"common"}) or random effects model (\code{"random"})
#'   should be plotted. Can be abbreviated.
#' @param show A character string indicating which comparisons should
#'   be printed (see Details).
#' @param overall A logical indicating whether network meta-analysis
#'   estimates should be printed.
#' @param direct A logical indicating whether direct estimates should
#'   be printed.
#' @param indirect A logical indicating whether indirect estimates
#'   should be printed.
#' @param prediction A logical indicating whether prediction intervals
#'   should be printed.
#' @param only.reference A logical indicating whether only comparisons
#'   with the reference group should be printed.
#' @param sortvar An optional vector used to sort comparisons (must be
#'   of same length as the total number of comparisons).
#' @param subset An optional logical vector specifying a subset of
#'   comparisons to print (must be of same length as the total number of
#'   comparisons) .
#' @param subgroup A character string indicating which layout should
#'   be used in forest plot: subgroups by comparisons
#'   (\code{"comparison"}) or subgroups by estimates
#'   (\code{"estimate"}). Can be abbreviated.
#' @param text.overall A character string used in the plot to label
#'   the network estimates.
#' @param text.direct A character string used in the plot to label the
#'   direct estimates.
#' @param text.indirect A character string used in the plot to label
#'   the indirect estimates.
#' @param text.predict A character string used in the plot to label
#'   the prediction interval.
#' @param type.overall A character string specifying how to plot
#'   treatment effects and confidence intervals for the overall
#'   network evidence.
#' @param type.direct A character string specifying how to plot
#'   treatment effects and confidence intervals for the direct
#'   evidence.
#' @param type.indirect A character string specifying how to plot
#'   treatment effects and confidence intervals for the indirect
#'   evidence.
#' @param col.square The colour for squares.
#' @param col.square.lines The colour for the outer lines of squares.
#' @param col.inside The colour for results and confidence limits if
#'   confidence limits are completely within squares squares.
#' @param col.diamond The colour of diamonds.
#' @param col.diamond.lines The colour of the outer lines of diamonds.
#' @param col.predict Background colour of prediction intervals.
#' @param col.predict.lines Colour of outer lines of prediction
#'   intervals.
#' @param col.subgroup The colour to print information on subgroups, i.e.,
#'   pairwise comparisons.
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
#' If direct estimates are included in the forest plot (\code{direct =
#' TRUE}, default), the following columns will be printed on the left
#' side of the forest plot: the comparisons (column \code{"studlab"}
#' in \code{\link[meta]{forest.meta}}), number of pairwise comparisons
#' (\code{"k"}), direct evidence proportion (\code{"k"}), and
#' I\eqn{^2} from pairwise comparison (\code{"I2"}).
#' 
#' If direct estimates are not included in the forest plot
#' (\code{direct = FALSE}), only the comparisons (\code{"studlab"})
#' are printed on the left side of the forest plot.
#' 
#' For more information see help page of \code{\link[meta]{forest.meta}}
#' function.
#' 
#' Argument \code{show} determines which comparisons are printed:
#' \tabular{ll}{
#' \dQuote{all} \tab All comparisons \cr
#' \dQuote{both} \tab Only comparisons contributing both direct and
#'   indirect evidence \cr
#' \dQuote{with.direct} \tab Comparisons providing direct evidence \cr
#' \dQuote{direct.only} \tab Comparisons providing only direct
#'   evidence \cr
#' \dQuote{indirect.only} \tab Comparisons providing only indirect
#'   evidence
#' }
#'
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link[meta]{forest.meta}}
#' 
#' @keywords hplot
#' 
#' @examples
#' # Examples: example(netsplit)
#' 
#' @method forest netsplit
#' @export

forest.netsplit <- function(x,
                            pooled = ifelse(x$x$random, "random", "common"),
                            show = x$show,
                            #
                            subgroup = "comparison",
                            #
                            overall = x$overall,
                            direct = x$direct,
                            indirect = x$indirect,
                            prediction = x$prediction,
                            #
                            #test = x$test,
                            only.reference = x$only.reference,
                            #
                            sortvar = NULL,
                            subset = NULL,
                            ##
                            text.overall = "Network estimate",
                            text.direct = "Direct estimate",
                            text.indirect = "Indirect estimate",
                            text.predict = "Prediction interval",
                            ##
                            type.overall,
                            type.direct,
                            type.indirect,
                            ##
                            col.square = "gray",
                            col.square.lines = col.square,
                            col.inside = "white",
                            col.diamond = "gray",
                            col.diamond.lines = "black",
                            col.predict = "red",
                            col.predict.lines = "black",
                            col.subgroup = "black",
                            ##
                            equal.size = TRUE,
                            ##
                            leftcols,
                            leftlabs,
                            rightcols = c("effect", "ci"),
                            rightlabs = NULL,
                            ##
                            digits = gs("digits.forest"),
                            digits.prop = max(gs("digits.pval") - 2, 2),
                            ##
                            backtransf = x$backtransf,
                            lab.NA = "",
                            smlab,
                            ...) {


  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  chkclass(x, "netsplit")
  x <- updateversion(x)
  ##
  pooled <- setchar(pooled, c("common", "random", "fixed"))
  pooled[pooled == "fixed"] <- "common"
  ##
  subgroup <- setchar(subgroup, c("comparison", "estimate"))
  ##
  chklogical(overall)
  chklogical(direct)
  chklogical(indirect)
  chklogical(prediction)
  ##
  missing.only.reference <- missing(only.reference)
  if (!missing.only.reference)
    chklogical(only.reference)
  #
  # Catch sortvar and subset from data:
  #
  mc <- match.call()
  sfsp <- sys.frame(sys.parent())
  #
  error <- try(sortvar.x <- catch("sortvar", mc, x, sfsp), silent = TRUE)
  if (!any(class(error) == "try-error"))
    sortvar <- sortvar.x
  ##
  if (!is.null(sortvar)) {
    if (length(dim(sortvar)) == 2) {
      if (dim(sortvar)[1] != length(x$comparison))
        stop("Argument 'sortvar' must be of length ",
             length(x$comparison), ".",
             call. = FALSE)
      ##
      ## Set proportions to 0 or 1
      ##
      if (is.numeric(sortvar)) {
        sortvar[is_zero(abs(sortvar), n = 1000)] <- 0
        sortvar[is_zero(1 - abs(sortvar), n = 1000)] <-
          1 * sign(sortvar)[is_zero(1 - abs(sortvar), n = 1000)]
      }
      sortvar <- order(do.call(order, as.list(as.data.frame(sortvar))))
    }
    else
      chklength(sortvar, length(x$comparison),
                text = paste0("Argument 'sortvar' must be of length ",
                              length(x$comparison), "."))
    ##
    if (!is.numeric(sortvar))
      sortvar <- setchar(sortvar, x$comparison)
  }
  #
  error <- try(subset.x <- catch("subset", mc, x, sfsp), silent = TRUE)
  if (!any(class(error) == "try-error"))
    subset <- subset.x
  ##
  if (!is.null(subset)) {
    chklength(subset, length(x$comparison),
              text = paste0("Argument 'subset' must be of length ",
                            length(x$comparison), "."))
    if (!is.logical(subset))
      stop("Argument 'subset' must be a logical vector.")
  }
  ##
  chkchar(text.overall)
  chkchar(text.direct)
  chkchar(text.indirect)
  chkchar(text.predict)
  ##
  missing.type.overall <- missing(type.overall)
  if (missing.type.overall)
    type.overall <- "diamond"
  else
    type.overall <- setchar(type.overall, c("diamond", "square"))
  ##
  if (missing(type.direct))
    type.direct <- "square"
  else
    type.direct <- setchar(type.direct, c("diamond", "square"))
  if (missing(type.indirect))
    type.indirect <- "square"
  else
    type.indirect <- setchar(type.indirect, c("diamond", "square"))
  ##
  chkcolor(col.square)
  chkcolor(col.square.lines)
  chkcolor(col.inside)
  chkcolor(col.diamond)
  chkcolor(col.diamond.lines)
  chkcolor(col.predict)
  chkcolor(col.predict.lines)
  chkcolor(col.subgroup)
  ##
  chklogical(equal.size)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.prop, min = 0, length = 1)
  chklogical(backtransf)
  ##
  chkchar(lab.NA)
  ##
  if (pooled == "common") {
    if (!(missing(prediction)) & prediction)
      warning("Prediction intervals not shown for estimates ",
              "from common effects model.")
    prediction <- FALSE
  }
  ##
  if (!any(c(overall, direct, indirect)))
    stop("At least, one of the following estimates ",
         "must be included in forest plot:\n",
         "- network estimates (argument 'overall')\n",
         "- direct estimates (argument 'direct')\n",
         "- indirect estimates (argument 'indirect')")
  ##
  missing.leftcols <- missing(leftcols)
  if (missing.leftcols)
    if (direct)
      leftcols <- c("studlab", "k", "prop", "I2")
    else
      leftcols <- "studlab"
  ##
  missing.leftlabs <- missing(leftlabs)
  if (missing.leftlabs) {
    leftlabs <- rep(NA, length(leftcols))
    leftlabs[leftcols == "studlab"] <- "Comparison"
    leftlabs[leftcols == "k"] <- "Number of\nStudies"
    leftlabs[leftcols == "n"] <- "Number of\nParticipants"
    leftlabs[leftcols == "prop"] <- "Direct\nEvidence"
    leftlabs[leftcols == "I2"] <- "I2"
    leftlabs[leftcols == "tau2"] <- "tau2"
    leftlabs[leftcols == "tau"] <- "tau"
    leftlabs[leftcols == "Q"] <- "Q"
  }
  ##
  n.subgroup <- direct + indirect + overall + prediction
  missing.smlab <- missing(smlab)
  ##
  if (n.subgroup == 1 & overall & missing.type.overall)
    type.overall <- "square"
  ##
  if (missing(text.predict))
    if (!(length(x$level.predict) == 0) &&
        x$level.ma != x$level.predict)
      text.predict <- paste0(text.predict, " (",
                             round(x$level.predict * 100), "%-PI)")
  ##
  if (overall & n.subgroup > 1) {
    if (text.overall == text.predict)
      stop("Text must be different for arguments 'text.overall' and ",
           "'text.predict'.")
    if (text.overall == text.direct)
      stop("Text must be different for arguments 'text.overall' and ",
           "'text.direct'.")
    if (text.overall == text.indirect)
      stop("Text must be different for arguments 'text.overall' and ",
           "'text.indirect'.")
  }
  ##
  if (prediction & n.subgroup > 1) {
    if (text.predict == text.overall)
      stop("Text must be different for arguments 'text.predict' and ",
           "'text.overall'.")
    if (text.predict == text.direct)
      stop("Text must be different for arguments 'text.predict' and ",
           "'text.direct'.")
    if (text.predict == text.indirect)
      stop("Text must be different for arguments 'text.predict' and ",
           "'text.indirect'.")
  }
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  ## Check whether first argument is a list. In this case only use
  ## this list as input.
  if (length(args) > 0 && is.list(args[[1]]))
    args <- args[[1]]
  ##
  additional.arguments <- names(args)
  ##
  if (length(additional.arguments) > 0) {
    if (!is.na(charmatch("showa", additional.arguments)))
      if (!missing(show))
        warning("Deprecated argument 'showall' ignored as ",
                "argument 'show' is also provided.")
      else {
        warning("Deprecated argument 'showall' has been replaced by ",
                "argument 'show'.")
        show <- args[[charmatch("showa", additional.arguments)]]
        if (show)
          show <- "all"
        else
          show <- "both"
      }
  }
  ##
  show <- setchar(show, c("all", "both", "with.direct",
                          "direct.only", "indirect.only",
                          "reference.only"))
  ##
  if (show == "reference.only") {
    warning("Argument 'show = \"reference.only\" replaced with ",
            "'only.reference = TRUE'.",
            call. = FALSE)
    show <- "both"
    if (missing.only.reference)
      only.reference <- TRUE
  }
  
  
  ##
  ##
  ## (2) Extract results for common and random effects model
  ##
  ##
  if (pooled == "common") {
    dat.direct <- x$direct.common
    ##
    dat.indirect <- x$indirect.common
    dat.indirect$Q <- dat.indirect$tau2 <-
      dat.indirect$tau <- dat.indirect$I2 <- NA
    ##
    dat.overall <- x$common
    dat.overall$Q <- dat.overall$tau2 <-
      dat.overall$tau <- dat.overall$I2 <- NA
    ##
    dat.direct$prop <- formatPT(x$prop.common, digits = digits.prop)
    dat.indirect$prop <- NA
    dat.overall$prop <- NA
    ##
    if (missing.smlab)
      smlab <- "Common Effects Model"
  }
  else {
    dat.direct <- x$direct.random
    ##
    dat.indirect <- x$indirect.random
    dat.indirect$Q <- dat.indirect$tau2 <-
      dat.indirect$tau <- dat.indirect$I2 <- NA
    ##
    dat.overall <- x$random
    dat.overall$Q <- dat.overall$tau2 <-
      dat.overall$tau <- dat.overall$I2 <- NA
    ##
    dat.direct$prop <- formatPT(x$prop.random, digits = digits.prop)
    dat.indirect$prop <- NA
    dat.overall$prop <- NA
    ##
    if (missing.smlab)
      smlab <- "Random Effects Model"
  }
  ##
  if (missing.smlab & n.subgroup == 1)
    smlab <- paste0(if (direct)
      paste0(text.direct, "\n"),
      if (indirect)
        paste0(text.indirect, "\n"),
      if (overall)
        paste0(text.overall, "\n"),
      "(",
      tolower(smlab),
      ")")
  ##
  dat.predict <- x$predict
  ##
  if ( (subgroup == "estimate") | (prediction & !overall) )
    dat.predict$TE <- dat.overall$TE
  else
    dat.predict$TE <- NA
  ##
  dat.predict$seTE <- dat.predict$statistic <- dat.predict$p <-
    dat.predict$Q <- NA
  ##
  dat.predict <- dat.predict[, c("comparison", "TE", "seTE",
                                 "lower", "upper", "statistic", "p", "Q")]
  dat.predict$tau2 <- NA
  dat.predict$tau <- NA
  dat.predict$I2 <- NA
  dat.predict$prop <- NA
  ##
  dat.direct$comps <- dat.indirect$comps <-
    dat.overall$comps <- dat.predict$comps <- x$comparison
  ##
  dat.direct$k <- x$k
  dat.indirect$k <- dat.overall$k <- dat.predict$k <- NA
  ##
  dat.indirect$n <- dat.overall$n <- dat.predict$n <- NA
  ##
  dat.direct$evidence   <- text.direct
  dat.indirect$evidence <- text.indirect
  dat.overall$evidence  <- text.overall
  dat.predict$evidence  <- text.predict
  ##
  dat.direct$type.study <- type.direct
  dat.indirect$type.study <- type.indirect
  dat.overall$type.study <- type.overall
  dat.predict$type.study <- "predict"
  ##
  dat.direct$col.estimate <- if (type.direct == "square")
                               col.square
                             else
                               col.diamond
  dat.indirect$col.estimate <- if (type.indirect == "square")
                                 col.square
                               else
                                 col.diamond
  dat.overall$col.estimate <- if (type.overall == "square")
                                col.square
                              else
                                col.diamond
  ##
  dat.direct$col.lines <- if (type.direct == "square")
                            col.square.lines
                          else
                            col.diamond.lines
  dat.indirect$col.lines <- if (type.indirect == "square")
                              col.square.lines
                            else
                              col.diamond
  dat.overall$col.lines <- if (type.overall == "square")
                             col.square.lines
                           else
                             col.diamond.lines
  ##
  dat.predict$col.estimate <- col.predict
  dat.predict$col.lines <- col.predict.lines
  ##
  ## col.square.lines = col.square,
  ## col.inside = "white",
  ## col.diamond.lines = "black",
  ## col.predict.lines = "black",
  
  
  ##
  ##
  ## (3) Select treatment comparisons to show in forest plot
  ##
  ##
  if (show == "all")
    sel <- rep_len(TRUE, length(x$direct.common$TE))
  else if (show == "with.direct")
    sel <- (!is.na(x$direct.common$TE) & !is.na(x$direct.random$TE))
  else if (show == "both")
    sel <- (!is.na(x$direct.common$TE)  & !is.na(x$indirect.common$TE) &
            !is.na(x$direct.random$TE) & !is.na(x$indirect.random$TE))
  else if (show == "direct.only")
    sel <- (!is.na(x$direct.common$TE)  & is.na(x$indirect.common$TE) &
            !is.na(x$direct.random$TE) & is.na(x$indirect.random$TE))
  else if (show == "indirect.only")
    sel <- (is.na(x$direct.common$TE)  & !is.na(x$indirect.common$TE) &
            is.na(x$direct.random$TE) & !is.na(x$indirect.random$TE))
  ##
  if (only.reference) {
    if (x$reference.group == "") {
      warning("First treatment used as reference as argument ",
              "'reference.group' was unspecified in netsplit().",
              call. = FALSE)
      x$reference.group <-
        compsplit(x$comparison, x$sep.trts)[[1]][1]
    }
    ##
    sel.ref <-
      apply(!is.na(sapply(compsplit(x$comparison, x$sep.trts),
                          match, x$reference.group)), 2, sum) >= 1
    ##
    sel <- sel & sel.ref
  }
  ##
  dat.direct <- dat.direct[sel, ]
  dat.indirect <- dat.indirect[sel, ]
  dat.overall <- dat.overall[sel, ]
  dat.predict <- dat.predict[sel, ]
  #
  if (!is.null(sortvar))
    sortvar <- sortvar[sel]
  #
  if (!is.null(subset))
    subset <- subset[sel]
  #
  if (!is.null(sortvar)) {
    o <- order(sortvar)
    ##
    dat.direct <- dat.direct[o, ]
    dat.indirect <- dat.indirect[o, ]
    dat.overall <- dat.overall[o, ]
    dat.predict <- dat.predict[o, ]
    #
    if (!is.null(subset))
      subset <- subset[o]
  }
  #
  if (!is.null(subset)) {
    dat.direct <- dat.direct[subset, ]
    dat.indirect <- subset(dat.indirect, subset)
    dat.overall <- subset(dat.overall, subset)
    dat.predict <- subset(dat.predict, subset)
  }
  ##
  if (direct & all(is.na(dat.direct$I2)) & missing.leftcols) {
    leftlabs <- leftlabs[leftcols != "I2"]
    leftcols <- leftcols[leftcols != "I2"]
  }
  
  
  ##
  ##
  ## (4) Forest plot
  ##
  ##
  if (subgroup == "comparison") {
    dat <- rbind(
      if (direct) dat.direct,
      if (indirect) dat.indirect,
      if (overall) dat.overall,
      if (prediction) dat.predict
    )
    ##
    if (nrow(dat) == 0) {
      warning("No comparison(s) selected. Consider using argument ",
              "'show = \"all\"'.")
      return(invisible(NULL))
    }
    ##
    if (n.subgroup > 1)
      m <-
        suppressWarnings(metagen(dat$TE, dat$seTE,
                                 studlab = dat$evidence, data = dat,
                                 sm = x$sm,
                                 common = FALSE, random = FALSE,
                                 method.tau = "DL", method.tau.ci = "",
                                 subgroup = dat$comps,
                                 print.subgroup.name = FALSE))
    else
      m <-
        suppressWarnings(metagen(dat$TE, dat$seTE,
                                 studlab = dat$comps, data = dat, sm = x$sm,
                                 common = FALSE, random = FALSE,
                                 method.tau = "DL", method.tau.ci = ""))
    ##
    if (overall) {
      m$w.common[m$studlab == text.overall] <- max(m$w.common, na.rm = TRUE)
      m$w.random[m$studlab == text.overall] <- max(m$w.random, na.rm = TRUE)
    }
    ##
    if (prediction) {
      m$lower[m$studlab == text.predict] <- dat.predict$lower
      m$upper[m$studlab == text.predict] <- dat.predict$upper
      ##
      m$w.common[m$studlab == text.predict] <- max(m$w.common, na.rm = TRUE)
      m$w.random[m$studlab == text.predict] <- max(m$w.random, na.rm = TRUE)
    }
    ##
    forest(m,
           digits = digits,
           common = FALSE, random = FALSE,
           hetstat = FALSE, test.subgroup = FALSE,
           leftcols = leftcols,
           leftlabs = leftlabs,
           rightcols = rightcols,
           rightlabs = rightlabs,
           lab.NA = lab.NA,
           smlab = smlab,
           backtransf = backtransf,
           type.study = dat$type.study,
           col.square = dat$col.estimate,
           col.square.lines = dat$col.lines,
           weight.study = if (equal.size) "same" else "common",
           ...)
  }
  else {
    dat <- rbind(
      if (direct) dat.direct,
      if (indirect) dat.indirect,
      if (overall) dat.overall,
      if (prediction) dat.predict
    )
    ##
    if (nrow(dat) == 0) {
      warning("No comparison(s) selected. Consider using argument ",
              "'show = \"all\"'.")
      return(invisible(NULL))
    }
    ##
    if (n.subgroup > 1)
      m <-
        suppressWarnings(metagen(dat$TE, dat$seTE,
                                 studlab = dat$comps, data = dat,
                                 sm = x$sm,
                                 common = FALSE, random = FALSE,
                                 method.tau = "DL", method.tau.ci = "",
                                 subgroup = dat$evidence,
                                 print.subgroup.name = FALSE))
    else
      m <-
        suppressWarnings(metagen(dat$TE, dat$seTE,
                                 studlab = dat$comps, data = dat, sm = x$sm,
                                 common = FALSE, random = FALSE,
                                 method.tau = "DL", method.tau.ci = ""))
    ##
    if (overall) {
      m$w.common[m$subgroup == text.overall] <- max(m$w.common, na.rm = TRUE)
      m$w.random[m$subgroup == text.overall] <- max(m$w.random, na.rm = TRUE)
    }
    ##
    if (prediction) {
      m$lower[m$subgroup == text.predict] <- dat.predict$lower
      m$upper[m$subgroup == text.predict] <- dat.predict$upper
      ##
      m$w.common[m$subgroup == text.predict] <- max(m$w.common, na.rm = TRUE)
      m$w.random[m$subgroup == text.predict] <- max(m$w.random, na.rm = TRUE)
    }
    ##
    forest(m,
           digits = digits,
           ##
           overall = FALSE, common = FALSE, random = FALSE,
           hetstat = FALSE, test.subgroup = FALSE,
           ##
           subgroup.hetstat = FALSE,
           prediction.subgroup = FALSE,
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
           type.study = dat$type.study,
           col.square = dat$col.estimate,
           col.square.lines = dat$col.lines,
           ##
           weight.study = if (equal.size) "same" else "common",
           ...)
  }


  res <- list(dat = dat,
              leftcols = leftcols,
              leftlabs = leftlabs,
              rightcols = rightcols,
              rightlabs = rightlabs,
              lab.NA = lab.NA,
              backtransf = backtransf,
              smlab = smlab,
              weight.study = if (equal.size) "same" else "common",
              ##
              args = list(
                pooled = pooled,
                show = show,
                subgroup = subgroup,
                overall = overall,
                direct = direct,
                indirect = indirect,
                prediction = prediction,
                only.reference = only.reference,
                sortvar = sortvar,
                subset = subset,
                text.overall = text.overall,
                text.direct = text.direct,
                text.indirect = text.indirect,
                text.predict = text.predict,
                type.overall = type.overall,
                type.direct = type.direct,
                type.indirect = type.indirect,
                col.square = col.square,
                col.square.lines = col.square.lines,
                col.inside = col.inside,
                col.diamond = col.diamond,
                col.diamond.lines = col.diamond.lines,
                col.predict = col.predict,
                col.predict.lines = col.predict.lines,
                col.subgroup = col.subgroup,
                equal.size = equal.size,
                digits = digits,
                digits.prop = digits.prop)
              )

  
  invisible(res)
}


#' @rdname forest.netsplit
#' @method plot netsplit
#' @export

plot.netsplit <- function(x, ...)
  forest(x, ...)
