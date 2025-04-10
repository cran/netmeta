#' Print method for objects of class netcomb
#' 
#' @description
#' Print method for objects of class \code{netcomb}.
#' 
#' @param x An object of class \code{netcomb} or
#'   \code{summary.netcomb}.
#' @param common A logical indicating whether results for the common
#'   effects model should be printed.
#' @param random A logical indicating whether results for the random
#'   effects model should be printed.
#' @param overall.hetstat A logical indicating whether to print heterogeneity
#'   measures.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots. If
#'   \code{backtransf = TRUE}, results for \code{sm = "OR"} are
#'   presented as odds ratios rather than log odds ratios, for
#'   example.
#' @param nchar.comps A numeric defining the minimum number of
#'   characters used to create unique names for components.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.stat Minimal number of significant digits for z- or
#'   t-value, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall treatment effect, see \code{print.default}.
#' @param digits.pval.Q Minimal number of significant digits for
#'   p-value of heterogeneity tests, see \code{print.default}.
#' @param digits.Q Minimal number of significant digits for
#'   heterogeneity statistics, see \code{print.default}.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance, see \code{print.default}.
#' @param digits.tau Minimal number of significant digits for square
#'   root of between-study variance, see \code{print.default}.
#' @param digits.I2 Minimal number of significant digits for I-squared
#'   statistic, see \code{print.default}.
#' @param big.mark A character used as thousands separator.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   overall effect should be printed according to JAMA reporting
#'   standards.
#' @param print.tau2 A logical specifying whether between-study
#'   variance \eqn{\tau^2} should be printed.
#' @param print.tau A logical specifying whether \eqn{\tau}, the
#'   square root of the between-study variance \eqn{\tau^2}, should be
#'   printed.
#' @param print.Q A logical value indicating whether to print the
#'   results of the test of heterogeneity.
#' @param print.I2 A logical specifying whether heterogeneity
#'   statistic I\eqn{^2} should be printed.
#' @param print.I2.ci A logical specifying whether confidence interval for
#'   heterogeneity statistic I\eqn{^2} should be printed.
#' @param text.tau2 Text printed to identify between-study variance
#'   \eqn{\tau^2}.
#' @param text.tau Text printed to identify \eqn{\tau}, the square
#'   root of the between-study variance \eqn{\tau^2}.
#' @param text.I2 Text printed to identify heterogeneity statistic
#'   I\eqn{^2}.
#' @param details.methods A logical specifying whether details on statistical
#'   methods should be printed.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
#'
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netcomb}}, \code{\link{discomb}}
#' 
#' @keywords print
#' 
#' @examples
#' # Examples: example(netcomb)
#' 
#' @method print netcomb
#' @export

print.netcomb <- function(x,
                          common = x$common,
                          random = x$random,
                          overall.hetstat = x$overall.hetstat,
                          backtransf = x$backtransf,
                          nchar.comps = x$nchar.comps,
                          ##
                          digits = gs("digits"),
                          digits.stat = gs("digits.stat"),
                          digits.pval = gs("digits.pval"),
                          digits.pval.Q = max(gs("digits.pval.Q"), 2),
                          digits.Q = gs("digits.Q"),
                          digits.tau2 = gs("digits.tau2"),
                          digits.tau = gs("digits.tau"),
                          digits.I2 = gs("digits.I2"),
                          #
                          big.mark = gs("big.mark"),
                          scientific.pval = gs("scientific.pval"),
                          zero.pval = gs("zero.pval"),
                          JAMA.pval = gs("JAMA.pval"),
                          #
                          print.tau2 = gs("print.tau2"),
                          print.tau = gs("print.tau"),
                          print.Q = gs("print.Q"),
                          print.I2 = gs("print.I2"),
                          print.I2.ci = gs("print.I2.ci"),
                          #
                          text.tau2 = gs("text.tau2"),
                          text.tau = gs("text.tau"),
                          text.I2 = gs("text.I2"),
                          #
                          details.methods = gs("details"),
                          legend = gs("legend"),
                          ##
                          warn.deprecated = gs("warn.deprecated"),
                          ##
                          ...) {
  
  
  ##
  ##
  ## (1) Check for netcomb object and upgrade object
  ##
  ##
  chkclass(x, "netcomb")
  x <- updateversion(x)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  chklogical(overall.hetstat)
  chklogical(backtransf)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  chknumeric(digits.pval.Q, min = 1, length = 1)
  chknumeric(digits.Q, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  #
  chkchar(big.mark, length = 1)
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  #
  chklogical(print.tau2)
  chklogical(print.tau)
  chklogical(print.Q)
  chklogical(print.I2)
  chklogical(print.I2.ci)
  #
  chkchar(text.tau2)
  chkchar(text.tau)
  chkchar(text.I2)
  #
  chklogical(details.methods)
  chklogical(legend)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "comb.fixed",
                      warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  ##
  random <- deprecated(random, missing(random), args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  ##
  nchar.comps <-
    deprecated(nchar.comps, missing(nchar.comps), args, "nchar.trts")
  
  
  ##
  ##
  ## (3) Print results
  ##
  ##
  I2 <- round(100 * x$I2, digits.I2)
  lower.I2 <- round(100 * x$lower.I2, digits.I2)
  upper.I2 <- round(100 * x$upper.I2, digits.I2)
  ##
  if (common | random) {
    cat("Number of studies: k = ", x$k, "\n", sep = "")
    cat("Number of pairwise comparisons: m = ", x$m, "\n", sep = "")
    cat("Number of treatments: n = ", x$n, "\n", sep = "")
    cat("Number of active components: c = ", x$c, "\n", sep = "")
    if (!is.null(x$d))
      cat("Number of designs: d = ", x$d, "\n", sep = "")
    if (inherits(x, "discomb"))
      cat("Number of subnetworks: s = ", x$s, "\n", sep = "")
    ##
    cat("\n")
    ##
    ## (a) Results for comparisons
    ##
    sm <- x$sm
    reference.group <- x$reference.group
    baseline.reference <- x$baseline.reference
    ##
    comps <- sort(c(x$comps, x$inactive))
    comps.abbr <- treats(comps, nchar.comps)
    sep.comps <- x$sep.comps
    ##
    if (!backtransf & (is_relative_effect(sm) | sm == "VE"))
      sm.lab <- paste0("log", if (sm == "VE") "VR" else sm)
    else
      sm.lab <- sm
    ##
    ci.lab <- paste0(round(100 * x$level, 1), "%-CI")
    ##
    TE.common <- x$TE.common
    lowTE.common <- x$lower.common
    uppTE.common <- x$upper.common
    statistic.common <- replaceNULL(x$statistic.common, x$zval.common)
    pval.common <- x$pval.common
    ##
    TE.random <- x$TE.random
    lowTE.random <- x$lower.random
    uppTE.random <- x$upper.random
    statistic.random <- replaceNULL(x$statistic.random, x$zval.random)
    pval.random <- x$pval.random
    ##
    if (!is.null(x$seq)) {
      TE.common <- TE.common[x$seq, x$seq]
      lowTE.common <- lowTE.common[x$seq, x$seq]
      uppTE.common <- uppTE.common[x$seq, x$seq]
      statistic.common <- statistic.common[x$seq, x$seq]
      pval.common <- pval.common[x$seq, x$seq]
      ##
      if (!all(is.na(TE.random))) {
        TE.random <- TE.random[x$seq, x$seq]
        lowTE.random <- lowTE.random[x$seq, x$seq]
        uppTE.random <- uppTE.random[x$seq, x$seq]
        statistic.random <- statistic.random[x$seq, x$seq]
        pval.random <- pval.random[x$seq, x$seq]
      }
    }
    ##
    if (reference.group != "") {
      if (all(colnames(TE.common) != reference.group))
        stop("Argument 'reference.group' must match any of ",
             "the following values: ",
             paste(paste0("'", colnames(TE.common), "'"),
                   collapse = " - "))
      ##
      trts <- rownames(TE.common)
      trts.abbr <- compos(trts, comps, comps.abbr, sep.comps)
      ##
      if (baseline.reference)
        comptext <- paste0("comparison: ",
                           if (x$n == 2)
                             paste0("'",
                                    trts.abbr[trts != reference.group],
                                    "'")
                           else
                             "other treatments",
                           " vs '",
                           trts.abbr[trts == reference.group],
                           "'")
      else
        comptext <- paste0("comparison: '",
                           trts.abbr[trts == reference.group],
                           "' vs ",
                           if (x$n == 2)
                             paste0("'",
                                    trts.abbr[trts != reference.group],
                                    "'")
                           else
                             "other treatments")
      ##
      noeffect <- 1L * (backtransf & is_relative_effect(sm))
      #
      if (backtransf) {
        TE.common    <- backtransf(TE.common, sm)
        lowTE.common <- backtransf(lowTE.common, sm)
        uppTE.common <- backtransf(uppTE.common, sm)
        ##
        TE.random    <- backtransf(TE.random, sm)
        lowTE.random <- backtransf(lowTE.random, sm)
        uppTE.random <- backtransf(uppTE.random, sm)
        #
        # Switch lower and upper limit for VE if results have been
        # backtransformed
        #
        if (sm == "VE") {
          tmp.l <- lowTE.common
          lowTE.common <- uppTE.common
          uppTE.common <- tmp.l
          #
          tmp.l <- lowTE.random
          lowTE.random <- uppTE.random
          uppTE.random <- tmp.l
        }
      }
      ##
      TE.common <- round(TE.common, digits)
      lowTE.common <- round(lowTE.common, digits)
      uppTE.common <- round(uppTE.common, digits)
      statistic.common <- round(statistic.common, digits.stat)
      ##
      TE.random    <- round(TE.random, digits)
      lowTE.random <- round(lowTE.random, digits)
      uppTE.random <- round(uppTE.random, digits)
      statistic.random <- round(statistic.random, digits.stat)
      ##
      if (baseline.reference) {
        TE.common.b <- TE.common[, trts == reference.group]
        lowTE.common.b <-
          lowTE.common[, trts == reference.group]
        uppTE.common.b <-
          uppTE.common[, trts == reference.group]
        statistic.common.b <-
          statistic.common[, trts == reference.group]
        pval.common.b <-
          pval.common[, trts == reference.group]
      }
      else {
        TE.common.b <- TE.common[trts == reference.group, ]
        lowTE.common.b <-
          lowTE.common[trts == reference.group, ]
        uppTE.common.b <-
          uppTE.common[trts == reference.group, ]
        statistic.common.b <-
          statistic.common[trts == reference.group, ]
        pval.common.b <-
          pval.common[trts == reference.group, ]
      }
      ##
      dat1.c <-
        cbind(formatN(TE.common.b, digits, text.NA = "NA",
                      big.mark = big.mark),
              formatCI(formatN(round(lowTE.common.b, digits),
                               digits, "NA", big.mark = big.mark),
                       formatN(round(uppTE.common.b, digits),
                               digits, "NA", big.mark = big.mark)),
              formatN(statistic.common.b, digits.stat, text.NA = "NA",
                      big.mark = big.mark),
              formatPT(pval.common.b, digits = digits.pval,
                       scientific = scientific.pval,
                       zero = zero.pval, JAMA = JAMA.pval)
              )
      ##
      dat1.c[trts == reference.group, ] <- rep(".", ncol(dat1.c))
      dimnames(dat1.c) <-
        list(trts, c(sm.lab, ci.lab, "z", "p-value"))
      ##
      if (!is.na(TE.common.b[trts == reference.group]) &&
          TE.common.b[trts == reference.group] == noeffect)
        dat1.c[trts == reference.group, ] <- "."
      rownames(dat1.c) <- trts.abbr
      ##
      if (baseline.reference) {
        TE.random.b <- TE.random[, trts == reference.group]
        lowTE.random.b <-
          lowTE.random[, trts == reference.group]
        uppTE.random.b <-
          uppTE.random[, trts == reference.group]
        statistic.random.b <-
          statistic.random[, trts == reference.group]
        pval.random.b <-
          pval.random[, trts == reference.group]
      }
      else {
        TE.random.b <- TE.random[trts == reference.group, ]
        lowTE.random.b <-
          lowTE.random[trts == reference.group, ]
        uppTE.random.b <-
          uppTE.random[trts == reference.group, ]
        statistic.random.b <-
          statistic.random[trts == reference.group, ]
        pval.random.b <-
          pval.random[trts == reference.group, ]
      }
      ##
      dat1.r <-
        cbind(formatN(TE.random.b, digits, text.NA = "NA",
                      big.mark = big.mark),
              formatCI(formatN(round(lowTE.random.b, digits),
                               digits, "NA", big.mark = big.mark),
                       formatN(round(uppTE.random.b, digits),
                               digits, "NA", big.mark = big.mark)),
              formatN(statistic.random.b, digits.stat, text.NA = "NA",
                      big.mark = big.mark),
              formatPT(pval.random.b, digits = digits.pval,
                       scientific = scientific.pval,
                       zero = zero.pval, JAMA = JAMA.pval)
              )
      ##
      dat1.r[trts == reference.group, ] <- rep(".", ncol(dat1.r))
      dimnames(dat1.r) <-
        list(colnames(TE.random), c(sm.lab, ci.lab, "z", "p-value"))
      ##
      if (!is.na(TE.random.b[trts == reference.group]) &&
          TE.random.b[trts == reference.group] == noeffect)
        dat1.r[trts == reference.group, ] <- "."
      rownames(dat1.r) <- trts.abbr
    }
    ##
    ## (b) Results for combinations
    ##
    ci.comb.c <- data.frame(TE = x$Comb.common,
                            seTE = x$seComb.common,
                            lower = x$lower.Comb.common,
                            upper = x$upper.Comb.common,
                            statistic = x$statistic.Comb.common,
                            p = x$pval.Comb.common,
                            stringsAsFactors = FALSE)
    rownames(ci.comb.c) <- x$trts
    ## Only consider combinations with identifiable components
    ci.comb.c <- ci.comb.c[!is.na(ci.comb.c$TE), ]
    ## Drop result for inactive component (if available)
    if (!is.null(x$inactive))
      ci.comb.c <- subset(ci.comb.c, rownames(ci.comb.c) != x$inactive)
    ## Drop combinations consisting of single components
    ci.comb.c <-
      subset(ci.comb.c, grepl(sep.comps, rownames(ci.comb.c), fixed = TRUE))
    ##
    dat2.c <- formatCC(ci.comb.c,
                       backtransf, sm, x$level,
                       comps, comps.abbr, sep.comps,
                       digits, digits.stat, digits.pval,
                       scientific.pval, zero.pval, JAMA.pval,
                       big.mark)
    ##
    ci.comb.r <- data.frame(TE = x$Comb.random,
                            seTE = x$seComb.random,
                            lower = x$lower.Comb.random,
                            upper = x$upper.Comb.random,
                            statistic = x$statistic.Comb.random,
                            p = x$pval.Comb.random,
                            stringsAsFactors = FALSE)
    rownames(ci.comb.r) <- x$trts
    ## Only consider combinations with identifiable components
    ci.comb.r <- ci.comb.r[!is.na(ci.comb.r$TE), ]
    ## Drop result for inactive component (if available)
    if (!is.null(x$inactive))
      ci.comb.r <- subset(ci.comb.r, rownames(ci.comb.r) != x$inactive)
    ## Drop combinations consisting of single components
    ci.comb.r <-
      subset(ci.comb.r, grepl(sep.comps, rownames(ci.comb.r), fixed = TRUE))
    ##
    dat2.r <- formatCC(ci.comb.r,
                       backtransf, sm, x$level,
                       comps, comps.abbr, sep.comps,
                       digits, digits.stat, digits.pval,
                       scientific.pval, zero.pval, JAMA.pval,
                       big.mark)
    ##
    ## (c) Results for components
    ##
    ci.comp.c <- data.frame(TE = x$Comp.common,
                            seTE = x$seComp.common,
                            lower = x$lower.Comp.common,
                            upper = x$upper.Comp.common,
                            statistic = x$statistic.Comp.common,
                            p = x$pval.Comp.common,
                            stringsAsFactors = FALSE)
    rownames(ci.comp.c) <- x$comps
    ## Only consider components with identifiable components
    sel.comb.c <- !is.na(ci.comp.c$TE)
    ##
    dat3.c <- formatCC(ci.comp.c,
                       backtransf, sm, x$level,
                       comps, comps.abbr, sep.comps,
                       digits, digits.stat, digits.pval,
                       scientific.pval, zero.pval, JAMA.pval,
                       big.mark)
    dat3.c[!sel.comb.c, 2:4] <- "."
    ##
    ci.comp.r <- data.frame(TE = x$Comp.random,
                            seTE = x$seComp.random,
                            lower = x$lower.Comp.random,
                            upper = x$upper.Comp.random,
                            statistic = x$statistic.Comp.random,
                            p = x$pval.Comp.random,
                            stringsAsFactors = FALSE)
    rownames(ci.comp.r) <- x$comps
    ## Only consider components with identifiable components
    sel.comb.r <- !is.na(ci.comp.r$TE)
    ##
    dat3.r <- formatCC(ci.comp.r,
                       backtransf, sm, x$level,
                       comps, comps.abbr, sep.comps,
                       digits, digits.stat, digits.pval,
                       scientific.pval, zero.pval, JAMA.pval,
                       big.mark)
    dat3.r[!sel.comb.r, 2:4] <- "."
    ##
    if (common) {
      if (reference.group != "") {
        cat("Common effects model",
            if (!is.null(x$inactive))
              paste0(" (inactive component: '", x$inactive, "')"),
            "\n\n",
            sep = "")
        ##
        cat("Treatment estimate (sm = '", sm.lab, "', ", comptext, "):\n",
            sep = "")
        prmatrix(dat1.c, quote = FALSE, right = TRUE)
        cat("\n")
      }
      ##
      if (nrow(dat2.c) >= 1) {
        cat("Incremental effect for existing combinations:\n")
        print(dat2.c)
        cat("\n")
      }
      ##
      cat("Incremental effect for components:\n")
      print(dat3.c)
    }
    ##
    if (random) {
      if (common)
        cat("\n")
      if (reference.group != "") {
        cat("Random effects model",
            if (!is.null(x$inactive))
              paste0(" (inactive component: '", x$inactive, "')"),
            "\n\n",
            sep = "")
        ##
        cat("Treatment estimate (sm = '", sm.lab,
            "', ", comptext, "):\n", sep = "")
        prmatrix(dat1.r, quote = FALSE, right = TRUE)
        cat("\n")
      }
      ##
      if (nrow(dat2.r) >= 1) {
        cat("Incremental effect for existing combinations:\n")
        print(dat2.r)
        cat("\n")
      }
      ##
      cat("Incremental effect for components:\n")
      print(dat3.r)
    }
    ##
    ## (d) Heterogeneity / inconsistency
    ##
    print.I2 <- print.I2 & !is.na(I2)
    #
    if (overall.hetstat) {
      if (print.tau2 | print.tau | print.I2 | print.Q)
        cat("\n")
      #
      if (print.tau2 | print.tau | print.I2) {
        print.I2.ci <- print.I2.ci & print.I2
        #
        text.hetstat <- "Quantifying heterogeneity / inconsistency:\n"
        #
        if (print.tau2)
          text.hetstat <- paste0(
            text.hetstat,
            formatPT(x$tau^2,
                     lab = TRUE, labval = text.tau2,
                     digits = digits.tau2,
                     lab.NA = "NA", big.mark = big.mark))
        #
        if (print.tau)
          text.hetstat <- paste0(
            text.hetstat,
            if (print.tau2) "; ",
            formatPT(x$tau,
                     lab = TRUE, labval = text.tau,
                     digits = digits.tau,
                     lab.NA = "NA", big.mark = big.mark))
        #
        if (print.I2)
          text.hetstat <- paste0(
            text.hetstat,
            if (print.tau2 | print.tau) "; ",
            paste0(text.I2, " = ", round(I2, digits.I2), "%"),
            if (print.I2.ci & (!(is.na(lower.I2) | is.na(upper.I2))))
              pasteCI(lower.I2, upper.I2,
                      digits.I2, big.mark, unit = "%"))
        #
        text.hetstat <- paste0(text.hetstat, "\n")
        #
        cat(text.hetstat)
      }
      #
      if (print.Q) {
        if (print.tau2 | print.tau | print.I2)
          cat("\n")
        #
        if (inherits(x, "discomb") && x$s > 1) {
          cat("Heterogeneity statistic:\n")
          #
          hetdat <- 
            data.frame(Q = formatN(x$Q.additive, digits.Q),
                       df.Q = formatN(x$df.Q.additive, 0),
                       pval = formatPT(x$pval.Q.additive,
                                       digits = digits.pval.Q,
                                       scientific = scientific.pval,
                                       zero = zero.pval, JAMA = JAMA.pval),
                       row.names = "")
          #
          names(hetdat) <- c("Q", "df", "p-value")
          #
          print(hetdat)
        }
        else {
          cat("Heterogeneity statistics:\n")
          #
          hetdat <- 
            data.frame(Q = formatN(c(x$Q.additive,
                                     x$Q.standard,
                                     x$Q.diff),
                                   digits.Q),
                       df.Q = formatN(c(x$df.Q.additive,
                                        x$df.Q.standard,
                                        x$df.Q.diff), 0),
                       pval = formatPT(c(x$pval.Q.additive,
                                         x$pval.Q.standard,
                                         x$pval.Q.diff),
                                       digits = digits.pval.Q,
                                       scientific = scientific.pval,
                                       zero = zero.pval, JAMA = JAMA.pval),
                       row.names = c("Additive model", "Standard model",
                                     "Difference"))
          #
          names(hetdat) <- c("Q", "df", "p-value")
          #
          print(hetdat)
        }
      }
    }
    #
    # Print details of component network meta-analysis methods
    #
    if (details.methods)
      cat(textmeth(x, random, print.tau2, print.tau,
                   text.tau2, text.tau, digits.tau2, digits.tau,
                   print.I2, text.I2, big.mark))
    #
    # Add legend with abbreviated component labels
    #
    if (legend) {
      diff.comps <- comps != comps.abbr
      if (any(diff.comps)) {
        tmat <- data.frame(comps.abbr, comps)
        tmat <- tmat[diff.comps, ]
        names(tmat) <- c("Abbreviation", " Component name")
        tmat <- tmat[order(tmat$Abbreviation), ]
        ##
        cat("\nLegend:\n")
        prmatrix(tmat, quote = FALSE, right = TRUE,
                 rowlab = rep("", length(comps.abbr))) 
      }
    }
  }
  
  invisible(NULL)
}
