#' Calculate comparison effects of two arbitrary complex interventions
#' in component network meta-analysis
#' 
#' @description
#' Calculate comparison effects of two arbitrary complex interventions
#' (i.e., combinations of several components) in component network
#' meta-analysis.
#' 
#' @param x An object of class \code{netcomb} or \code{netcomparison}
#'   (print function).
#' @param treat1 A character vector defining the first complex
#'   intervention(s).
#' @param treat2 A character vector defining the second complex
#'   intervention(s).
#' @param common A logical indicating whether results for common
#'   effects model should be conducted.
#' @param random A logical indicating whether results for random
#'   effects model should be conducted.
#' @param level The level used to calculate confidence intervals for
#'   combinations of components.
#' @param nchar.comps A numeric defining the minimum number of
#'   characters used to create unique names for components (see
#'   Details).
#' @param backtransf A logical indicating whether printed results
#'   should be back transformed. If \code{backtransf=TRUE}, results
#'   for \code{sm="OR"} are printed as odds ratios rather than log
#'   odds ratios.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.stat Minimal number of significant digits for z-value
#'   of test for overall effect, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for
#'   p-values, see \code{print.default}.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   combination effect should be printed according to JAMA reporting
#'   standards.
#' @param big.mark A character used as thousands separator.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
#'
#' @details
#' R functions \code{\link{netcomb}} and \code{\link{discomb}}
#' calculate effects for individual components and complex
#' interventions present in the component network meta-analysis
#' (CNMA). This function can be used to calculate the effect for
#' comparisons of two arbitrary complex interventions defined by
#' arguments \code{treat1} and \code{treat2}.
#'
#' All complex interventions occuring in the network are considered
#' for the first complex intervention if argument \code{treat1} is
#' missing. The reference group defined in the (C)NMA is used as
#' second complex intervention if argument \code{treat2} is
#' missing. The first complex intervention in the (C)NMA is used if
#' the reference group is not defined.
#' 
#' The following matrices are needed to calculate comparison effects
#' of arbitrary complex interventions, (Rücker et al., 2020, Section
#' 3.2):
#' \itemize{
#' \item B matrix describing how comparisons are composed by complex
#'   intervetions,
#' \item C matrix describing how the complex interventions are
#'   composed by the components.
#' }
#' Internally, both matrices are constructed based on arguments
#' \code{x}, \code{treat1} and \code{treat2}.
#' 
#' By default, component names are not abbreviated in
#' printouts. However, in order to get more concise printouts,
#' argument \code{nchar.comps} can be used to define the minimum
#' number of characters for abbreviated component names (see
#' \code{\link{abbreviate}}, argument \code{minlength}). R function
#' \code{\link{treats}} is utilised internally to create abbreviated
#' component names.
#' 
#' @note
#' R function \code{\link{netcomplex}} can be used to calculate the
#' effect for arbitrary complex interventions in a component network
#' meta-analysis.
#' 
#' @return
#' A list is returned by the function \code{netcomparison} with the
#' following elements:
#' \item{comparison}{Comparison.}
#' \item{TE.common, TE.random}{A vector of comparison effects (common
#'   and random effects model).}
#' \item{seTE.common, seTE.random}{A vector with corresponding standard
#'   errors (common and random effects model).}
#' \item{lower.common, lower.random}{A vector with lower confidence
#'   limits for comparisons (common and random effects model).}
#' \item{upper.common, upper.random}{A vector with upper confidence
#'   limits for comparisons (common and random effects model).}
#' \item{statistic.common, statistic.random}{A vector with z-values for
#'   the overall effect of comparisons (common and random effects
#'   model).}
#' \item{pval.common, pval.random}{A vector with p-values for the
#'   overall effect of comparisons (common and random effects model).}
#' \item{trts}{Treatments included in comparisons.}
#' \item{comps}{Components included in comparisons.}
#' \item{treat1, treat2}{A defined above.}
#' \item{common, random}{A defined above.}
#' \item{level, nchar.comps, backtransf, x}{A defined above.}
#' \item{B.matrix}{B matrix.}
#' \item{C.matrix}{C matrix.}
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netcomb}}, \code{\link{discomb}},
#'   \code{\link{netcomplex}}
#' 
#' @references
#' Rücker G, Petropoulou M, Schwarzer G (2020):
#' Network meta-analysis of multicomponent interventions.
#' \emph{Biometrical Journal},
#' \bold{62}, 808--21
#' 
#' @examples
#' data(Linde2016)
#' 
#' # Only consider studies including Face-to-face PST (to reduce
#' # runtime of example)
#' #
#' face <- subset(Linde2016, id %in% c(16, 24, 49, 118))
#' 
#' # Conduct random effects network meta-analysis
#' #
#' net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
#'   data = face, ref = "placebo", sm = "OR", common = FALSE)
#' 
#' # Additive model for treatment components (with placebo as inactive
#' # treatment)
#' #
#' nc1 <- netcomb(net1, inactive = "placebo")
#' 
#' # Result for comparison Face-to-face PST vs TCA
#' netcomparison(nc1, "Face-to-face PST", "TCA", nchar.comps = 4)
#' netcomparison(nc1, "F", "T", nchar.comps = 4)
#'
#' # Result for comparison Face-to-face PST vs TCA + Placebo
#' netcomparison(nc1, "Face-to-face PST", "TCA + Plac", nchar.comps = 4)
#' 
#' \donttest{
#'
#' # Forest plots
#' #
#' t1 <- c("F + TCA", "F + Plac", "SSRI + Plac + TCA")
#' t2 <- c("UC", "Plac", "UC")
#' #
#' netcomparison(nc1, t1, t2)
#' #
#' forest(netcomparison(nc1, t1, t2))
#' forest(netcomparison(nc1, t1, t2), nchar.comps = 4)
#' forest(netcomparison(nc1, c("F", "TCA"), "UC"), nchar.comps = 4)
#' 
#' # Artificial example
#' #
#' treat1 <- rep("A", 3)
#' treat2 <- c("B+C", "A+C", "C+D")
#' TE <- c(0, 1, 0)
#' seTE <- rep(1, 3)
#' # Conduct (C)NMA
#' net2 <- netmeta(TE, seTE, treat1, treat2, random = FALSE)
#' nc2 <- netcomb(net2)
#'
#' # Result for comparison A vs B + D
#' netcomparison(nc2, "A", "B + D")
#' # Same results
#' netcomparison(nc2, "A", "B+D")
#' netcomparison(nc2, "A", "D+B")
#' netcomparison(nc2, "a", "d+b")
#' 
#' # Generated B matrix
#' netcomparison(nc2, "A", "B + D")$C.matrix
#' # Generated B matrix
#' netcomparison(nc2, "A", "B + D")$B.matrix
#' }
#' 
#' @rdname netcomparison
#' @export
#' @export netcomparison

netcomparison <- function(x, treat1, treat2,
                          common = x$common,
                          random = x$random,
                          level = x$level.ma,
                          nchar.comps = x$nchar.comps,
                          backtransf = x$backtransf,
                          warn.deprecated = gs("warn.deprecated"),
                          ...) {
  
  
  chkclass(x, "netcomb")
  x <- updateversion(x)
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  common <- deprecated(common, missing(common), args, "fixed",
                       warn.deprecated)
  chklogical(common)
  ##
  chklogical(random)
  chklevel(level)
  ##
  chknumeric(nchar.comps, min = 1, length = 1)
  chklogical(backtransf)
  
  
  missing.treat1 <- missing(treat1)
  missing.treat2 <- missing(treat2)
  ##
  if (missing(treat1)) {
    treat1.orig <- NULL
    treat2.orig <- NULL
    treat1 <- x$trts
  }
  else
    treat1.orig <- treat1
  ##
  if (missing(treat2)) {
    treat2.orig <- NULL
    treat2 <- x$reference.group
    if (treat2 == "") {
      treat2 <- treat1[1]
      treat1 <- treat1[-1]
    }
    else
      treat1 <- treat1[treat1 != treat2]
  }
  else
    treat2.orig <- treat2
  ##
  if (length(treat1) == 1 & length(treat2) > 1)
    treat1 <- rep(treat1, length(treat2))
  else if (length(treat1) > 1 & length(treat2) == 1)
    treat2 <- rep(treat2, length(treat1))
  ##
  n.comparisons <- length(treat1)
  ##
  chklength(treat2, n.comparisons, "treat1")
  ##
  comps <- x$comps
  
  
  ##
  ## Check components
  ##
  comps1.list <- compsplit(treat1, x$sep.comps)
  comps2.list <- compsplit(treat2, x$sep.comps)
  ##
  comps1 <- setref(unique(unlist(comps1.list)), c(comps, x$inactive),
                   error.text = "component names in argument 'treat1'",
                   length = 0)
  comps2 <- setref(unique(unlist(comps2.list)), c(comps, x$inactive),
                   error.text = "component names in argument 'treat2'",
                   length = 0)
  ##
  comps1.list <- lapply(comps1.list, setref, c(comps, x$inactive), length = 0)
  comps2.list <- lapply(comps2.list, setref, c(comps, x$inactive), length = 0)
  ##
  add1 <- add2 <- rep("", n.comparisons)
  ##
  for (i in seq_len(n.comparisons)) {
    add1[i] <-
      if (attr(compsplit(treat1[i], x$sep.comps), "withspace")) " " else ""
    treat1[i] <-
      paste(comps1.list[[i]], collapse = paste0(add1[i], x$sep.comps, add1[i]))
  }
  ##
  for (i in seq_len(n.comparisons)) {
    add2[i] <-
      if (attr(compsplit(treat2[i], x$sep.comps), "withspace")) " " else ""
    treat2[i] <-
      paste(comps2.list[[i]], collapse = paste0(add2[i], x$sep.comps, add2[i]))
  }
  ##
  trts <- sort(unique(c(treat1, treat2)))
  
  
  ##
  ## Extract comparisons
  ##
  comparison <- rep("", n.comparisons)
  ##
  for (i in seq_len(n.comparisons)) {
    sel1.i <- !comps1.list[[i]] %in% comps2.list[[i]]
    sel2.i <- !comps2.list[[i]] %in% comps1.list[[i]]
    #
    if (any(sel1.i) | any(sel2.i)) {
      comb1 <- paste(comps1.list[[i]][sel1.i],
                     collapse = paste0(add1[i], x$sep.comps, add1[i]))
      comb2 <- paste(comps2.list[[i]][sel2.i],
                     collapse = paste0(add2[i], x$sep.comps, add2[i]))
      #
      if (comb1 == "")
        comparison[i] <- comb2
      else if (comb2 == "")
        comparison[i] <- comb1
      else
        comparison[i] <- paste0(comb1, x$sep.trts, comb2)
    }
  }
  
  
  ##
  ## Generate C matrix
  ##
  C.matrix <- matrix(0, nrow = length(trts), ncol = length(comps))
  rownames(C.matrix) <- trts
  colnames(C.matrix) <- comps
  ##
  for (i in seq_len(n.comparisons)) {
    C.matrix[treat1[i], ] <- 1L * colnames(C.matrix) %in% comps1.list[[i]]
    C.matrix[treat2[i], ] <- 1L * colnames(C.matrix) %in% comps2.list[[i]]
  }
  
  
  ##
  ## Generate B matrix
  ##
  B.matrix <- matrix(0, nrow = n.comparisons, ncol = length(trts))
  rownames(B.matrix) <- seq_len(n.comparisons)
  colnames(B.matrix) <- trts
  ##
  for (i in seq_len(n.comparisons)) {
    ##
    B.matrix[i, ] <-
      1L * colnames(B.matrix) %in% treat1[i] +
      -1L * colnames(B.matrix) %in% treat2[i]
  }
  
  
  ##
  ## Calculate estimates for comparisons
  ##
  X.matrix <- B.matrix %*% C.matrix
  #
  # Identify inestimable components and drop them from the design matrix
  #
  compNA <- is.na(x$Comp.common)
  inestimable <- apply(abs(X.matrix[, compNA, drop = FALSE]), 1, sum) > 0
  X.matrix <- X.matrix[, !compNA, drop = FALSE]
  #
  TE.common <- as.vector(X.matrix %*% x$Comp.common[!compNA])
  seTE.common <-
    sqrt(
      diag(X.matrix %*% x$Lplus.matrix.common[!compNA, !compNA] %*%
             t(X.matrix)))
  #
  TE.random <- as.vector(X.matrix %*% x$Comp.random[!compNA])
  seTE.random <-
    sqrt(diag(X.matrix %*% x$Lplus.matrix.random[!compNA, !compNA] %*%
                t(X.matrix)))
  #
  # Set comparisons with inestimable components to NA
  #
  TE.common[inestimable] <- NA
  seTE.common[inestimable] <- NA
  TE.random[inestimable] <- NA
  seTE.random[inestimable] <- NA
  #
  # Calculate confidence intervals
  #
  ci.f <- ci(TE.common, seTE.common, level = level)
  ci.r <- ci(TE.random, seTE.random, level = level)
  
  
  res <- list(comparison = comparison,
              treat1 = treat1, treat2 = treat2,
              ##
              TE.common = ci.f$TE,
              seTE.common = ci.f$seTE,
              lower.common = ci.f$lower,
              upper.common = ci.f$upper,
              statistic.common = ci.f$statistic,
              pval.common = ci.f$p,
              ##
              TE.random = ci.r$TE,
              seTE.random = ci.r$seTE,
              lower.random = ci.r$lower,
              upper.random = ci.r$upper,
              statistic.random = ci.r$statistic,
              pval.random = ci.r$p,
              ##
              common = common,
              random = random,
              level = level,
              ##
              trts = trts,
              comps = colnames(C.matrix)[apply(C.matrix, 2, sum) > 0],
              inactive = x$inactive,
              nchar.comps = nchar.comps,
              backtransf = backtransf,
              ##
              B.matrix = B.matrix,
              C.matrix = C.matrix,
              ##
              x = x,
              ##
              add1 = add1,
              add2 = add2,
              ##
              treat1.orig = treat1.orig,
              treat2.orig = treat2.orig,
              ##
              version = packageDescription("netmeta")$Version
              )
  ##
  ## Backward compatibility
  ##
  res$fixed <- res$common
  ##
  res$TE.fixed <- res$TE.common
  res$seTE.fixed <- res$seTE.common
  res$lower.fixed <- res$lower.common
  res$upper.fixed <- res$upper.common
  res$statistic.fixed <- res$statistic.common
  res$pval.fixed <- res$pval.common
  ##
  class(res) <- c("netcomparison", class(res))
  
  res
}


#' @rdname netcomparison
#' @method print netcomparison
#' @export

print.netcomparison <- function(x,
                                ##
                                common = x$common,
                                random = x$random,
                                backtransf = x$backtransf,
                                ##
                                nchar.comps = x$nchar.comps,
                                ##
                                digits = gs("digits"),
                                digits.stat = gs("digits.stat"),
                                digits.pval = gs("digits.pval"),
                                ##
                                scientific.pval = gs("scientific.pval"),
                                zero.pval = gs("zero.pval"),
                                JAMA.pval = gs("JAMA.pval"),
                                big.mark = gs("big.mark"),
                                ##
                                legend = gs("legend"),
                                warn.deprecated = gs("warn.deprecated"),
                                ...) {
  
  chkclass(x, "netcomparison")
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  common <- deprecated(common, missing(common), args, "fixed",
                       warn.deprecated)
  chklogical(common)
  ##
  chklogical(random)
  chklogical(backtransf)
  ##
  chknumeric(nchar.comps, min = 1, length = 1)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  ##
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  ##
  chklogical(legend)
  
  
  ##
  ## Abbreviated component and treatment labels
  ##
  n.comparisons <- length(x$treat1)
  treat1 <- rep("", n.comparisons)
  treat2 <- rep("", n.comparisons)
  ##
  if (common | random) {
    comps <- x$comps
    if (!is.null(x$inactive) && any(grepl(x$inactive, x$comparison)))
      comps <- c(comps, x$inactive)
    comps.abbr <- treats(comps, nchar.comps)
    ##
    for (i in seq_len(n.comparisons))
      treat1[i] <- compos(x$treat1[i], comps, comps.abbr,
                          x$x$sep.comps, x$add1[i] == " ")
    for (i in seq_len(n.comparisons))
      treat2[i] <- compos(x$treat2[i], comps, comps.abbr,
                          x$x$sep.comps, x$add2[i] == " ")
    ##
    comps1.list <- compsplit(x$treat1, x$x$sep.comps)
    comps2.list <- compsplit(x$treat2, x$x$sep.comps)
    ##
    comps1.list <- lapply(comps1.list, charfac,
                          levels = comps, labels = comps.abbr)
    comps2.list <- lapply(comps2.list, charfac,
                          levels = comps, labels = comps.abbr)
    ##
    comparison <- rep("", n.comparisons)
    ##
    for (i in seq_len(n.comparisons)) {
      sel1.i <- !comps1.list[[i]] %in% comps2.list[[i]]
      sel2.i <- !comps2.list[[i]] %in% comps1.list[[i]]
      ##
      if (any(sel1.i) | any(sel2.i)) {
        comb1 <- paste(comps1.list[[i]][sel1.i],
                       collapse = paste0(x$add1[i], x$x$sep.comps, x$add1[i]))
        comb2 <- paste(comps2.list[[i]][sel2.i],
                       collapse = paste0(x$add2[i], x$x$sep.comps, x$add2[i]))
        #
        if (comb1 == "")
          comparison[i] <- comb2
        else if (comb2 == "")
          comparison[i] <- comb1
        else
          comparison[i] <- paste0(comb1, x$x$sep.trts, comb2)
      }
    }
  }
  
  
  sm <- sm.lab <- x$x$sm
  ##
  relative <- is_relative_effect(sm) | sm == "VE"
  #
  if (!backtransf & relative)
    sm.lab <- paste0("log", if (sm == "VE") "VR" else sm)
  ##  
  ci.lab <- paste0(round(100 * x$level, 1), "%-CI")
  
  
  if (common) {
    TE.common <- x$TE.common
    lower.common <- x$lower.common
    upper.common <- x$upper.common
    ##
    if (backtransf) {
      TE.common <- backtransf(TE.common, sm)
      lower.common <- backtransf(lower.common, sm)
      upper.common <- backtransf(upper.common, sm)
      #
      # Switch lower and upper limit for VE if results have been
      # backtransformed
      #
      if (sm == "VE") {
        tmp.l <- lower.common
        lower.common <- upper.common
        upper.common <- tmp.l
      }
    }
  }
  ##
  if (random) {
    TE.random <- x$TE.random
    lower.random <- x$lower.random
    upper.random <- x$upper.random
    ##
    if (backtransf) {
      TE.random <- backtransf(TE.random, sm)
      lower.random <- backtransf(lower.random, sm)
      upper.random <- backtransf(upper.random, sm)
      #
      # Switch lower and upper limit for VE if results have been
      # backtransformed
      #
      if (sm == "VE") {
        tmp.l <- lower.random
        lower.random <- upper.random
        upper.random <- tmp.l
      }
    }
  }
  
  
  if (common) {
    pval.f <- formatPT(x$pval.common, digits = digits.pval,
                       scientific = scientific.pval,
                       zero = zero.pval, JAMA = JAMA.pval,
                       lab.NA = "")
    ##
    res.c <-
      cbind(comparison = comparison,
            treat1 = treat1,
            treat2 = treat2,
            Comb = formatN(TE.common, digits = digits,
                           "NA", big.mark = big.mark),
            CI = formatCI(formatN(round(lower.common, digits),
                                  digits, "NA",
                                  big.mark = big.mark),
                          formatN(round(upper.common, digits),
                                  digits, "NA",
                                  big.mark = big.mark)),
            zval = formatN(x$statistic.common, digits = digits.stat,
                           "NA", big.mark = big.mark),
            pval = pval.f)
    ##
    dimnames(res.c) <- list(rep("", nrow(res.c)),
                            c("comparison", "treat1", "treat2",
                              sm.lab, ci.lab, "z", "p-value"))
    ##
    res.c[res.c[, "treat1"] == res.c[, "treat2"], 4:7] <- "."
    ##
    cat("Results for comparisons (additive CNMA model,",
        "common effects model):\n\n")
    ##
    prmatrix(res.c, quote = FALSE, right = TRUE, na.print = "--")
    ##
    if (random)
      cat("\n")
  }
  
  
  if (random) {
    ##
    pval.r <- formatPT(x$pval.random, digits = digits.pval,
                       scientific = scientific.pval,
                       zero = zero.pval, JAMA = JAMA.pval,
                       lab.NA = "")
    ##
    res.r <-
      cbind(comparison = comparison,
            treat1 = treat1,
            treat2 = treat2,
            Comb = formatN(TE.random, digits = digits,
                           "NA", big.mark = big.mark),
            CI = formatCI(formatN(round(lower.random, digits),
                                  digits, "NA",
                                  big.mark = big.mark),
                          formatN(round(upper.random, digits),
                                  digits, "NA",
                                  big.mark = big.mark)),
            zval = formatN(x$statistic.random, digits = digits.stat,
                           "NA", big.mark = big.mark),
            pval = pval.r)
    ##
    dimnames(res.r) <- list(rep("", nrow(res.r)),
                            c("comparison", "treat1", "treat2",
                              sm.lab, ci.lab, "z", "p-value"))
    ##
    res.r[res.r[, "treat1"] == res.r[, "treat2"], 4:7] <- "."
    ##
    cat("Results for comparisons (additive CNMA model,",
        "random effects model):\n\n")
    ##
    prmatrix(res.r, quote = FALSE, right = TRUE, na.print = "--")
  }
  
  
  if (legend && (common | random)) {
    diff.comps <- comps != comps.abbr
    any.comps <- any()
    ##
    if (any(diff.comps)) {
      cat("\nLegend:\n")
      ##
      tmat <- data.frame(comps.abbr, comps)
      tmat <- tmat[diff.comps, ]
      names(tmat) <- c("Abbreviation", " Component name")
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      prmatrix(tmat, quote = FALSE, right = TRUE,
               rowlab = rep("", length(comps.abbr))) 
    }
  }
  
  invisible(NULL)
}
