formatCC <- function(x,
                     backtransf, sm, level,
                     comps, comps.abbr, sep.comps,
                     digits, digits.stat, digits.pval.Q,
                     scientific.pval, zero.pval, JAMA.pval,
                     big.mark) {
  
  
  relative <- is.relative.effect(sm)
  
  sm.lab <- sm
  ##
  if (sm != "") {
    sm.lab <- paste0("i", sm)
    if (!backtransf & relative)
      sm.lab <- paste0("log(", sm.lab, ")")
  }
  ##  
  ci.lab <- paste(round(100 * level, 1), "%-CI", sep = "")
  
  
  ## First column contains row names
  ##
  res <- x
  ##
  rownames(res) <- compos(rownames(res), comps, comps.abbr, sep.comps)
  ##
  if (backtransf & relative) {
    res$TE <- exp(res$TE)
    res$lower <- exp(res$lower)
    res$upper <- exp(res$upper)
  }
  ##
  res$TE <- formatN(res$TE, digits, "NA", big.mark)
  res$lower <- formatCI(formatN(round(res$lower, digits),
                                digits, "NA", big.mark),
                        formatN(round(res$upper, digits),
                                digits, "NA", big.mark))
  res$statistic <- formatN(round(res$statistic, digits.stat),
                           digits.stat, big.mark = big.mark)
  res$p <- formatPT(res$p,
                    digits = digits.pval.Q,
                    scientific = scientific.pval,
                    zero = zero.pval,
                    JAMA = JAMA.pval)
  ##
  res$seTE <- res$upper <- res$z <- NULL
  ##
  sel <- names(res) == "TE"
  names(res)[sel] <- sm.lab
  ##
  sel <- names(res) == "lower"
  names(res)[sel] <- ci.lab
  ##
  sel <- names(res) == "statistic"
  names(res)[sel] <- "z"
  ##
  sel <- names(res) == "p"
  names(res)[sel] <- "p-value"
  
  
  res
}
