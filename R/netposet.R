#' Partial order of treatments in network meta-analysis
#' 
#' @description
#' Partial order of treatments in network meta-analysis. The set of
#' treatments in a network is called a partially ordered set (in
#' short, a \emph{poset}), if different outcomes provide different
#' treatment ranking lists.
#' 
#' @aliases netposet print.netposet
#' 
#' @param \dots See details.
#' @param outcomes A character vector with outcome names.
#' @param treatments A character vector with treatment names.
#' @param small.values See details.
#' @param common A logical indicating whether to show results for the
#'   common effects model.
#' @param random A logical indicating whether to show results for the
#'   random effects model.
#' @param x An object of class \code{netposet}.
#' @param pooled A character string indicating whether Hasse diagram
#'   should be drawn for common (\code{"common"}) or random effects
#'   model (\code{"random"}). Can be abbreviated.
#' @param fixed Ignored deprecated argument (replaced by
#'   \code{common}).
#' @param comb.fixed Ignored deprecated argument (replaced by
#'   \code{common}).
#' @param comb.random Ignored deprecated argument (replaced by
#'   \code{random}).
#' 
#' @details
#' In network meta-analysis, frequently different outcomes are
#' considered which may each provide a different ordering of
#' treatments. The concept of a partially ordered set (in short, a
#' \emph{poset}, Carlsen & Bruggemann, 2014) of treatments can be used
#' to gain further insights in situations with apparently conflicting
#' orderings. This implementation for rankings in network meta-analyis
#' is described in Rücker & Schwarzer (2017).
#' 
#' In function \code{netposet}, argument \code{\dots{}} can be any of the
#' following:
#' \itemize{
#' \item arbitrary number of \code{netrank} objects providing
#'   P-scores;
#' \item arbitrary number of \code{netmeta} objects;
#' \item single ranking matrix with each column providing P-scores
#'   (Rücker & Schwarzer 2015) or SUCRA values (Salanti et al. 2011)
#'   for an outcome and rows corresponding to treatments.
#' }
#' Note, albeit in general a ranking matrix is not constrained to have
#' values between 0 and 1, \code{netposet} stops with an error in this
#' case as this function expects a matrix with P-scores or SUCRA
#' values.
#' 
#' Argument \code{outcomes} can be used to label outcomes. If argument
#' \code{outcomes} is missing,
#' \itemize{
#' \item column names of the ranking matrix are used as outcome labels
#'   (if first argument is a ranking matrix and column names are
#'   available);
#' \item capital letters 'A', 'B', \dots{} are used as outcome labels
#'   and a corresponding warning is printed.
#' }
#'
#' Argument \code{treatments} can be used to provide treatment labels
#' if the first argument is a ranking matrix. If argument
#' \code{treatment} is missing,
#' \itemize{
#' \item row names of the ranking matrix are used as treatment labels
#'   (if available);
#' \item letters 'a', 'b', \dots{} are used as treatment labels and a
#'   corresponding warning is printed.
#' }
#' If argument \code{\dots{}} consists of \code{netmeta} objects,
#' \code{netrank} is called internally to calculate P-scores. In this
#' case, argument \code{small.values} can be used to specify for each
#' outcome whether small values are beneficial (\code{"desirable"}) or
#' harmfull (\code{"undesirable"}); see \code{\link{netrank}}. This
#' argument is ignored for a ranking matrix and \code{netrank}
#' objects.
#' 
#' Arguments \code{common} and \code{random} can be used to define
#' whether results should be printed and plotted for common and random
#' effects model. If netmeta and netrank objects are provided in
#' argument \code{\dots{}}, values for \code{common} and \code{random}
#' within these objects are considered; if these values are not
#' unique, argument \code{common} or \code{random} are set to
#' \code{TRUE}.
#' 
#' In function \code{print.netposet}, argument \code{\dots{}} is
#' passed on to the printing function.
#'
#' @return
#' An object of class \code{netposet} with corresponding \code{print},
#' \code{plot}, and \code{hasse} functions. The object is a list
#' containing the following components:
#' \item{P.common}{Ranking matrix with rows corresponding to treatments
#'   and columns corresponding to outcomes (common effects model).}
#' \item{M0.common}{Hasse matrix skipping unnecessary paths (common
#'   effects model).}
#' \item{M.common}{"Full" Hasse matrix (common effects model).}
#' \item{O.common}{Matrix with information about partial ordering
#'   (common effects model).}
#' \item{P.random}{Ranking matrix with rows corresponding to
#'   treatments and columns corresponding to outcomes (random effects
#'   model).}
#' \item{M0.random}{Hasse matrix skipping unnecessary paths (random
#'   effects model).}
#' \item{M.random}{"Full" Hasse matrix (random effects model).}
#' \item{O.random}{Matrix with information about partial ordering
#'   (random effects model).}
#' \item{small.values, common, random}{As.defined above.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package netmeta used to create object.}
#' 
#' @author Gerta Rücker \email{gerta.ruecker@@uniklinik-freiburg.de}, Guido
#'   Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{netrank}},
#'   \code{\link{plot.netrank}}, \code{\link{hasse}},
#'   \code{\link{plot.netposet}}, \code{\link[metadat]{dat.linde2015}}
#' 
#' @references
#' Carlsen L, Bruggemann R (2014):
#' Partial order methodology: a valuable tool in chemometrics.
#' \emph{Journal of Chemometrics},
#' \bold{28}, 226--34
#' 
#' Rücker G, Schwarzer G (2015):
#' Ranking treatments in frequentist network meta-analysis works
#' without resampling methods.
#' \emph{BMC Medical Research Methodology},
#' \bold{15}, 58
#' 
#' Rücker G, Schwarzer G (2017):
#' Resolve conflicting rankings of outcomes in network meta-analysis:
#' Partial ordering of treatments.
#' \emph{Research Synthesis Methods},
#' \bold{8}, 526--36
#' 
#' Salanti G, Ades AE, Ioannidis JP (2011):
#' Graphical methods and numerical summaries for presenting results
#' from multiple-treatment meta-analysis: an overview and tutorial.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{64}, 163--71
#' 
#' @examples
#' \donttest{
#' # Define order of treatments in depression data set dat.linde2015
#' #
#' trts <- c("TCA", "SSRI", "SNRI", "NRI",
#'   "Low-dose SARI", "NaSSa", "rMAO-A", "Hypericum", "Placebo")
#' 
#' # Outcome labels
#' #
#' outcomes <- c("Early response", "Early remission")
#' 
#' # (1) Early response
#' #
#' pw1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(resp1, resp2, resp3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#' #
#' net1 <- netmeta(pw1, common = FALSE,
#'   seq = trts, ref = "Placebo", small.values = "undesirable")
#' 
#' # (2) Early remission
#' #
#' pw2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(remi1, remi2, remi3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#' #
#' net2 <- netmeta(pw2, common = FALSE,
#'   seq = trts, ref = "Placebo", small.values = "undesirable")
#' 
#' # Partial order of treatment rankings (two outcomes)
#' #
#' po12 <- netposet(netrank(net1), netrank(net2), outcomes = outcomes)
#' 
#' # Scatter plot
#' #
#' plot(po12)
#' 
#' # Biplot (same scatter plot as only two outcomes considered)
#' #
#' plot(po12, "biplot")
#' 
#' # Hasse diagram for outcomes early response and early remission
#' #
#' \dontrun{
#' if (requireNamespace("Rgraphviz", quietly = TRUE))
#'   hasse(po12)
#' }
#' 
#' # Scatter plot for outcomes early response and early remission
#' #
#' oldpar <- par(pty = "s")
#' plot(po12)
#' par(oldpar)
#'
#' 
#' # Consider five outcomes
#' #
#' # Outcome labels
#' #
#' outcomes <- c("Early response", "Early remission",
#'   "Lost to follow-up", "Lost to follow-up due to AEs",
#'    "Adverse events (AEs)")
#' 
#' # (3) Loss to follow-up
#' #
#' pw3 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss1, loss2, loss3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#' #
#' net3 <- netmeta(pw3, common = FALSE,
#'   seq = trts, ref = "Placebo", small.values = "desirable")
#' 
#' # (4) Loss to follow-up due to adverse events
#' #
#' pw4 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss.ae1, loss.ae2, loss.ae3), n = list(n1, n2, n3),
#'   studlab = id, data = subset(dat.linde2015, id != 55), sm = "OR")
#' #
#' net4 <- netmeta(pw4, common = FALSE,
#'   seq = trts, ref = "Placebo", small.values = "desirable")
#' 
#' # (5) Adverse events
#' #
#' pw5 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(ae1, ae2, ae3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#' #
#' net5 <- netmeta(pw5, common = FALSE,
#'   seq = trts, ref = "Placebo", small.values = "desirable")
#' 
#' # Partial order of treatment rankings (based on netrank() objects)
#' #
#' po.ranks <- netposet(netrank(net1), netrank(net2),
#'   netrank(net3), netrank(net4), netrank(net5), outcomes = outcomes)
#' 
#' # Same result (based on netmeta() objects)
#' #
#' po.nets <- netposet(net1, net2, net3, net4, net5,
#'   outcomes = outcomes)
#' #
#' all.equal(po.ranks, po.nets)
#' 
#' # Print matrix with P-scores (random effects model)
#' #
#' po.nets$P.random
#' 
#' \dontrun{
#' # Hasse diagram for all outcomes (random effects model)
#' #
#' if (requireNamespace("Rgraphviz", quietly = TRUE))
#'   hasse(po.ranks)
#' }
#' }
#' 
#' # Example using ranking matrix with P-scores
#' #
#' # Ribassin-Majed L, Marguet S, Lee A.W., et al. (2017):
#' # What is the best treatment of locally advanced nasopharyngeal
#' # carcinoma? An individual patient data network meta-analysis.
#' # Journal of Clinical Oncology, 35, 498-505
#' #
#' outcomes <- c("OS", "PFS", "LC", "DC")
#' treatments <- c("RT", "IC-RT", "IC-CRT", "CRT",
#'   "CRT-AC", "RT-AC", "IC-RT-AC")
#' #
#' # P-scores (from Table 1)
#' #
#' pscore.os  <- c(15, 33, 63, 70, 96, 28, 45) / 100
#' pscore.pfs <- c( 4, 46, 79, 52, 94, 36, 39) / 100
#' pscore.lc  <- c( 9, 27, 47, 37, 82, 58, 90) / 100
#' pscore.dc  <- c(16, 76, 95, 48, 72, 32, 10) / 100
#' #
#' pscore.matrix <- data.frame(pscore.os, pscore.pfs, pscore.lc, pscore.dc)
#' rownames(pscore.matrix) <- treatments
#' colnames(pscore.matrix) <- outcomes
#' pscore.matrix
#' #
#' po <- netposet(pscore.matrix)
#' po12 <- netposet(pscore.matrix[, 1:2])
#' po
#' po12
#' #
#' \dontrun{
#' if (requireNamespace("Rgraphviz", quietly = TRUE)) {
#'   hasse(po)
#'   hasse(po12)
#' }
#' }
#' #
#' oldpar <- par(pty = "s")
#' plot(po12)
#' par(oldpar)
#' 
#' @rdname netposet
#' @export netposet

netposet <- function(..., outcomes, treatments, small.values,
                     common, random, fixed, comb.fixed, comb.random) {
  
  
  args <- list(...)
  
  
  if (!missing(common))
    chklogical(common)
  if (!missing(random))
    chklogical(random)
  ##
  missing.small.values <- missing(small.values)
  
  
  any.netmeta <- any.netrank <- FALSE
  
  
  ##
  ## First argument is expected to be a ranking matrix if only a
  ## single argument is provided
  ##
  if (length(args) == 1) {
    ##
    ## (1) Check arguments
    ##
    if (!(is.matrix(args[[1]]) | is.data.frame(args[[1]])))
      stop("First argument must be a matrix or ",
           "data frame if only a single argument is provided.",
           call. = FALSE)
    ##
    if (!missing.small.values)
      warning("Argument 'small.values' is ignored as first argument ",
              "is a ranking matrix.",
              call. = FALSE)
    ##
    ranking.matrix <- args[[1]]
    ##
    if (any(ranking.matrix[!is.na(ranking.matrix)] > 1) |
        any(ranking.matrix[!is.na(ranking.matrix)] < 0))
      stop("All elements of ranking matrix must be between 0 and 1.",
           call. = FALSE)
    ##
    n.outcomes <- ncol(ranking.matrix)
    n.treatments <- nrow(ranking.matrix)
    ##
    if (n.outcomes == 1)
      stop("Minimum number of two outcomes ",
           "(equal to number of columns of ranking matrix).",
           call. = FALSE)
    ##
    if (!missing(outcomes)) {
      if (length(outcomes) != n.outcomes)
        stop("Number of outcomes provided by argument 'outcomes' ",
             "is different from number of columns of ranking matrix.",
             call. = FALSE)
    }
    else
      outcomes <- colnames(ranking.matrix)
    ##
    if (!missing(treatments)) {
      if (length(treatments) != n.treatments)
        stop("Number of treatments provided by argument 'treatments' ",
             "is different from number of rows of ranking matrix.",
             call. = FALSE)
    }
    else
      treatments <- rownames(ranking.matrix)
    ##
    ## (2) P-Score matrices
    ##
    if (is.null(outcomes)) {
      warning("Outcomes are labelled 'A' to '", LETTERS[n.outcomes],
              "' as argument 'outcomes' is missing and column names are NULL.",
              call. = FALSE)
      outcomes <- LETTERS[1:n.outcomes]
    }
    ##
    if (is.null(treatments)) {
      warning("Treatments are labelled 'a' to '", letters[n.treatments],
              "' as row names are NULL.",
              call. = FALSE)
      treatments <- letters[1:n.treatments]
    }
    ##
    rownames(ranking.matrix) <- treatments
    colnames(ranking.matrix) <- outcomes
    ##
    ranking.matrix.common <- ranking.matrix
    ranking.matrix.random <- ranking.matrix
    ##
    if (missing(common))
      common <- FALSE
    if (missing(random))
      random <- FALSE
  }
  ##
  ## More than one element provided in '...' (netmeta or netrank objects)
  ##
  else {
    ##
    ## (1) Check arguments
    ##
    n.outcomes <- length(args)
    ##
    if (!missing(treatments))
      warning("Argument 'treatments' ignored (only needed if ",
              "first argument is a ranking matrix).",
              call. = FALSE)
    ##
    if (missing(outcomes)) {
      warning("Outcomes are labelled 'A' to '", LETTERS[n.outcomes],
              "' as argument 'outcomes' is missing.",
              call. = FALSE)
      outcomes <- LETTERS[1:n.outcomes]  
    }
    ##
    if (length(outcomes) != n.outcomes)
      stop("Number of outcomes provided by argument 'outcomes' ",
           "and number of rankings differ.",
           call. = FALSE)
    ##
    if (!missing.small.values) {
      if (length(small.values) != n.outcomes)
        stop("Number of values provided by argument 'small.values' ",
             "and number of rankings differ.",
             call. = FALSE)
    }
    else
      small.values <- rep("", n.outcomes)
    ##
    for (i in seq_along(args)) {
      if (inherits(args[[i]], "netmeta"))
        any.netmeta <- TRUE
      else if (inherits(args[[i]], "netrank"))
        any.netrank <- TRUE
      else
        stop("All elements of argument '...' must be of class ",
             "'netmeta' or 'netrank'.",
             call. = FALSE)
      ##
      args[[i]] <- updateversion(args[[i]])
      ##
      if (!missing.small.values)
        small.values[i] <- setsv(small.values[i])
      else {
        small.values[i] <-
          ifelse(is.null(args[[i]]$small.values),
                 "desirable", args[[i]]$small.values)
      }
      ##
      if (length(small.values) != n.outcomes)
        stop("Length of argument 'small.values' must be equal to ",
             "number of outcomes.",
             call. = FALSE)
    }
    ##
    ## (2) Extract P-Scores
    ##
    ranking.list.common <- ranking.list.random <- list()
    commons <- randoms <- rep_len(NA, length(args))
    ##
    for (i in seq_along(args)) {
      ##
      args.i <- args[[i]]
      ##
      if (inherits(args.i, "netmeta")) {
        ranking.list.common[[i]] <- netrank(args.i,
                                            small.values =
                                              small.values[i])$ranking.common
        ranking.list.random[[i]] <- netrank(args.i,
                                            small.values =
                                              small.values[i])$ranking.random
        commons[i]  <- args.i$common
        randoms[i] <- args.i$random
      }
      else if (inherits(args.i, "netrank")) {
        ranking.list.common[[i]] <- args.i$ranking.common
        ranking.list.random[[i]] <- args.i$ranking.random
        commons[i]  <- args.i$x$common
        randoms[i] <- args.i$x$random
      }
    }
    
    
    ranking.treatments <- lapply(ranking.list.common, names)
    n.treatments <- unlist(lapply(ranking.treatments, length))
    ##
    if (length(unique(n.treatments)) != 1) {
      sel.max <- seq_along(n.treatments)[n.treatments == max(n.treatments)][1]
      ##
      treatments <- ranking.treatments[[sel.max]]
      ##
      ## Act on rankings with missing treatments
      ##
      for (j in seq_along(n.treatments)[n.treatments < max(n.treatments)]) {
        treatments.j <- ranking.treatments[[j]]
        missing.j <- !(treatments %in% treatments.j)
        ##
        if (any(treatments.j != treatments[!missing.j]))
          stop("Treatment names of all rankings must be in same order.",
               call. = FALSE)
        ##
        ranking.j <- ranking.list.common[[j]]
        ranking.list.common[[j]] <- ranking.list.common[[sel.max]]
        ranking.list.common[[j]][treatments[!missing.j]] <- ranking.j
        ranking.list.common[[j]][treatments[missing.j]]  <- NA
        ##
        ranking.j <- ranking.list.random[[j]]
        ranking.list.random[[j]] <- ranking.list.random[[sel.max]]
        ranking.list.random[[j]][treatments[!missing.j]] <- ranking.j
        ranking.list.random[[j]][treatments[missing.j]]  <- NA
        ##
        ranking.treatments[[j]] <- treatments
      }
    }
    else
      treatments <- ranking.treatments[[1]]
    ##
    for (i in seq_along(ranking.treatments))
      if (any(ranking.treatments[[i]] != treatments)) {
        if (all(sort(ranking.treatments[[i]]) == sort(treatments)))
          stop("Different order of treatments provided:\n ",
               paste(paste0("'", treatments, "'"), collapse = " - "),
               "\n ",
               paste(paste0("'", ranking.treatments[[i]], "'"),
                     collapse = " - "),
               call. = FALSE)
        else
          stop("Different treatments provided:\n ",
               paste(paste0("'", treatments, "'"), collapse = " - "),
               "\n ",
               paste(paste0("'", ranking.treatments[[i]], "'"),
                     collapse = " - "),
               call. = FALSE)
      }
    ##
    ranking.matrix.common <- matrix(unlist(ranking.list.common,
                                          use.names = FALSE),
                                   ncol = length(outcomes), byrow = FALSE)
    ranking.matrix.random <- matrix(unlist(ranking.list.random,
                                           use.names = FALSE),
                                    ncol = length(outcomes), byrow = FALSE)
    rownames(ranking.matrix.common) <- treatments
    rownames(ranking.matrix.random) <- treatments
    colnames(ranking.matrix.common) <- outcomes
    colnames(ranking.matrix.random) <- outcomes
    ##
    text.netmeta.netrank <- if (any.netmeta) "netmeta" else ""
    text.netmeta.netrank <- if (any.netmeta & any.netrank) "netmeta and netrank"
    text.netmeta.netrank <- if (!any.netmeta & any.netrank) "netrank"
    text.netmeta.netrank <- paste(text.netmeta.netrank, "objects.")
    ##
    if (missing(common)) {
      common <- unique(commons)
      if (length(common) != 1) {
        warning("Argument 'common' set to TRUE as different values are ",
                "available in ", text.netmeta.netrank,
                call. = FALSE)
        common <- TRUE
      }
    }
    ##
    if (missing(random)) {
      random <- unique(randoms)
      if (length(random) != 1) {
        warning("Argument 'random' set to TRUE as different values are ",
                "available in ", text.netmeta.netrank,
                call. = FALSE)
        random <- TRUE
      }
    }
  } 
  
  
  n <- nrow(ranking.matrix.common)
  o <- ncol(ranking.matrix.common)
  
  
  Pos.common <- M.common <- matrix(0, nrow = n, ncol = n)
  ##
  rownames(Pos.common) <- colnames(Pos.common) <-
    rownames(ranking.matrix.common)
  rownames(M.common) <- colnames(M.common) <- rownames(ranking.matrix.common)
  ##
  Pos.random <- M.random <- matrix(0, nrow = n, ncol = n)
  ##
  rownames(Pos.random) <- colnames(Pos.random) <-
    rownames(ranking.matrix.random)
  rownames(M.random) <- colnames(M.random) <- rownames(ranking.matrix.random)
  
  
  ## Pos[i, j] counts how many rankings judge treatment i superior to
  ## treatment j (number between 0 and o)
  ##
  for (i in 1:n)
    for (j in 1:n)
      if (i != j)
        for (k in 1:o)
          if (!is.na(ranking.matrix.common[i, k]) &
              !is.na(ranking.matrix.common[j, k]) &
              ranking.matrix.common[i, k] >= ranking.matrix.common[j, k])
            Pos.common[i, j] <- Pos.common[i, j] + 1
  ##
  for (i in 1:n)
    for (j in 1:n)
      if (i != j)
        for (k in 1:o)
          if (!is.na(ranking.matrix.random[i, k]) &
              !is.na(ranking.matrix.random[j, k]) &
              ranking.matrix.random[i, k] >= ranking.matrix.random[j, k])
            Pos.random[i, j] <- Pos.random[i, j] + 1
  
  
  ## Calculate "full" Hasse matrix M
  ##
  M.common[Pos.common == o] <- 1
  M.random[Pos.random == o] <- 1
  
  
  ## Matrix with information about partial ordering
  ##
  PO.common <- M.common[1:2, ]
  PO.common[1, ] <- rowSums(M.common)
  PO.common[2, ] <- colSums(M.common)
  rownames(PO.common) <- c("inferior", "superior")
  ##
  PO.random <- M.random[1:2, ]
  PO.random[1, ] <- rowSums(M.random)
  PO.random[2, ] <- colSums(M.random)
  rownames(PO.random) <- c("inferior", "superior")
  
  
  ## Skipping each direct path where a path of length 2 is found
  ##
  M0.common <- M.common - sign(M.common %*% M.common)
  ##
  M0.random <- M.random - sign(M.random %*% M.random)
  
  
  if (missing.small.values | !any.netmeta | length(args) == 1)
    small.values <- NULL
  
  
  res <- list(P.common = ranking.matrix.common,
              M0.common = M0.common,
              M.common = M.common,
              O.common = PO.common,
              P.random = ranking.matrix.random,
              M0.random = M0.random,
              M.random = M.random,
              O.random = PO.random,
              small.values = small.values,
              common = common,
              random = random,
              call = match.call(),
              version = packageDescription("netmeta")$Version)
  ##
  class(res) <- "netposet"
  
  res
}


#' @rdname netposet
#' @method print netposet
#' @export

print.netposet <- function(x,
                           pooled = ifelse(x$random, "random", "common"),
                           ...) {
  
  chkclass(x, "netposet")
  x <- updateversion(x)
  
  
  pooled <- setchar(pooled, c("common", "random", "fixed"))
  pooled[pooled == "fixed"] <- "common"
  
  if (pooled == "common") {
    if (!(!x$common & !x$random))
      cat("Common effects model\n")
    print(x$M0.common, ...)
    if (x$random)
      cat("\n")
  }
  ##
  else {
    cat("Random effects model\n")
    print(x$M0.random, ...)
  }
  
  ## diag(M0) <- NA
  ## M0[is.na(M0)] <- "."
  ## prmatrix(M0, quote = FALSE, right = TRUE)
  
  invisible(NULL)
}
