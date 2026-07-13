#' Signal Extraction for DIVASCHOIR Categorical Blocks
#'
#' Extracts categorical signal spaces from DIVASCHOIR/MCA-centered categorical
#' data blocks. Raw categorical blocks are first transformed with
#' `MatCenterDIVASCHOIR`; numeric matrices are treated as already centered.
#'
#' @param categorical_datablock A list of categorical blocks or MCA-centered
#'   numeric matrices. Blocks should have matched samples.
#' @param ranks Optional vector of ranks to keep for each block. If `NULL`, the
#'   numerical rank is estimated from the singular values.
#' @param tol Numerical tolerance used for rank estimation when `ranks = NULL`.
#'
#' @return A list containing:
#'   \describe{
#'     \item{VBars}{List of sample-space right singular vector matrices.}
#'     \item{UBars}{List of category-space left singular vector matrices.}
#'     \item{rBars}{Vector of retained ranks for each categorical block.}
#'     \item{singVals}{List of singular values for each categorical block.}
#'     \item{mcaBlocks}{List of MCA-centered matrices decomposed by SVD.}
#'   }
#'
#' @export
DIVASCHOIRSignalExtract <- function(categorical_datablock, ranks = NULL, tol = sqrt(.Machine$double.eps)) {
  categorical_datablock <- if (is.null(categorical_datablock)) list() else categorical_datablock

  nb <- length(categorical_datablock)
  dataname <- names(categorical_datablock)
  if (is.null(dataname) && nb > 0) {
    dataname <- paste0("CategoricalDatablock", 1:nb)
  }

  VBars <- vector("list", nb)
  UBars <- vector("list", nb)
  singVals <- vector("list", nb)
  rBars <- numeric(nb)
  mcaBlocks <- vector("list", nb)

  for (ib in seq_len(nb)) {
    cat(sprintf("Categorical signal estimation for %s\n", dataname[ib]))

    X <- categorical_datablock[[ib]]
    if (!is.matrix(X) || !is.numeric(X)) {
      X <- MatCenterDIVASCHOIR(X)
    }

    svdResult <- La.svd(X)
    svals <- svdResult$d

    if (is.null(ranks)) {
      rank_tol <- tol * max(dim(X)) * max(svals)
      r <- sum(svals > rank_tol)
    } else {
      r <- ranks[ib]
    }
    r <- min(r, length(svals))

    UBars[[ib]] <- svdResult$u[, seq_len(r), drop = FALSE]
    VBars[[ib]] <- t(svdResult$vt)[, seq_len(r), drop = FALSE]
    singVals[[ib]] <- svals
    rBars[ib] <- r
    mcaBlocks[[ib]] <- X
  }

  return(list(
    VBars = VBars,
    UBars = UBars,
    rBars = rBars,
    singVals = singVals,
    mcaBlocks = mcaBlocks
  ))
}
