#' Construct MCA-Style SVD Input for Categorical DIVAS Blocks
#'
#' `MatCenterDIVASCHOIR` converts a categorical data block into the centered
#' and scaled indicator matrix used for MCA-style singular value decomposition.
#'
#' For a block with `n` samples and total one-hot dimension `J`, this constructs
#' the matrix
#'
#' \deqn{\tilde{M}_k^T = D_{\hat{p}_k}^{-1/2}
#' ((X_k^{obs})^T - \hat{p}_k 1_n^T)}
#'
#' where `X_obs` is the `n x J` stacked one-hot matrix,
#' `p_hat = n^{-1} t(X_obs) 1_n`, and `D_p_hat = diag(p_hat)`.
#'
#' @param block A data.frame, matrix, or vector containing categorical variables
#'   measured across samples. Rows are samples and columns are categorical
#'   variables. A vector is treated as one categorical variable.
#' @param return_details Logical. If `FALSE`, return only the MCA-style SVD input
#'   matrix with categories/features in rows and samples in columns. If `TRUE`,
#'   return a list containing the SVD input, the one-hot matrix, and category
#'   masses.
#'
#' @return If `return_details = FALSE`, a numeric matrix with categories/features
#'   in rows and samples in columns. If `return_details = TRUE`, a list with:
#'   \describe{
#'     \item{svd_input}{The centered and scaled MCA-style SVD input matrix.}
#'     \item{indicator}{The `n x J` one-hot indicator matrix `X_obs`.}
#'     \item{p_hat}{The category masses used for centering and scaling.}
#'   }
#'
#' @export
MatCenterDIVASCHOIR <- function(block, return_details = FALSE) {
  if (is.null(block)) {
    stop("Input block cannot be NULL.")
  }

  if (is.vector(block) && !is.list(block)) {
    block <- data.frame(variable = block)
  } else {
    block <- as.data.frame(block)
  }

  if (nrow(block) < 1 || ncol(block) < 1) {
    stop("Input block must contain at least one sample and one categorical variable.")
  }

  if (anyNA(block)) {
    stop("Input block contains NA values. Please impute or remove missing categorical values before calling MatCenterDIVASCHOIR.")
  }

  block[] <- lapply(block, function(x) {
    if (is.factor(x)) {
      x
    } else {
      factor(x)
    }
  })

  indicator_list <- Map(function(x, variable_name) {
    mm <- stats::model.matrix(~ x - 1)
    level_names <- sub("^x", "", colnames(mm))
    colnames(mm) <- paste0(variable_name, level_names)
    mm
  }, block, names(block))
  indicator <- do.call(cbind, indicator_list)
  p_hat <- colMeans(indicator)

  if (any(p_hat <= 0)) {
    stop("At least one category has zero mass after one-hot encoding.")
  }

  svd_input <- sweep(t(indicator), 1, p_hat, "-")
  svd_input <- sweep(svd_input, 1, sqrt(p_hat), "/")

  if (!is.null(rownames(block))) {
    colnames(svd_input) <- rownames(block)
  }

  if (return_details) {
    return(list(
      svd_input = svd_input,
      indicator = indicator,
      p_hat = p_hat
    ))
  }

  svd_input
}
