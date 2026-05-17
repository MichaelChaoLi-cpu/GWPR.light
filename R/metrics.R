#' @title Metrics Module
#' @description Functions for computing evaluation metrics for linear and
#'   logistic panel regression models.
#' @name metrics
NULL

#' Compute linear regression metrics
#'
#' Computes R2, MSE, RMSE, and MAE for a linear model.
#'
#' @param y Numeric vector of observed values.
#' @param y_hat Numeric vector of predicted values. May contain NA values,
#'   which are excluded from all metric calculations.
#'
#' @return A named list with elements: \code{R2}, \code{MSE}, \code{RMSE},
#'   \code{MAE}. Returns \code{NA_real_} for each metric when no complete
#'   observations remain after removing NAs.
#'
#' @examples
#' y <- c(1, 2, 3, 4, 5)
#' y_hat <- c(1, 2, 3, 4, 5)
#' compute_linear_metrics(y, y_hat)
#'
#' @export
compute_linear_metrics <- function(y, y_hat) {
  keep <- !is.na(y_hat) & !is.na(y)
  y     <- y[keep]
  y_hat <- y_hat[keep]

  if (length(y) == 0L) {
    return(list(R2 = NA_real_, MSE = NA_real_, RMSE = NA_real_, MAE = NA_real_))
  }

  ss_res <- sum((y - y_hat)^2)
  ss_tot <- sum((y - mean(y))^2)
  R2   <- if (ss_tot == 0) NA_real_ else 1 - ss_res / ss_tot
  MSE  <- mean((y - y_hat)^2)
  RMSE <- sqrt(MSE)
  MAE  <- mean(abs(y - y_hat))

  list(R2 = R2, MSE = MSE, RMSE = RMSE, MAE = MAE)
}

#' Compute logistic regression metrics
#'
#' Computes log_loss, accuracy, precision, recall, and f1_score for a binary
#' classification model.
#'
#' @param y Integer or numeric vector of true binary labels (0/1).
#' @param prob Numeric vector of predicted probabilities for class 1.
#' @param threshold Numeric scalar; classification threshold (default 0.5).
#'
#' @return A named list with elements: \code{log_loss}, \code{accuracy},
#'   \code{precision}, \code{recall}, \code{f1_score}. \code{precision} and
#'   \code{recall} are \code{NA_real_} when their denominators are zero.
#'
#' @examples
#' y <- c(1, 0, 1, 0)
#' prob <- c(0.9, 0.1, 0.8, 0.2)
#' compute_logistic_metrics(y, prob)
#'
#' @export
compute_logistic_metrics <- function(y, prob, threshold = 0.5) {
  ll         <- safe_log_loss(y, prob)
  class_pred <- as.integer(prob >= threshold)
  counts     <- safe_confusion_counts(y, class_pred)

  TP <- counts$TP
  FP <- counts$FP
  FN <- counts$FN
  TN <- counts$TN
  n  <- TP + FP + FN + TN

  accuracy  <- if (n == 0L) NA_real_ else (TP + TN) / n
  precision <- counts$precision
  recall    <- counts$recall

  f1_score <- if (is.na(precision) || is.na(recall)) {
    NA_real_
  } else {
    denom <- precision + recall
    if (denom == 0) NA_real_ else 2 * precision * recall / denom
  }

  list(
    log_loss  = ll,
    accuracy  = accuracy,
    precision = precision,
    recall    = recall,
    f1_score  = f1_score
  )
}

#' Numerically safe log loss
#'
#' Clips predicted probabilities to \code{[eps, 1 - eps]} before computing
#' binary cross-entropy, preventing infinite log-loss values.
#'
#' @param y Numeric vector of true binary labels (0/1).
#' @param prob Numeric vector of predicted probabilities.
#' @param eps Numeric scalar; clipping bound (default \code{1e-15}).
#'
#' @return A single numeric value: the mean binary cross-entropy.
#'
#' @examples
#' safe_log_loss(c(1, 0, 1), c(0.99, 0.01, 0.5))
#'
#' @export
safe_log_loss <- function(y, prob, eps = 1e-15) {
  prob <- pmax(pmin(prob, 1 - eps), eps)
  -mean(y * log(prob) + (1 - y) * log(1 - prob))
}

#' Compute confusion matrix counts
#'
#' Returns TP, FP, FN, TN and derived precision and recall. Issues a
#' \code{warning()} and returns \code{NA_real_} when precision or recall
#' denominators are zero.
#'
#' @param y Integer or numeric vector of true binary labels (0/1).
#' @param class_pred Integer or numeric vector of predicted binary labels (0/1).
#'
#' @return A named list: \code{TP}, \code{FP}, \code{FN}, \code{TN},
#'   \code{precision}, \code{recall}.
#'
#' @examples
#' safe_confusion_counts(c(1, 0, 1, 0), c(1, 1, 0, 0))
#'
#' @export
safe_confusion_counts <- function(y, class_pred) {
  TP <- sum(y == 1L & class_pred == 1L)
  FP <- sum(y == 0L & class_pred == 1L)
  FN <- sum(y == 1L & class_pred == 0L)
  TN <- sum(y == 0L & class_pred == 0L)

  precision_denom <- TP + FP
  recall_denom    <- TP + FN

  if (precision_denom == 0L) {
    warning("precision denominator (TP + FP) is zero; returning NA_real_")
    precision <- NA_real_
  } else {
    precision <- TP / precision_denom
  }

  if (recall_denom == 0L) {
    warning("recall denominator (TP + FN) is zero; returning NA_real_")
    recall <- NA_real_
  } else {
    recall <- TP / recall_denom
  }

  list(TP = TP, FP = FP, FN = FN, TN = TN,
       precision = precision, recall = recall)
}
