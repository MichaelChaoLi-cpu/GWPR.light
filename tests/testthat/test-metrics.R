test_that("compute_linear_metrics: perfect predictions give R2=1, MSE=0, RMSE=0, MAE=0", {
  y     <- c(1, 2, 3, 4, 5)
  y_hat <- c(1, 2, 3, 4, 5)
  m <- compute_linear_metrics(y, y_hat)
  expect_equal(m$R2,   1)
  expect_equal(m$MSE,  0)
  expect_equal(m$RMSE, 0)
  expect_equal(m$MAE,  0)
})

test_that("compute_linear_metrics: all-wrong predictions give correct metrics", {
  y     <- c(0, 0, 0, 0)
  y_hat <- c(1, 1, 1, 1)
  m <- compute_linear_metrics(y, y_hat)
  expect_equal(m$MSE,  1)
  expect_equal(m$RMSE, 1)
  expect_equal(m$MAE,  1)
  # ss_tot == 0 because all y equal => R2 is NA
  expect_true(is.na(m$R2))
})

test_that("compute_linear_metrics: known non-trivial values", {
  y     <- c(3, 3, 3, 3)     # all same
  y_hat <- c(2, 4, 2, 4)
  m <- compute_linear_metrics(y, y_hat)
  expect_equal(m$MSE,  1)
  expect_equal(m$RMSE, 1)
  expect_equal(m$MAE,  1)
})

test_that("compute_linear_metrics: NA in y_hat are dropped", {
  y     <- c(1, 2, 3, 4)
  y_hat <- c(1, 2, NA, 4)
  m <- compute_linear_metrics(y, y_hat)
  # Only indices 1,2,4 are used: errors = 0,0,0 => perfect
  expect_equal(m$MSE, 0)
  expect_equal(m$MAE, 0)
})

test_that("compute_linear_metrics: all NA y_hat returns NA metrics", {
  y     <- c(1, 2, 3)
  y_hat <- c(NA, NA, NA)
  m <- compute_linear_metrics(y, y_hat)
  expect_true(is.na(m$R2))
  expect_true(is.na(m$MSE))
  expect_true(is.na(m$RMSE))
  expect_true(is.na(m$MAE))
})

# ---------------------------------------------------------------------------
# safe_log_loss
# ---------------------------------------------------------------------------

test_that("safe_log_loss: perfect predictions give log_loss near 0", {
  y    <- c(1, 0, 1, 0)
  prob <- c(1 - 1e-10, 1e-10, 1 - 1e-10, 1e-10)
  ll <- safe_log_loss(y, prob)
  expect_true(is.finite(ll))
  expect_lt(ll, 1e-8)
})

test_that("safe_log_loss: probabilities at 0/1 remain finite (no Inf)", {
  y    <- c(1, 0, 1, 0)
  prob <- c(1, 0, 1, 0)   # would produce Inf without clipping
  ll <- safe_log_loss(y, prob)
  expect_true(is.finite(ll))
})

test_that("safe_log_loss: probabilities near 0/1 remain finite", {
  y    <- c(1, 0)
  prob <- c(1e-300, 1 - 1e-300)
  ll <- safe_log_loss(y, prob)
  expect_true(is.finite(ll))
})

# ---------------------------------------------------------------------------
# safe_confusion_counts
# ---------------------------------------------------------------------------

test_that("safe_confusion_counts: basic correctness", {
  y          <- c(1, 1, 0, 0)
  class_pred <- c(1, 0, 1, 0)
  cc <- safe_confusion_counts(y, class_pred)
  expect_equal(cc$TP, 1L)
  expect_equal(cc$FP, 1L)
  expect_equal(cc$FN, 1L)
  expect_equal(cc$TN, 1L)
  expect_equal(cc$precision, 0.5)
  expect_equal(cc$recall,    0.5)
})

test_that("safe_confusion_counts: precision NA when no positive predictions", {
  y          <- c(1, 0, 1, 0)
  class_pred <- c(0, 0, 0, 0)   # TP + FP == 0
  expect_warning(
    cc <- safe_confusion_counts(y, class_pred),
    regexp = "precision denominator"
  )
  expect_true(is.na(cc$precision))
  expect_false(is.na(cc$recall))
})

test_that("safe_confusion_counts: recall NA when no positive ground-truth", {
  y          <- c(0, 0, 0, 0)   # TP + FN == 0
  class_pred <- c(1, 0, 1, 0)
  expect_warning(
    cc <- safe_confusion_counts(y, class_pred),
    regexp = "recall denominator"
  )
  expect_true(is.na(cc$recall))
  # precision is defined here (TP=0, FP=2)
  expect_false(is.na(cc$precision))
})

# ---------------------------------------------------------------------------
# compute_logistic_metrics
# ---------------------------------------------------------------------------

test_that("compute_logistic_metrics: perfect predictions give accuracy=1, log_loss~0", {
  y    <- c(1, 0, 1, 0)
  prob <- c(0.999, 0.001, 0.999, 0.001)
  m <- compute_logistic_metrics(y, prob)
  expect_true(is.finite(m$log_loss))
  expect_lt(m$log_loss, 0.01)
  expect_equal(m$accuracy, 1)
  expect_equal(m$precision, 1)
  expect_equal(m$recall, 1)
  expect_equal(m$f1_score, 1)
})

test_that("compute_logistic_metrics: all-wrong predictions", {
  y    <- c(1, 1, 0, 0)
  prob <- c(0.01, 0.01, 0.99, 0.99)   # everything wrong
  m <- compute_logistic_metrics(y, prob)
  expect_equal(m$accuracy, 0)
  # precision: TP=0, FP=2 => 0; recall: TP=0, FN=2 => 0
  expect_equal(m$precision, 0)
  expect_equal(m$recall, 0)
  expect_true(is.na(m$f1_score) || m$f1_score == 0)
})

test_that("compute_logistic_metrics: threshold change alters classification", {
  y    <- c(1, 1, 0, 0)
  prob <- c(0.6, 0.4, 0.6, 0.4)

  m_low  <- suppressWarnings(compute_logistic_metrics(y, prob, threshold = 0.3))
  m_high <- suppressWarnings(compute_logistic_metrics(y, prob, threshold = 0.7))

  # At threshold 0.3 all are predicted positive
  expect_equal(m_low$accuracy, 0.5)   # 2 correct out of 4

  # At threshold 0.7 all are predicted negative
  expect_equal(m_high$accuracy, 0.5)  # 2 correct out of 4

  # Metrics should differ between thresholds
  expect_false(identical(m_low$precision, m_high$precision) &&
                 identical(m_low$recall,    m_high$recall))
})

test_that("compute_logistic_metrics: precision NA does not error", {
  y    <- c(1, 1, 1, 1)   # no negatives
  prob <- c(0.4, 0.4, 0.4, 0.4)  # all predicted negative at default threshold
  # TP + FP = 0 => precision NA
  expect_warning(
    m <- compute_logistic_metrics(y, prob),
    regexp = "precision"
  )
  expect_true(is.na(m$precision))
  expect_true(is.na(m$f1_score))
})

test_that("compute_logistic_metrics: recall NA does not error", {
  y    <- c(0, 0, 0, 0)   # no positives
  prob <- c(0.6, 0.6, 0.6, 0.6)  # all predicted positive
  # TP + FN = 0 => recall NA
  expect_warning(
    m <- compute_logistic_metrics(y, prob),
    regexp = "recall"
  )
  expect_true(is.na(m$recall))
  expect_true(is.na(m$f1_score))
})

test_that("compute_logistic_metrics: f1_score NA when both precision and recall are NA", {
  # All same class, all predicted wrong class => both NA
  y    <- c(0, 0, 0, 0)
  prob <- c(0.4, 0.4, 0.4, 0.4)  # all predicted 0 (negative)
  # TP=0, FP=0 => precision NA; TP=0, FN=0 => recall NA
  suppressWarnings(m <- compute_logistic_metrics(y, prob))
  expect_true(is.na(m$f1_score))
})
