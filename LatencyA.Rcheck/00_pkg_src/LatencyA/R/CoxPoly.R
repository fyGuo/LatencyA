#' This is a function to implement polynomial Cox regression for latency analysis
#' @importFrom survival coxph coxph.control
#' @param data A data frame containing the data
#' @param time_start start time
#' @param time_end end time
#' @param status survival status. 1 for death and 0 for censored
#' @param exposure Name of the time-varying exposure variable. From the most recent to the furthest in time.
#' @param degree degree of the polynomial
#' @param latency prespecified latency
#' @param adjusted_variable A vector of adjusted variables
#' @param adjusted_model A vector of adjusted models
#' @param id A string to specify the name of the id column
#' @return A Cox model object
#' @export
#' @examples
#' time_start <- "age_start"
#' time_end <- "age_end"
#' status <- "failure"
#' exposure <- paste0("lag", 0:19)
#' degree <- 3
#' latency <- 20
#' fit <- CoxPoly(sim_data, time_start, time_end, status, exposure, degree, latency)
#' fit

CoxPoly <- function(data, time_start, time_end, status, exposure, degree, latency, adjusted_variable = NULL,
                    adjusted_model = NULL, id = "id"){

  # rename the id column to "id"
  names(data)[names(data) == id] <- "id"

  # extract time-varying exposures
  X <- data[,exposure]
  X <- as.matrix(X)

  # create a vector from 0 to latency - 1 to denote time
  v <- 0:(latency - 1)

  # create a matrix for the transformation information
  Sum <- rowSums(X)
  Sum_k <- matrix(0, nrow = dim(X)[1], ncol = degree)
  for (k in 1:degree) {
    Sum_k[,k] <- X %*% (v^(k))
  }
  Sum_k <- data.frame(Sum_k)

  regression_data <- data.frame(data[,c(time_start, time_end, status,adjusted_variable, "id")], Sum, Sum_k)

  names(regression_data)[  names(regression_data) == "Sum"] <- "X0"

  # run Cox regression.
  model <- paste0("Surv(", time_start, ",", time_end, ",", status, ")")
  # include the cumulative exposure name
  model <- paste0(model, "~", "X0", "+", paste0(colnames(Sum_k), collapse = "+"))

  # include the adjusted variables if given
  if (length(adjusted_model) > 0) {model <- paste0(model, "+", adjusted_model)}

  fit <- coxph(as.formula(model),
               data = regression_data,
               control = coxph.control(timefix = FALSE),
               cluster = id)
  return(fit)
}

#' This is a function to extract log HR based on the resuls from CoxPoly
#' @param fit A model from CoxPoly
#' @param lag A lag time, or a vector of lag times
#' @return A data frame with one row per lag, giving the estimated log HR and its
#'   variance at each lag time
#' @export
#' @examples
#' time_start <- "age_start"
#' time_end <- "age_end"
#' status <- "failure"
#' exposure <- paste0("lag", 0:19)
#' degree <- 3
#' latency <- 20
#' fit <- CoxPoly(sim_data, time_start, time_end, status, exposure, degree, latency)
#' extract_CoxPoly(fit, 10)
#' extract_CoxPoly(fit, c(0, 5, 10))
extract_CoxPoly <- function(fit, lag) {
  # extract the coefficients
  coef <- coef(fit)

  # extract the variance-covariance matrix
  vcov <- vcov(fit)

  coef <- coef[stringr::str_detect(names(coef), "\\bX\\d+\\b")]


  vcov <- vcov[names(coef),
               names(coef)]

  K <- length(coef) - 1

  # one polynomial-basis row per requested lag: row i is c(1, lag[i], ..., lag[i]^K).
  # outer() keeps it a matrix even for a single lag, so the scalar case is just
  # the 1-row case.
  B_lag <- outer(lag, 0:K, "^")

  # log HR at each lag, and the per-lag variance from the diagonal of
  # B_lag %*% vcov %*% t(B_lag), computed without forming the full product
  log_HR <- as.vector(B_lag %*% coef)
  log_HR_var <- rowSums((B_lag %*% vcov) * B_lag)

  return(data.frame(lag = lag, log_HR = log_HR, log_HR_var = log_HR_var))
}






