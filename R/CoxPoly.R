#' This is a function to implement polynomial Cox regression for latency analysis
#' @import survival
#' @param data A data frame containing the data
#' @param time_start start time
#' @param time_end end time
#' @param status survival status. 1 for death and 0 for censored
#' @param exposure Name of the time-varying exposure variable. From the most recent to the furthest in time.
#' @param degree degree of the polynomial
#' @param latency prespecified latency
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
                    adjusted_model = NULL){
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

  regression_data <- data.frame(data[,c(time_start, time_end, status,adjusted_variable)], Sum, Sum_k)

  names(regression_data)[  names(regression_data) == "Sum"] <- "X0"

  # run Cox regression.
  model <- paste0("Surv(", time_start, ",", time_end, ",", status, ")")
  # include the cumulative exposure name
  model <- paste0(model, "~", "X0", "+", paste0(colnames(Sum_k), collapse = "+"))

  # include the adjusted variables if given
  if (length(adjusted_model) > 0) {model <- paste0(model, "+", adjusted_model)}

  fit <- coxph(as.formula(model),
               data = regression_data,
               control = coxph.control(timefix = FALSE))
  return(fit)
}

#' This is a function to extract log HR based on the resuls from CoxPoly
#' @param fit A model from CoxPoly
#' @param lag A value of lag time
#' @return The log HR estimated at that lag time
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
extract_CoxPoly <- function(fit, lag) {
  # extract the coefficients
  coef <- fit$coefficients
  coef <- coef[stringr::str_detect(names(coef), "\\bX\\d+\\b")]
  K <- length(coef) - 1

  log_HR <- coef["X0"] + as.matrix(coef[2:length(coef)] %*% lag^(1:K))
  return(log_HR)
}






