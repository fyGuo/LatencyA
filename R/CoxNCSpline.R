#' This function implments natural cubic splines in Cox regression for latency analysis
#' @import survival
#' @import Hmisc
#' @import dplyr
#' @import stats
#' @param data A data frame containing the data
#' @param time_start start time
#' @param time_end end time
#' @param status survival status. 1 for death and 0 for censored
#' @param exposure Name of the time-varying exposure variable. From the most recent to the furthest in time.
#' @param knots_number Pre-specified number of knots
#' @param latency prespecified latency
#' @param adjusted_variable A vector of adjusted variables
#' @param adjusted_model A vector of adjusted models
#' @param id A string to specify the name of the id column
#' @return A list. The first element is the best knots with the lowest AIC. The second element is a Cox model
#' @export
#' @examples
#' time_start <- "age_start"
#' time_end <- "age_end"
#' status <- "failure"
#' exposure <- paste0("lag", 0:19)
#' knots_number <- 3
#' latency <- 20
#' fit <- CoxNCSpline(sim_data, time_start, time_end, status, exposure, knots_number, latency)



CoxNCSpline <- function(data, time_start, time_end, status, exposure,knots_number, latency,
                        adjusted_variable = NULL, adjusted_model = NULL, id = "id"){
  # rename the id column to "id"
  names(data)[names(data) == id] <- "id"

  # extract time-varying exposures
  X <- data[,exposure]
  X <- as.matrix(X)

  # create a vector from 0 to latency - 1 to denote time
  v <- 0:(latency - 1)


  # generate the B matrix
   B <- matrix(NA, nrow = latency, ncol = knots_number)

   # the first column of B matrix is the intercept, and thus, is 1
   B[,1] <- 1

   # automatically set knots based on the quntiles
   knots <- c(1, quantile(1:(latency-2), 1:(knots_number-2) * 1/(knots_number-1)), (latency-2))
   spline <- Hmisc::rcspline.eval(v, knots = knots, inclx = TRUE)

   # set the other columns of the B matrix
   B[, 2:knots_number] <- as.matrix(spline)

   # record the knots and output it
   knots <- attributes(spline)$knots

   # make a summation matrix
   Sum_mat <- matrix(NA, nrow = dim(X)[1], ncol = knots_number)
   Sum_mat <- X %*% B
   Sum_mat <- as.data.frame(Sum_mat)


   # specify the full model
   # then make a data frame to fit the regression
   regression_data <- cbind(data[,c(time_start, time_end, status,adjusted_variable, "id")],  Sum_mat )

   model <- paste0("Surv(", time_start, ",", time_end, ",", status, ")")
   # include the cumulative exposure name
   model <- paste0(model, "~", paste0(colnames(Sum_mat), collapse = "+"))

   # include the adjusted variables if given
   if (length(adjusted_model) > 0) {model <- paste0(model, "+", adjusted_model)}
  # run Cox regression.
  fit <- coxph(as.formula(model),
               data = regression_data,
               control = coxph.control(timefix = FALSE),
               cluster = id)
  list(knots, fit) %>% return()
}


#' This is a function to extract log HR based on the resuls from CoxNCSpline
#' @param fit A model from CoxNCSpline
#' @param lag A value of lag time
#' @param latency prespecified latency
#' @return The log HR estimated at that lag time
#' @export
#' @examples
#' time_start <- "age_start"
#' time_end <- "age_end"
#' status <- "failure"
#' exposure <- paste0("lag", 0:19)
#' knots_number <- 3
#' latency <- 20
#' fit <- CoxNCSpline(sim_data, time_start, time_end, status, exposure, knots_number, latency)
#' extract_CoxNCSpline(fit, 10, latency = latency)
extract_CoxNCSpline  <- function(fit, lag, latency) {

  #extract the knots
  knots <- fit[[1]]

  # extract the coefficients
  coef <- coef(fit[[2]])

  # extract the variance-covariance matrix
  vcov <- vcov(fit[[2]])

  # first we want to extract the coefficient for the exposure
  # we extract the coefficients with "VX", where X stands for numbers
  coef <- coef[stringr::str_detect(names(coef), "\\bV\\d+\\b")]

  vcov <- vcov[names(coef),
               names(coef)]
  # make a function
  B <- matrix(NA, nrow = latency, ncol = length(knots) )
  # the first column of B matrix is the intercept, and thus, is 1
  B[,1] <- 1
  spline <- rcspline.eval(0:(latency - 1), inclx = TRUE, knots = knots)
  B[, 2:length(knots)] <- as.matrix(spline)

  # the second column of B is the lag time
  B_lag <- B[which(B[,2] == lag), ]

  log_HR_estimate <- t(coef) %*% B_lag
  log_HR_var <- t(B_lag) %*% vcov %*% B_lag

  data.frame(log_HR = log_HR_estimate, log_HR_var = log_HR_var) %>% return()
}

