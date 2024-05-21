#' This function implments GridSearch of natural cubic splines in Cox regression for latency analysis
#' This is GreadSearching by knots
#' @import survival
#' @import Hmisc
#' @import dplyr
#' @import stats
#' @param data A data frame containing the data
#' @param time_start start time
#' @param time_end end time
#' @param status survival status. 1 for death and 0 for censored
#' @param exposure Name of the time-varying exposure variable. From the most recent to the furthest in time.
#' @param  knots_number A vector of potential number of knots
#' @param latency pre-specified latency
#' @return A list. The first element is the best knots with the lowest AIC. The second element is a Cox model
#' @export
#' @examples
#' time_start <- "age_start"
#' time_end <- "age_end"
#' status <- "failure"
#' exposure <- paste0("lag", 0:19)
#' knots_number <- c(3,4,5)
#' latency <- 20
#' fit <- CoxKnotsearch(sim_data, time_start, time_end, status, exposure, knots_number, latency)
CoxKnotsearch <- function(data, time_start, time_end, status, exposure,
                          knots_number, latency, adjusted_variable = NULL, adjusted_model = NULL) {
  X <- data[,exposure]
  X <- as.matrix(X)

  v <- 0:(latency - 1)
  mini_AIC <- 99999999
  fit_best <- NA
  knots_best <- NA
  # here we may put multiple possible number of knots as  knots_number
  for (j in knots_number){
    # generate the B matrix
    B <- matrix(NA, nrow = latency, ncol = j )
    # the first column of B matrix is the intercept, and thus, is 1
    B[,1] <- 1

    # Make a list of potential knots
    knot_list <- vector(mode = "list", length = j)
    for (i in 1:j){
      cut_start <- quantile(v, 0.1) %>% round()
      cut_end <- quantile(v, 0.9) %>% round()
      if (i == 1) {
        knot_list[[i]] <- min(v):cut_start
      } else if (i == j){
        knot_list[[i]] <- (cut_end +1):max(v)
      } else {
        knot_list[[i]] <- (round(quantile(cut_start:cut_end, (i-2)/(j-2)))+1) :
          round(quantile(cut_start:cut_end, (i-1)/(j-2)))
      }
    }

    knot_list <- expand.grid(knot_list)
    for (i in 1:dim(knot_list)[1]){
      knots <- as.numeric(knot_list[i,])
      spline <- rcspline.eval(0:(latency-1), knots = knots, inclx = TRUE)
      B[, 2:(j)] <- as.matrix(spline)

      Sum_mat <- matrix(NA, nrow = dim(X)[1], ncol = j)
      Sum_mat <- X %*% B
      Sum_mat <- as.data.frame(Sum_mat)
      # then make a data frame to fit the regression
      regression_data <- cbind(data[,c(time_start, time_end, status,adjusted_variable)],  Sum_mat )

      model <- paste0("Surv(", time_start, ",", time_end, ",", status, ")")
      # include the cumulative exposure name
      model <- paste0(model, "~", paste0(colnames(Sum_mat), collapse = "+"))

      # include the adjusted variables if given
      if (length(adjusted_model) > 0) {model <- paste0(model, "+", adjusted_model)}


      fit <- coxph(as.formula(model),
                   data = regression_data,
                   control = coxph.control(timefix = FALSE))
      AIC <- extractAIC(fit)[2]
      if (AIC < mini_AIC) {
        mini_AIC <- AIC
        fit_best <- fit
        knots_best <- knots
      }
    }
  }
  list(knots_best, coef(fit_best)) %>% return()
}


#' This is a function to extract log HR based on the resuls from CoxKnotsearch
#' @param fit A model from CoxKnotsearch
#' @param lag A value of lag time
#' @param latency prespecified latency
#' @return The log HR estimated at that lag time
#' @export
#' @examples
#' time_start <- "age_start"
#' time_end <- "age_end"
#' status <- "failure"
#' exposure <- paste0("lag", 0:19)
#' knots_number <- c(3,4,5)
#' latency <- 20
#' fit <- CoxKnotsearch(sim_data, time_start, time_end, status, exposure, knots_number, latency)
#' extract_CoxKnotsearch(fit, 10, latency = latency)
extract_CoxKnotsearch  <- function(fit, lag, latency) {

  #extract the knots
  knots <- fit[[1]]

  # extract the coefficients
  coef <- fit[[2]]

  # make a function
  B <- matrix(NA, nrow = latency, ncol = length(knots) )
  # the first column of B matrix is the intercept, and thus, is 1
  B[,1] <- 1
  spline <- rcspline.eval(0:(latency-1), inclx = TRUE, knots = knots)
  B[, 2:length(knots)] <- as.matrix(spline)

  # the second column of B is the lag time
  B_lag <- B[which(B[,2] == lag), ]

  log_HR <- t(coef) %*% B_lag

  return(log_HR)
}

