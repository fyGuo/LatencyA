#' This function implments Stepwise model selection by terms of natural cubic splines in Cox regression for latency analysis
#' This is GreadSearching by knots
#' @import survival
#' @import Hmisc
#' @import dplyr
#' @import stats
#' @importFrom MASS stepAIC
#' @param data A data frame containing the data
#' @param time_start start time
#' @param time_end end time
#' @param status survival status. 1 for death and 0 for censored
#' @param exposure Name of the time-varying exposure variable. From the most recent to the furthest in time.
#' @param  knots A vector of prespecified knots
#' @param latency prespecified latency
#' @return A list. The first element is the best knots with the lowest AIC. The second element is a Cox model
#' @export
#' @examples
#' time_start <- "age_start"
#' time_end <- "age_end"
#' status <- "failure"
#' exposure <- paste0("lag", 0:19)
#' knots <- 0:19
#' latency <- 20
#' fit <- CoxTermsearch(sim_data, time_start, time_end, status, exposure, knots = knots, latency)


CoxTermsearch <- function(data, time_start, time_end, status, exposure,
                          knots, latency) {
  X <- data[,exposure]
  X <- as.matrix(X)

  v <- 0:(latency - 1)
  # generate the B matrix
  B <- matrix(NA, nrow = latency, ncol = length(knots)  )
  # the first column of B matrix is the intercept, and thus, is 1
  B[,1] <- 1

  spline <- rcspline.eval(0:(latency-1), knots = knots, inclx = TRUE)
  B[, 2:dim(B)[2]] <- as.matrix(spline)

  # generate the sum matrix after transformation
  Sum_mat <- matrix(NA, nrow = dim(X)[1], ncol = dim(B)[2])
  Sum_mat <- X %*% B


  # specifiy the full model

  temp <- data.frame(time_start = data[,time_start],
                     time_end = data[,time_end],
                     status = data[,status],
                     Sum_mat)

  full_model <- coxph(Surv(temp$time_start, temp$time_end, temp$status)~Sum_mat,
                      control = coxph.control(timefix = FALSE))

  fit_best <- stepAIC(full_model, direction = "both", trace = FALSE)

  return(fit_best)
}

#' This is a function to extract log HR based on the resuls from CoxTermsearch
#' @param fit A model from CoxKnotsearch
#' @param lag A value of lag time
#' @param latency prespecified latency
#' @param knots A vector of prespecified knots
#' @return The log HR estimated at that lag time
#' @import stringr
#' @export
#' @examples
#' time_start <- "age_start"
#' time_end <- "age_end"
#' status <- "failure"
#' exposure <- paste0("lag", 0:19)
#' knots <- 0:19
#' latency <- 20
#' fit <- CoxTermsearch(sim_data, time_start, time_end, status, exposure, knots = knots, latency)
#' extract_CoxTermsearch(fit, 10, latency = latency, knot = knots)

extract_CoxTermsearch  <- function(fit, lag, latency, knots) {
  coef <- coef(fit)
  # replace NA with 0 which means the coefficient has no effects
  coef[!is.na(coef)] <- 0

  # if the selected model contains no predictors, it means that no effect and logHR is 0
  if (is.null(coef)){
    log_HR = 0
  } else {

    # check terms that is used in the final model
    s <- stringr::str_extract(names(coef), "\\d+") %>% as.numeric()

    # make the B matrix
    v <- 0:(latency - 1)
    # generate the B matrix
    B <- matrix(NA, nrow = latency, ncol = length(knots)  )
    # the first column of B matrix is the intercept, and thus, is 1
    B[,1] <- 1

    spline <- rcspline.eval(0:(latency-1), knots = knots, inclx = TRUE)
    B[, 2:dim(B)[2]] <- as.matrix(spline)

    # select B columns only selected into the final model
    B_lag <- B[which(B[,2] == lag), s]

    log_HR <- t(coef) %*% B_lag
  }



  return(log_HR)
}
