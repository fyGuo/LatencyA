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
#' @param adjusted_variable A vector of adjusted variables
#' @param adjusted_model A vector of adjusted models
#' @return A list. The first element is the best knots with the lowest AIC. The second element is a Cox model
#' @export
#' @examples
#' time_start <- "age_start"
#' time_end <- "age_end"
#' status <- "failure"
#' exposure <- paste0("lag", 0:15)
#' knots <- 0:15
#' latency <- 16
#' adjusted_variable <- "L"
#' adjusted_model <- "L"
#' fit <- CoxTermsearch(sim_data, time_start, time_end, status, exposure, knots = knots, latency)

#' #In this example, we will generate a random covariate L and adjust for it in our model
#' sim_data$L <- rbinom(nrow(sim_data), 1, 0.5)
#' adjusted_variable <- "L"
#' adjusted_model <- "L+I(L^2)"
#' fit <- CoxTermsearch(sim_data, time_start, time_end, status, exposure, knots = knots, latency,
#' adjusted_variable = adjusted_variable, adjusted_model = adjusted_model)
#' #L is selected into the final model but L^2 was not
CoxTermsearch <- function(data, time_start, time_end, status, exposure,
                          knots, latency, adjusted_variable = NULL, adjusted_model = NULL) {
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
  Sum_mat <- as.data.frame(Sum_mat)


  # specifiy the full model
  # then make a data frame to fit the regression
  regression_data <- cbind(data[,c(time_start, time_end, status,adjusted_variable)],  Sum_mat )

  model <- paste0("Surv(", time_start, ",", time_end, ",", status, ")")
  # include the cumulative exposure name
  model <- paste0(model, "~", paste0(colnames(Sum_mat), collapse = "+"))

  # include the adjusted variables if given
  if (length(adjusted_model) > 0) {model <- paste0(model, "+", adjusted_model)}


  full_model <- coxph(as.formula(model),
                      data = regression_data,
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


#' #In this example, we will generate a random covariate L and adjust for it in our model
#' sim_data$L <- rbinom(nrow(sim_data), 1, 0.5)
#' adjusted_variable <- "L"
#' adjusted_model <- "L+I(L^2)"
#' fit <- CoxTermsearch(sim_data, time_start, time_end, status, exposure, knots = knots, latency,
#' adjusted_variable = adjusted_variable, adjusted_model = adjusted_model)
#' #L is selected into the final model but L^2 was not
#' extract_CoxTermsearch(fit, 10, latency = latency, knot = knots)

extract_CoxTermsearch  <- function(fit, lag, latency, knots) {
  coef <- coef(fit)

  # first we want to extract the coefficient for the exposure
  # we extract the coefficients with "VX", where X stands for numbers
  coef <- coef[stringr::str_detect(names(coef), "\\bV\\d+\\b")]

  # replace NA with 0 which means the coefficient has no effects
  coef[is.na(coef)] <- 0

  # if the selected model contains no predictors, it means that no effect and logHR is 0
  if (is.null(coef)){
    log_HR = 0
  } else {

    # check terms that is used in the final model
    s <- stringr::str_extract(names(coef), "\\d+\\b") %>% as.numeric()

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

#' This function conduct bootstraps for @CoxTermsearch
#' @import survival
#' @import Hmisc
#' @import dplyr
#' @import stats
#' @import future
#' @import furrr
#' @importFrom MASS stepAIC
#' @param data A data frame containing the data
#' @param time_start start time
#' @param time_end end time
#' @param status survival status. 1 for death and 0 for censored
#' @param exposure Name of the time-varying exposure variable. From the most recent to the furthest in time.
#' @param knots A vector of prespecified knots
#' @param latency prespecified latency
#' @param adjusted_variable A vector of adjusted variables
#' @param adjusted_model A vector of adjusted models
#' @param lag A value of lag time
#' @param parallel A boolean numeric to indicate whether to run the function in parallel
#' @param boot_iter Number of bootstrap iterations
#' @param id A string of the column name that contains the subject ID
#' @return A data.frame. The first column is the mean of log HR. The second column is the variance of log HR
#' @export
#' @examples
#' time_start <- "age_start"
#' time_end <- "age_end"
#' status <- "failure"
#' exposure <- paste0("lag", 0:15)
#' knots <- 0:15
#' latency <- 16
#' adjusted_variable <- "L"
#' adjusted_model <- "L"
#' lag <- 2
#' boot_iter <- 10
#' id <- "id"
#' parallel <- TRUE
#' result <- CoxTermsearch_boot(sim_data, time_start, time_end, status, exposure, knots = knots, latency,
#' lag = lag, parallel = parallel, boot_iter = boot_iter, id = id)

CoxTermsearch_boot  <- function(data, time_start, time_end, status, exposure,
                           knots, latency, adjusted_variable = NULL, adjusted_model = NULL,
                           lag, parallel = FALSE, boot_iter = 10, id = "id") {
  ids <- unique(sim_data[,id])
  if (parallel == FALSE) {
    log_HR <- numeric(boot_iter)
    for (i in 1:boot_iter) {
      boot_ids <- sample(ids, size = length(ids), replace = TRUE)
      data<-sim_data[ sim_data[,id]  %in% boot_ids,]
      fit <- CoxTermsearch(data, time_start, time_end, status, exposure, knots = knots, latency)
      log_HR[i] <- extract_CoxTermsearch(fit, lag, latency, knots)
    }

    data.frame(log_HR = mean(log_HR), log_HR_var = var(log_HR)) %>% return()
  } else{
    plan("multicore")
    log_HR <- furrr::future_map_dbl(1:boot_iter, ~{
      boot_ids <- sample(ids, size = length(ids), replace = TRUE)
      data<-sim_data[ sim_data[,id]  %in% boot_ids,]
      fit <- CoxTermsearch(data, time_start, time_end, status, exposure, knots = knots, latency)
      log_HR <- extract_CoxTermsearch(fit, lag, latency, knots)
      return(log_HR)
    },
    .options = furrr_options(seed = T))

    data.frame(log_HR = mean(log_HR), log_HR_var = var(log_HR)) %>% return()
  }
}

