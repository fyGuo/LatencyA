#' This function implments Stepwise model selection by terms of natural cubic splines in Cox regression for latency analysis
#' This is GreadSearching by knots
#' @importFrom survival coxph coxph.control
#' @importFrom Hmisc rcspline.eval
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
#' @param criteria A string to specify the criteria for model selection. It can be either "AIC" or "BIC". The default is "AIC"
#' @return A list. The first element is the best knots with the lowest AIC. The second element is a Cox model
#' @export
#' @examples
#' time_start <- "age_start"
#' time_end <- "age_end"
#' status <- "failure"
#' exposure <- paste0("lag", 0:15)
#' knots <- seq(0, 15, by = 3)
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
                          knots, latency, adjusted_variable = NULL, adjusted_model = NULL,
                          criteria = "AIC") {
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
  if (criteria == "AIC") {
    fit_best <- stepAIC(full_model, direction = "both", trace = FALSE)
  } else if (criteria == "BIC") {
    fit_best <- stepAIC(full_model, direction = "both", trace = FALSE, k = log(nrow(regression_data)))
  }

  return(fit_best)
}

#' This is a function to extract log HR based on the resuls from CoxTermsearch
#' @param fit A model from CoxKnotsearch
#' @param lag A lag time, or a vector of lag times
#' @param latency prespecified latency
#' @param knots A vector of prespecified knots
#' @return A data frame with one row per lag, giving the estimated log HR at each lag time
#' @import stringr
#' @export
#' @examples
#' time_start <- "age_start"
#' time_end <- "age_end"
#' status <- "failure"
#' exposure <- paste0("lag", 0:19)
#' knots <- seq(0, 15, by = 3)
#' latency <- 20
#' fit <- CoxTermsearch(sim_data, time_start, time_end, status, exposure, knots = knots, latency)
#' extract_CoxTermsearch(fit, 10, latency = latency, knot = knots)
#' extract_CoxTermsearch(fit, c(0, 5, 10), latency = latency, knot = knots)


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
    log_HR = rep(0, length(lag))
  } else {
    # check terms that is used in the final model
    s <- stringr::str_extract(names(coef), "\\d+\\b") |> as.numeric()

    # make the B matrix
    v <- 0:(latency - 1)
    # generate the B matrix
    B <- matrix(NA, nrow = latency, ncol = length(knots)  )
    # the first column of B matrix is the intercept, and thus, is 1
    B[,1] <- 1

    spline <- rcspline.eval(0:(latency-1), knots = knots, inclx = TRUE)
    B[, 2:dim(B)[2]] <- as.matrix(spline)

    # one row of the basis per requested lag, keeping only the columns selected
    # into the final model (drop = FALSE keeps it a matrix even for a single lag)
    B_lag <- B[match(lag, v), s, drop = FALSE]

    log_HR <- as.vector(B_lag %*% coef)
  }



  return(data.frame(lag = lag, log_HR = log_HR))
}

#' This function conduct bootstraps for @CoxTermsearch
#' @importFrom survival coxph coxph.control
#' @importFrom Hmisc rcspline.eval
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
#' @param parallel_plan A string to specify the future plan used when parallel = TRUE. It can be either "multicore" or "multisession". The default is "multicore"
#' @param boot_iter Number of bootstrap iterations
#' @param id A string of the column name that contains the subject ID
#' @param criteria A string to specify the criteria for model selection. It can be either "AIC" or "BIC". The default is "AIC"
#' @return A data.frame. The first column is the mean of log HR. The second column is the variance of log HR, the third and
#' fourth columns are the 2.5 and 97.5 percentiles of log HR if parallel = TRUE
#' @export
#' @examples
#' time_start <- "age_start"
#' time_end <- "age_end"
#' status <- "failure"
#' exposure <- paste0("lag", 0:15)
#' knots <- seq(0, 15, by = 3)
#' latency <- 16
#' adjusted_variable <- "L"
#' adjusted_model <- "L"
#' lag <- 2
#' boot_iter <- 1
#' id <- "id"
#' parallel <- TRUE
#' result <- CoxTermsearch_boot(sim_data, time_start, time_end, 
#' status, exposure, knots = knots, latency,
#' lag = lag, parallel = parallel, boot_iter = boot_iter, id = id)

CoxTermsearch_boot  <- function(data, time_start, time_end, status, exposure,
                           knots, latency, adjusted_variable = NULL, adjusted_model = NULL,
                           lag, parallel = FALSE, parallel_plan = c("multicore", "multisession"),
                           boot_iter = 10, id = "id", criteria = "AIC") {
  parallel_plan <- match.arg(parallel_plan)
  if (parallel == FALSE) {
    log_HR <- numeric(boot_iter)
    for (i in 1:boot_iter) {
      temp <- resample_clusters(data, id = id)
      fit <- CoxTermsearch(temp, time_start, time_end, status, exposure, knots = knots, latency,
                           adjusted_variable = adjusted_variable, adjusted_model = adjusted_model, criteria = criteria)
      log_HR[i] <- extract_CoxTermsearch(fit, lag, latency, knots)$log_HR
    }

    return(data.frame(median_log_HR = median(log_HR), 
                      log_HR_var = var(log_HR),
                     percentile_2.5 = quantile(log_HR, 0.025),
                     percentile_97.5 = quantile(log_HR, 0.975)))
  } else{
    plan(parallel_plan)
    log_HR <- furrr::future_map_dbl(1:boot_iter, ~{
      temp <- resample_clusters(data, id = id)
      fit <- CoxTermsearch(temp, time_start, time_end, status, exposure, knots = knots, latency,
                           adjusted_variable = adjusted_variable, adjusted_model = adjusted_model, criteria = criteria)
      log_HR <- extract_CoxTermsearch(fit, lag, latency, knots)$log_HR
      return(log_HR)
    },
    .options = furrr_options(seed = T))

    return(data.frame(median_log_HR = median(log_HR), 
                      log_HR_var = var(log_HR),
                     percentile_2.5 = quantile(log_HR, 0.025),
                     percentile_97.5 = quantile(log_HR, 0.975)))
  }
}

