#' This function implments GridSearch of natural cubic splines in Cox regression for latency analysis
#' This is GreadSearching by knots
#' @importFrom survival coxph coxph.control
#' @importFrom Hmisc rcspline.eval
#' @import stats
#' @param data A data frame containing the data
#' @param time_start start time
#' @param time_end end time
#' @param status survival status. 1 for death and 0 for censored
#' @param exposure Name of the time-varying exposure variable. From the most recent to the furthest in time.
#' @param  knots_number A vector of potential number of knots
#' @param latency pre-specified latency
#' @param adjusted_variable A vector of adjusted variables
#' @param adjusted_model A vector of adjusted models
#' @param criteria A string to specify the criteria for model selection. It can be either "AIC" or "BIC". The default is "AIC"
#' @return A list. The first element is the best knots with the lowest AIC. The second element is a Cox model
#' @export
#' @examples
#' time_start <- "age_start"
#' time_end <- "age_end"
#' status <- "failure"
#' exposure <- paste0("lag", 0:19)
#' # knots can be set as as a vector of potential number of knots, here we use 3 to save time
#' knots_number <- c(3)
#' latency <- 20
#' fit <- CoxKnotsearch(sim_data, time_start, time_end, status, exposure, knots_number, latency)
CoxKnotsearch <- function(data, time_start, time_end, status, exposure,
                          knots_number, latency, adjusted_variable = NULL, adjusted_model = NULL,
                          criteria = "AIC") {
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
      cut_start <- quantile(v, 0.1) |> round()
      cut_end <- quantile(v, 0.9) |> round()
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
      if (criteria == "AIC") {
        AIC <- extractAIC(fit)[2]
      } else if (criteria == "BIC") {
        AIC <- extractAIC(fit, k = log(nrow(regression_data)))[2]
      }
      if (AIC < mini_AIC) {
        mini_AIC <- AIC
        fit_best <- fit
        knots_best <- knots
      }
    }
  }
  return(list(knots_best, fit_best)) 
}


#' This is a function to extract log HR based on the resuls from CoxKnotsearch
#' @param fit A model from CoxKnotsearch
#' @param lag A lag time, or a vector of lag times
#' @param latency prespecified latency
#' @return A data frame with one row per lag, giving the estimated log HR at each lag time
#' @export
#' @examples
#' time_start <- "age_start"
#' time_end <- "age_end"
#' status <- "failure"
#' exposure <- paste0("lag", 0:19)
#' # knots can be set as as a vector of potential number of knots, here we use 3 to save time
#' knots_number <- c(3)
#' latency <- 20
#' fit <- CoxKnotsearch(sim_data, time_start, time_end, status, exposure, knots_number, latency)
#' extract_CoxKnotsearch(fit, 10, latency = latency)
#' extract_CoxKnotsearch(fit, c(0, 5, 10), latency = latency)
extract_CoxKnotsearch  <- function(fit, lag, latency) {

  #extract the knots
  knots <- fit[[1]]

  # extract the coefficients
  coef <- coef(fit[[2]])

  # first we want to extract the coefficient for the exposure
  # we extract the coefficients with "VX", where X stands for numbers
  coef <- coef[stringr::str_detect(names(coef), "\\bV\\d+\\b")]

  # make a function
  B <- matrix(NA, nrow = latency, ncol = length(knots) )
  # the first column of B matrix is the intercept, and thus, is 1
  B[,1] <- 1
  spline <- rcspline.eval(0:(latency-1), inclx = TRUE, knots = knots)
  B[, 2:length(knots)] <- as.matrix(spline)

  # one row of the basis per requested lag (drop = FALSE keeps it a matrix even
  # for a single lag, so the scalar case is just the 1-row case)
  B_lag <- B[match(lag, 0:(latency - 1)), , drop = FALSE]

  log_HR <- as.vector(B_lag %*% coef)

  return(data.frame(lag = lag, log_HR = log_HR))
}

#' This function conduct bootstraps for @CoxKnotsearch
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
#' @param knots_number A vector of prespecified knot number
#' @param latency prespecified latency
#' @param adjusted_variable A vector of adjusted variables
#' @param adjusted_model A vector of adjusted models
#' @param lag A value of lag time
#' @param parallel A boolean numeric to indicate whether to run the function in parallel
#' @param parallel_plan A string to specify the future plan used when parallel = TRUE. It can be either "multicore" or "multisession". The default is "multicore"
#' @param boot_iter Number of bootstrap iterations
#' @param id A string of the column name that contains the subject ID
#' @param criteria A string to specify the criteria for model selection. It can be either "AIC" or "BIC". The default is "AIC"
#' @return A data.frame. The first column is the mean of log HR. The second column is the variance of log HR, the third and fourth columns are the 2.5 and 97.5 percentiles of log HR if parallel = TRUE
#' @export
#' @examples
#' time_start <- "age_start"
#' time_end <- "age_end"
#' status <- "failure"
#' exposure <- paste0("lag", 0:15)
#' # knots can be set as as a vector of potential number of knots, here we use 3 to save time
#' knots_number <- c(3)
#' latency <- 16
#' adjusted_variable <- "L"
#' adjusted_model <- "L"
#' lag <- 2
#' boot_iter <- 1
#' id <- "id"
#' parallel <- TRUE
#' result <- CoxKnotsearch_boot(sim_data, time_start, time_end, status, 
#' exposure, knots_number = knots_number, latency,
#' lag = lag, parallel = parallel, boot_iter = boot_iter, id = id)


CoxKnotsearch_boot  <- function(data, time_start, time_end, status, exposure,
                                knots_number, latency, adjusted_variable = NULL, adjusted_model = NULL,
                                lag, parallel = FALSE, parallel_plan = c("multicore", "multisession"),
                                boot_iter = 10, id = "id", criteria = "AIC") {
  parallel_plan <- match.arg(parallel_plan)
  ids <- unique(data[,id])
  if (parallel == FALSE) {
    log_HR <- numeric(boot_iter)
    for (i in 1:boot_iter) {
      temp<-resample_clusters(data, id = id)
      fit <- CoxKnotsearch(temp, time_start, time_end, status, exposure,  knots_number =  knots_number, latency, criteria = criteria)
      log_HR[i] <- extract_CoxKnotsearch(fit, lag, latency)$log_HR
    }

    return(return(data.frame(median_log_HR = median(log_HR), 
                      log_HR_var = var(log_HR),
                     percentile_2.5 = quantile(log_HR, 0.025),
                     percentile_97.5 = quantile(log_HR, 0.975))))
  } else{
    plan(parallel_plan)
    log_HR <- furrr::future_map_dbl(1:boot_iter, ~{
      temp <- resample_clusters(data, id = id)
      fit <- CoxKnotsearch(temp, time_start, time_end, status, exposure, knots_number = knots_number, latency, criteria = criteria)
      log_HR <- extract_CoxKnotsearch(fit, lag, latency)$log_HR
      return(log_HR)
    },
    .options = furrr_options(seed = T))

    return(data.frame(median_log_HR = median(log_HR), 
                      log_HR_var = var(log_HR),
                     percentile_2.5 = quantile(log_HR, 0.025),
                     percentile_97.5 = quantile(log_HR, 0.975)))  
  }
}


