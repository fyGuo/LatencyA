pkgname <- "LatencyA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('LatencyA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CoxKnotsearch")
### * CoxKnotsearch

flush(stderr()); flush(stdout())

### Name: CoxKnotsearch
### Title: This function implments GridSearch of natural cubic splines in
###   Cox regression for latency analysis This is GreadSearching by knots
### Aliases: CoxKnotsearch

### ** Examples

time_start <- "age_start"
time_end <- "age_end"
status <- "failure"
exposure <- paste0("lag", 0:19)
# knots can be set as as a vector of potential number of knots, here we use 3 to save time
knots_number <- c(3)
latency <- 20
fit <- CoxKnotsearch(sim_data, time_start, time_end, status, exposure, knots_number, latency)



cleanEx()
nameEx("CoxKnotsearch_boot")
### * CoxKnotsearch_boot

flush(stderr()); flush(stdout())

### Name: CoxKnotsearch_boot
### Title: This function conduct bootstraps for @CoxKnotsearch
### Aliases: CoxKnotsearch_boot

### ** Examples

time_start <- "age_start"
time_end <- "age_end"
status <- "failure"
exposure <- paste0("lag", 0:15)
# knots can be set as as a vector of potential number of knots, here we use 3 to save time
knots_number <- c(3)
latency <- 16
adjusted_variable <- "L"
adjusted_model <- "L"
lag <- 2
boot_iter <- 1
id <- "id"
parallel <- TRUE
result <- CoxKnotsearch_boot(sim_data, time_start, time_end, status, 
exposure, knots_number = knots_number, latency,
lag = lag, parallel = parallel, boot_iter = boot_iter, id = id)



cleanEx()
nameEx("CoxNCSpline")
### * CoxNCSpline

flush(stderr()); flush(stdout())

### Name: CoxNCSpline
### Title: This function implments natural cubic splines in Cox regression
###   for latency analysis
### Aliases: CoxNCSpline

### ** Examples

time_start <- "age_start"
time_end <- "age_end"
status <- "failure"
exposure <- paste0("lag", 0:19)
knots_number <- 3
latency <- 20
fit <- CoxNCSpline(sim_data, time_start, time_end, status, exposure, knots_number, latency)



cleanEx()
nameEx("CoxPoly")
### * CoxPoly

flush(stderr()); flush(stdout())

### Name: CoxPoly
### Title: This is a function to implement polynomial Cox regression for
###   latency analysis
### Aliases: CoxPoly

### ** Examples

time_start <- "age_start"
time_end <- "age_end"
status <- "failure"
exposure <- paste0("lag", 0:19)
degree <- 3
latency <- 20
fit <- CoxPoly(sim_data, time_start, time_end, status, exposure, degree, latency)
fit



cleanEx()
nameEx("CoxTermsearch")
### * CoxTermsearch

flush(stderr()); flush(stdout())

### Name: CoxTermsearch
### Title: This function implments Stepwise model selection by terms of
###   natural cubic splines in Cox regression for latency analysis This is
###   GreadSearching by knots
### Aliases: CoxTermsearch

### ** Examples

time_start <- "age_start"
time_end <- "age_end"
status <- "failure"
exposure <- paste0("lag", 0:15)
knots <- seq(0, 15, by = 3)
latency <- 16
adjusted_variable <- "L"
adjusted_model <- "L"
fit <- CoxTermsearch(sim_data, time_start, time_end, status, exposure, knots = knots, latency)
#In this example, we will generate a random covariate L and adjust for it in our model
sim_data$L <- rbinom(nrow(sim_data), 1, 0.5)
adjusted_variable <- "L"
adjusted_model <- "L+I(L^2)"
fit <- CoxTermsearch(sim_data, time_start, time_end, status, exposure, knots = knots, latency,
adjusted_variable = adjusted_variable, adjusted_model = adjusted_model)
#L is selected into the final model but L^2 was not



cleanEx()
nameEx("CoxTermsearch_boot")
### * CoxTermsearch_boot

flush(stderr()); flush(stdout())

### Name: CoxTermsearch_boot
### Title: This function conduct bootstraps for @CoxTermsearch
### Aliases: CoxTermsearch_boot

### ** Examples

time_start <- "age_start"
time_end <- "age_end"
status <- "failure"
exposure <- paste0("lag", 0:15)
knots <- seq(0, 15, by = 3)
latency <- 16
adjusted_variable <- "L"
adjusted_model <- "L"
lag <- 2
boot_iter <- 1
id <- "id"
parallel <- TRUE
result <- CoxTermsearch_boot(sim_data, time_start, time_end, 
status, exposure, knots = knots, latency,
lag = lag, parallel = parallel, boot_iter = boot_iter, id = id)



cleanEx()
nameEx("extract_CoxKnotsearch")
### * extract_CoxKnotsearch

flush(stderr()); flush(stdout())

### Name: extract_CoxKnotsearch
### Title: This is a function to extract log HR based on the resuls from
###   CoxKnotsearch
### Aliases: extract_CoxKnotsearch

### ** Examples

time_start <- "age_start"
time_end <- "age_end"
status <- "failure"
exposure <- paste0("lag", 0:19)
# knots can be set as as a vector of potential number of knots, here we use 3 to save time
knots_number <- c(3)
latency <- 20
fit <- CoxKnotsearch(sim_data, time_start, time_end, status, exposure, knots_number, latency)
extract_CoxKnotsearch(fit, 10, latency = latency)
extract_CoxKnotsearch(fit, c(0, 5, 10), latency = latency)



cleanEx()
nameEx("extract_CoxNCSpline")
### * extract_CoxNCSpline

flush(stderr()); flush(stdout())

### Name: extract_CoxNCSpline
### Title: This is a function to extract log HR based on the resuls from
###   CoxNCSpline
### Aliases: extract_CoxNCSpline

### ** Examples

time_start <- "age_start"
time_end <- "age_end"
status <- "failure"
exposure <- paste0("lag", 0:19)
knots_number <- 3
latency <- 20
fit <- CoxNCSpline(sim_data, time_start, time_end, status, exposure, knots_number, latency)
extract_CoxNCSpline(fit, 10, latency = latency)
extract_CoxNCSpline(fit, c(0, 5, 10), latency = latency)



cleanEx()
nameEx("extract_CoxPoly")
### * extract_CoxPoly

flush(stderr()); flush(stdout())

### Name: extract_CoxPoly
### Title: This is a function to extract log HR based on the resuls from
###   CoxPoly
### Aliases: extract_CoxPoly

### ** Examples

time_start <- "age_start"
time_end <- "age_end"
status <- "failure"
exposure <- paste0("lag", 0:19)
degree <- 3
latency <- 20
fit <- CoxPoly(sim_data, time_start, time_end, status, exposure, degree, latency)
extract_CoxPoly(fit, 10)
extract_CoxPoly(fit, c(0, 5, 10))



cleanEx()
nameEx("extract_CoxTermsearch")
### * extract_CoxTermsearch

flush(stderr()); flush(stdout())

### Name: extract_CoxTermsearch
### Title: This is a function to extract log HR based on the resuls from
###   CoxTermsearch
### Aliases: extract_CoxTermsearch

### ** Examples

time_start <- "age_start"
time_end <- "age_end"
status <- "failure"
exposure <- paste0("lag", 0:19)
knots <- seq(0, 15, by = 3)
latency <- 20
fit <- CoxTermsearch(sim_data, time_start, time_end, status, exposure, knots = knots, latency)
extract_CoxTermsearch(fit, 10, latency = latency, knot = knots)
extract_CoxTermsearch(fit, c(0, 5, 10), latency = latency, knot = knots)
#In this example, we will generate a random covariate L and adjust for it in our model
sim_data$L <- rbinom(nrow(sim_data), 1, 0.5)
adjusted_variable <- "L"
adjusted_model <- "L+I(L^2)"
fit <- CoxTermsearch(sim_data, time_start, time_end, status, exposure, knots = knots, latency,
adjusted_variable = adjusted_variable, adjusted_model = adjusted_model)
#L is selected into the final model but L^2 was not
extract_CoxTermsearch(fit, 10, latency = latency, knot = knots)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
