#' This is a function to help resemple clusters for bootstrap. It is adapted from the function "clusterboot" in the R package "fpc".
#' @param data A data frame containing the data
#' @param id A string of the column name that contains the subject ID
#' @return A data frame with the same structure as the input data, but with resampled clusters

resample_clusters <- function(data, id = "id") {
  ids   <- unique(data[[id]])
  draw  <- sample(ids, length(ids), replace = TRUE)          # with replacement
  # rows for each id, looked up once
  by_id <- split(seq_len(nrow(data)), data[[id]])
  parts <- vector("list", length(draw))
  for (j in seq_along(draw)) {
    rows         <- data[by_id[[as.character(draw[j])]], , drop = FALSE]
    rows[[id]]   <- j                                        # fresh unique id
    parts[[j]]   <- rows
  }
  do.call(rbind, parts)
}