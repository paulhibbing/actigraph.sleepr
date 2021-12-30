#' Impute missing count values
#'
#' Trim leading and trailing NAs. Fill in the rest of the NAs using
#' cubic spline interpolation.
#' @param agdb A \code{tibble} (\code{tbl}) of activity data (at least)
#' an \code{epochlength} attribute.
#' @param ... Comma separated list of unquoted variables.
#' @return A \code{tibble} (\code{tbl}) of activity data. Each variable
#' in \code{...} is imputed.
#' @seealso \code{\link[zoo]{na.spline}}, \code{\link[zoo]{na.trim}}
#' @examples
#' library("dplyr")
#' data("gtxplus1day")
#'
#' gtxplus1day$axis1[5:10] <- NA
#' gtxplus1day %>%
#'   impute_epochs(axis1)
#' @export
impute_epochs <- function(agdb, ...) {
  selected <- tidyselect::vars_select(names(agdb), ...)
  if (length(selected) == 0) {
    return(agdb)
  }

  agdb %>%
    group_modify(
      ~ impute_epochs_(., selected)
    )
}
impute_epochs_ <- function(data, selected) {
  impute <- function(x) pmax(round(na.spline(x)), 0)

  data %>%
    select(timestamp, !!!selected) %>%
    na.trim() %>%
    inner_join(data,
      by = c("timestamp", selected)
    ) %>%
    mutate(
      across(all_of(selected), impute)
    )
}
#' Checks whether there are gaps in the time series
#'
#' The timestamps in the agd time series should run from
#' \code{first(timestamp)} to \code{last(timestamp)} in increments of
#' \code{epochlength} seconds. This function checks whether this holds
#' or not. If the data is grouped (e.g., by subject), the check is
#' performed for each group separately.
#' @param agdb A \code{tibble} (\code{tbl}) of activity data (at least)
#' an \code{epochlength} attribute.
#' @inheritParams get_epoch_length
#' @return True or false.
#' @export
has_missing_epochs <- function(agdb, TSname = "timestamp") {
  if (anyNA(agdb[[TSname]])) {
    return(TRUE)
  }

  agdb <- agdb %>%
    group_modify(
      function(x,y,TSname) has_missing_epochs_(x, TSname),
      TSname = TSname
    )
  any(agdb$missing)
}

has_missing_epochs_ <- function(data, TSname = "timestamp") {
  epoch_len <- get_epoch_length(data, TSname)

  epochs <- seq(first(data[[TSname]]), last(data[[TSname]]),
    by = epoch_len
  )
  tibble::tibble(missing = !identical(epochs, data[[TSname]]))
}
