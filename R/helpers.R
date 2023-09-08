
get_response_type <- function(y) {
  ret <- if (is.ordered(y))
    "ordered"
  else if (is.integer(y))
    "count"
  else if (inherits(y, "Surv"))
    "survival"
  else
    "continuous"
  ret
}

get_order <- function(response_type, y) {
  ret <- if (response_type == "ordered")
    nlevels(y) - 1L
  else
    10
  ret
}
