barplot2 <- function(xs, ys, ..., space = 0, horiz = FALSE) {
  mids <- diff(xs) / 2
  mids <- c(mids[1], mids) - space/2
  if (horiz) {
    rect(0, xs - mids, ys, xs + mids, ...)
  } else {
    rect(xs - mids, 0, xs + mids, ys, ...)
  }
}
