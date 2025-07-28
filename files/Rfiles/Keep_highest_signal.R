ThresholdCleaning <- function(l, threshold) {
  # removes the points with lowest intensity according to the threshold chosen by the user:
  threshold <- threshold/100
  l1 <- list()
  for (x in seq_along(l)) {
    # Filter the intensities according to threshold:
    temp <- cbind(l[[x]], "File" = rep(names(l)[x], nrow(l[[x]])))
    temp <- temp[order(temp[,3], decreasing = T),]
    temp <- temp[!is.na(temp[,3]),]
    thresh <- floor(threshold * nrow(temp))
    temp <- temp[c(1:thresh),]
    l1[[x]] <- temp
  }
  names(l1) <- names(l)
  return(l1)
}
