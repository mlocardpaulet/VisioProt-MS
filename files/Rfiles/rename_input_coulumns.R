RenameRoWinPro <- function(tab) {
  # Filter and rename the tables RoWinPro style
  names(tab)[3] <- "intensity"
  names(tab)[2] <- "Mass"
  names(tab)[1] <- "RT"
  tab
}

RenameBioPharma <- function(tab) {
  # Filter and rename the tables BioPharma style
  vec <- names(tab)
  vec[3] <- "intensity"
  vec[2] <- "Mass"
  vec[1] <- "RT"
  vec[4] <- "PeakStart"
  vec[5] <- "PeakStop" 
  names(tab) <- vec
  return(tab)
}
