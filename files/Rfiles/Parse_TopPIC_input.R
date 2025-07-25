TopPicMS1Parsing <- function(fname) {
  cat("== Start parsing TopPIC data MS1 ==\n")
  # Return a table in the style of RoWinPro tables for use in VisioProt.
  # fname is the path to the file to parse.
  allData <- readLines(fname)
  numline <- which(grepl("^[^P]+Parameters ", allData, perl = T))[2]
  allData <- allData[-(1:numline)]
  rep_ions_entries = which(allData=="BEGIN IONS")
  rep_ions_end = which(allData=="END IONS")
  IDs <- gsub("ID=", "", allData[grepl("^ID", allData)])
  SCANs <- gsub("SCANS=", "", allData[grepl("^SCANS", allData)])
  RT <- gsub("RETENTION_TIME=", "", allData[grepl("^RETENTION_TIME", allData)])
  # removeEntries <- c(rep_ions_entries,rep_ions_entries+1,rep_ions_entries+2,rep_ions_entries+3,rep_ions_entries[2:length(rep_ions_entries)]-1,rep_ions_entries[2:length(rep_ions_entries)]-2, length(allData)-1, length(allData))
  removeEntries <- which(!(substr(allData, 1, 1) %in% c(0:9))) # Remove all lines not starting with a number
  # Count the number of lines with comment per spectrum:
  numComments <- sum(!(substr(allData[rep_ions_entries[1]:rep_ions_end[1]], 1, 1) %in% c(0:9)))
  ions_per_scan <- sapply(seq_along(rep_ions_entries), function(x) {
    rep_ions_end[x] - rep_ions_entries[x] - numComments + 1
  })
  dat <- fread(paste(allData[-removeEntries], collapse = "\n"), sep = "\t")
  class(dat) <- "data.frame"
  names(dat) <- c("Mass", "intensity", "charge")
  dat$ID <- rep(IDs, ions_per_scan)
  dat$SCANs <- rep(SCANs, ions_per_scan)
  dat$RT <- rep(RT, ions_per_scan)
  dat <- dat[,c(6,1,2,3,5)] 
  c("Keep only the ions >= 5+\n")
  dat <- dat[dat[,4]>=5,]
  c("Change from seconds to minutes\n")
  dat[,1] <- as.numeric(dat[,1])/60
  c("Keep only the 100% highest intensities\n")
  dat <- dat[order(dat[,3], decreasing = T),]
  dat <- dat[!is.na(dat[,3]),]
  thresh <- floor(1 * nrow(dat))
  dat <- dat[c(1:thresh),]
  # For the functions to come (thresholding, renaming):
  dat[,4] <- rep(NA, nrow(dat))
  dat[,5] <- rep(NA, nrow(dat))
  cat("== End parsing TopPIC MS1 ==\n\n")
  return(dat)
  
}

TopPicMS2Parsing <- function(fname) {
  cat("== Start parsing TopPIC data MS2 ==\n")
  # Return a table in the style of RoWinPro tables for use in VisioProt.
  # fname is the path to the file to parse.
  allData <- readLines(fname)
  numline <- which(grepl("^[^P]+Parameters ", allData, perl = T))[2]
  allData <- allData[-(1:numline)]
  rep_ions_entries = which(allData=="BEGIN IONS")
  rep_ions_end = which(allData=="END IONS")
  IDs <- gsub("ID=", "", allData[grepl("^ID", allData)])
  SCANs <- gsub("SCANS=", "", allData[grepl("^SCANS", allData)])
  RT <- gsub("RETENTION_TIME=", "", allData[grepl("^RETENTION_TIME", allData)])
  Mass <- gsub("PRECURSOR_MASS=", "", allData[grepl("^PRECURSOR_MASS", allData)])
  intensity <- gsub("PRECURSOR_INTENSITY=", "", allData[grepl("^PRECURSOR_INTENSITY", allData)])
  charge <- gsub("PRECURSOR_CHARGE=", "", allData[grepl("^PRECURSOR_CHARGE", allData)])
  
  dat <- data.frame("RT"=RT, "Mass"=Mass, "intensity"=intensity, "Scan"=SCANs, stringsAsFactors = F)
  c("Change from seconds to minutes\n")
  dat[,1] <- as.numeric(dat[,1])/60
  cat("== End parsing TopPIC MS2 ==\n\n")
  return(dat)
  
}

TopPicIDParsing <- function(fname) {
  
  cat("== Start parsing TopPIC data ID ==\n")
  # Return a table in the style of RoWinPro tables for use in VisioProt.
  # fname is the path to the file to parse.
  allData <- readLines(fname)
  numline <- which(grepl("^[^P]+Parameters ", allData, perl = T))[2]
  allData <- allData[-(1:numline)]
  allData[1] <- gsub("#", "", allData[1])
  dat <- fread(paste(allData, collapse = "\n"), sep = "\t", header = T, stringsAsFactors = F)
  class(dat) <- "data.frame"
  cat("== End parsing TopPIC ID ==\n\n")
  return(dat)
  
}
