################################################################################
# 20180410 - Functions to parse files from TopPic 
################################################################################

require(data.table)


# MS1 files:

# fname <- "//GANDALF/mloca/RRelatedWork/VisioProt-MS/files/test/TopPic2/OFJMX160905_30.raw_ms1.msalign"

TopPicMS1Parsing <- function(fname) {
  # Return a table in the style of RoWinPro tables for use in VisioProt.
  # fname is the path to the file to parse.
  allData <- readLines(fname)
  allData <- allData[-(1:11)]
  rep_ions_entries = which(allData=="BEGIN IONS")
  IDs <- gsub("ID=", "", allData[rep_ions_entries+1])
  SCANs <- gsub("SCANS=", "", allData[rep_ions_entries+2])
  RT <- gsub("RETENTION_TIME=", "", allData[rep_ions_entries+3])
  removeEntries <- c(rep_ions_entries,rep_ions_entries+1,rep_ions_entries+2,rep_ions_entries+3,rep_ions_entries[2:length(rep_ions_entries)]-1,rep_ions_entries[2:length(rep_ions_entries)]-2, length(allData)-1, length(allData))
  ions_per_scan <- diff(rep_ions_entries) - 6
  ions_per_scan <- c(ions_per_scan, (length(allData) - rep_ions_entries[length(rep_ions_entries)] - 5))
  dat <- fread(paste(allData[-removeEntries], collapse = "\n"), sep = "\t")
  class(dat) <- "data.frame"
  names(dat) <- c("Mass", "intensity", "charge")
  dat$ID <- rep(IDs, ions_per_scan)
  dat$SCANs <- rep(SCANs, ions_per_scan)
  dat$RT <- rep(RT, ions_per_scan)
  dat <- dat[,c(6,1,2,3,5)]
  # Keep only the ions >= 5+
  dat <- dat[dat[,4]>=5,]
  # Change from seconds to minutes:
  dat[,1] <- as.numeric(dat[,1])/60
  # Keep only the 30% highest intensities:
  dat <- dat[order(dat[,3], decreasing = T),]
  dat <- dat[!is.na(dat[,3]),]
  thresh <- floor(0.3 * nrow(dat))
  dat <- dat[c(1:thresh),]
  return(dat)
}

# MS2 files:

# fname <- "//GANDALF/mloca/RRelatedWork/VisioProt-MS/files/test/TopPic2/OFJMX160905_30.raw_ms2.msalign"

TopPicMS1Parsing <- function(fname) {
  # Return a table in the style of RoWinPro tables for use in VisioProt.
  # fname is the path to the file to parse.
  allData <- readLines(fname)
  allData <- allData[-(1:11)]
  rep_ions_entries = which(allData=="BEGIN IONS")
  IDs <- gsub("ID=", "", allData[rep_ions_entries+1])
  SCANs <- gsub("SCANS=", "", allData[rep_ions_entries+2])
  RT <- gsub("RETENTION_TIME=", "", allData[rep_ions_entries+3])
  Mass <- gsub("PRECURSOR_MASS=", "", allData[rep_ions_entries+9])
  intensity <- gsub("PRECURSOR_INTENSITY=", "", allData[rep_ions_entries+10])
  charge <- gsub("PRECURSOR_CHARGE=", "", allData[rep_ions_entries+8])
  
  dat <- data.frame("RT"=RT, "Mass"=Mass, "intensity"=intensity, "Scan"=SCANs, stringsAsFactors = F)
  # Change from seconds to minutes:
  dat[,1] <- as.numeric(dat[,1])/60
  # Keep only the 30% highest intensities:
  #dat <- dat[order(dat[,3], decreasing = T),]
  #dat <- dat[!is.na(dat[,3]),]
  #thresh <- floor(0.3 * nrow(dat))
  #dat <- dat[c(1:thresh),]
  return(dat)
}
