RenameRoWinPro <- function(tab) {
  # Filter and rename the tables RoWinPro style
  names(tab)[3] <- "intensity"
  names(tab)[2] <- "Mass"
  names(tab)[1] <- "RT"
  return(tab)
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

RenamePDPrSM <- function(tab) {
  # Filter and rename the tables from PD
  # Change field names for compatibility between PSMs and PrSMs tables:
  if (sum(grepl("Master.Protein.Description", names(tab))) == 0) {
    names(tab)[names(tab) == "Protein.Accessions"] <- "Master.Protein.Descriptions"
  }
  if (sum(grepl("Precursor.MHplus.in.Da", names(tab))) == 0) {
    names(tab)[names(tab) == "Precursor.m.z..Da."] <- "Precursor.MHplus.in.Da"
    names(tab)[names(tab) == "mz.in.Da"] <- "Precursor.MHplus.in.Da"
  }
  names(tab)[names(tab) == "RT..min."] <- "RT.in.min"
  if (sum(names(tab) == "Modifications") > 0) {
    tab$Master.Protein.Descriptions <- paste0(tab$Master.Protein.Descriptions, "|", tab$Modifications)
  }
  return(tab)
}

RenamePDMSMS <- function(tab) {
  # Filter and rename the tables from PD
  names(tab)[names(tab) == "RT..min."] <- "RT.in.min"
  names(tab)[names(tab) == "Precursor.MH...Da."] <- "Precursor.MHplus.in.Da"
  # names(tab)[names(tab) == "Majority.Protein.Accessions"] <- "Protein.Accessions"
  return(tab)
}