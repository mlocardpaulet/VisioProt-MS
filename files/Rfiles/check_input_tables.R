

#---- MS trace input tables (deconvoluted) ----#

# Find what software the file comes from: --------------------------------------

is_biopharama_MS_input <- function(input = NULL) {
  #' @description
    #' this function returns a T/F indicating if all the column headers necessary
    #' for VisioProt-MS to plot a MS trace from a biopharmaFinder file are present
    #' To do so, it reads the first line of the file
  #' @param input path to a file
  
  grepl("Mass", readLines(input)[1]) & 
    grepl("Apex RT", readLines(input)[1]) & 
    grepl("Sum Intensity", readLines(input)[1]) & 
    grepl("Start Time (min)", readLines(input)[1], fixed = T) & 
    grepl("Stop Time (min)", readLines(input)[1], fixed = T)
  
}

is_Bruker_MS_input <- function(input = NULL) {
  #' @description
  #' this function returns a T/F indicating if all the column headers necessary
  #' for VisioProt-MS to plot a MS trace from a Bruker file are present
  #' To do so, it reads the second line of the file
  #' @param input path to a file
  
  substr(readLines(input)[2], 0, 13) == "Compound Name"
  
}

is_ProMex_MS_input <- function(input = NULL) {
  #' @description
  #' this function returns a T/F indicating if all the column headers necessary
  #' for VisioProt-MS to plot a MS trace from a ProMex file are present
  #' To do so, it reads the first line of the file
  #' @param input name of a file
  
  grepl(".ms1ft", input, fixed = T)
  
}

is_TopFD_input <- function(input = NULL) {
  #' @description
  #' this function returns a T/F indicating if all the column headers necessary
  #' for VisioProt-MS to plot a MS trace from a TopFD file are present
  #' To do so, it reads the first line of the file
  #' @param input name of a file
  
  grepl("ms1.msalign", input, fixed = T)
  
}

# Validataion and error messages: ----------------------------------------------

RoWinPro_MS_check <- function(input = NULL) { 
  #' @description
    #' this function checks that the input table has exactly 4 columns
  #' @param input path to the rowinpro file
  
  line_1 <- readLines(input)[1]
 
  validate(need(length(as.numeric(gregexpr("\t", line_1)[[1]]))==3, 
                paste0("Error in file format for plotting MS trace from Rowinpro:\n",
                       "You don't have the expected number of columns (4).")
  ))
 
}

Bruker_MS_check <- function(input = NULL) { 
  #' @description
    #' this function checks that the input table has the column "RT" and has no ","
    #' in the first line.
  #' @param input path to the rowinpro file
  
  validate(need(!grepl(",", readLines(input)[1]) & grepl(" RT", readLines(input)[2]),
           paste0("Error in file format for plotting MS trace from Bruker:\n",
                  "Check that there is no \",\" in the first line of the input file.")
  ))
  
}

ProMex_MS_check <- function(input = NULL) { 
  #' @description
  #' this function checks that the input table has the column "MinElutionTime" 
  #' in the first line.
  #' @param input path to the rowinpro file
  
  validate(need(grepl("MinElutionTime", readLines(input)[1]),
                paste0("Error in file format for plotting MS trace from ProMex:\n",
                       "Check that the column \"MinElutionTime\" is present in the input file.")
  ))
  
}

TopFD_MS_check <- function(input = NULL) { 
  #' @description
  #' this function checks that the input table has the column "MinElutionTime" 
  #' in the first line.
  #' @param input path to the rowinpro file
  
  validate(need(grepl("#TopFD", readLines(input)[1], fixed = T),
                paste0("Error in file format for plotting MS trace from TopFD:\n",
                       "Check that the column \"#TopFD\" is present in the input file.")
  ))
  
}

#---- output tables from top down analysis ----#

# Proteome Discoverer: ---------------------------------------------------------

PD_MS2_check <- function(PSM_tab = NULL, MSMS_tab = NULL) { # works for v3.0
  #' @description
    #' this function checks that the PrSM and MSMS tables outputed from PD contain
    #' all the information necessary for VisioProt-MS
  #' @param PSM_tab table PrSM from PD
  #' @param MSMS_tab table of MSMS from PD
  
  validate(need(sum(grepl("Master.Protein.Descriptions", names(PSM_tab))) == 1 | sum(grepl("Protein.Accessions", names(PSM_tab))) == 1, 
                paste0("Error in file format for plotting MS2 data:\n",
                       "You have uploaded a file from Proteome Discoverer.\n",
                       "One of the columns of the columns shoule be \"Master.Protein.Descriptions\" or \"Protein.Accessions\" in the PrSM table.")
  ))
  validate(need(sum(grepl("RT.in.min", names(MSMS_tab))) == 1 | sum(grepl("RT..min.", names(MSMS_tab))) == 1, 
                paste0("Error in file format for plotting MS2 data:\n",
                       "You have uploaded a file from Proteome Discoverer.\n",
                       "One of the columns of the columns shoule be \"RT.in.min\" or \"RT..min.\" in the MSMS table.")
  ))
  validate(need(sum(grepl("Spectrum.File", names(MSMS_tab))) == 1 & sum(grepl("First.Scan", names(PSM_tab))) == 1, 
                paste0("Error in file format for plotting MS2 data:\n",
                       "You have uploaded a file from Proteome Discoverer.\n",
                       "One of the columns of the columns shoule be \"Spectrum.File\" in the MSMS table and \"First.Scan\" in the PrSM table.")
  ))
}

# MSPathFinder: ----------------------------------------------------------------

MSPathFinder_MS2_check <- function(inputMSMS = NULL, inputMS = NULL) {
  validate(
    need(grepl("IcT", inputMSMS$name, fixed = T), 
         paste0("Error in file format for plotting MS2 data.\n",
                "You have to upload the \"IcTarget\" or \"IcTda\" output file from MSPathFinder",
                "associated with the deconvoluted MS weights uploaded as \"input file for MS\" (\".ms1ft\").")
         )
  )
  validate(
    need(grepl(".ms1ft", inputMS$name, fixed = T), 
         paste0("Error in file format for plotting MS2 data.\n",
                "You have to upload the \"IcTarget\" or \"IcTda\" output file from MSPathFinder",
                "associated with the deconvoluted MS weights uploaded as \"input file for MS\" (\".ms1ft\").")
    )
  )
  
  MS2PF <- read.table(inputMSMS$datapath, sep = "\t", header = F, skip = 1)
  
  # Validate that each MS feature has only one identification
  vec <- lapply(unique(MS2PF[,15]), function(x) {
    MS2PF[MS2PF[,15]==x,8]
  })
  vec <- lapply(vec, function(x) {
    length(unique(x))
  })
  vec <- unlist(vec)
  val <- max(vec)
  
  validate(
    need(val == 1, 
         "Several IDs (\"ProteinDesc\") have been attributed to the same MS1 feature in the MS2 output file."
         )
  )
}
