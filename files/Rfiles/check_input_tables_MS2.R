PD_MS2_check <- function(PSM_tab = NULL, MSMS_tab = NULL) {
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

