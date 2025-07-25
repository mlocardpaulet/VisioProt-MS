PD_MS2_check <- function(PSM_tab, MSMS_tab) {
  validate(need(sum(grepl("Master.Protein.Descriptions", names(PSM_tab))) == 1 | sum(grepl("Protein.Accessions", names(PSM_tab))) == 1, 
         "Error in file format for plotting MS2 data.\nYou have to upload the following files:\n- A MSMSSpectrumInfo.txt file from BioPharma Finder (in the \"MS/MS File\" field).\n- The corresponding PSMs.txt or PrSMs.txt file (in the \"PSM File\" field).\n\nOne of the columns of the columns shoule be \"Master.Protein.Descriptions\" or \"Protein.Accessions\" in the PrSM table.")
  )
  validate(need(sum(grepl("RT.in.min", names(MSMS_tab))) == 1 | sum(grepl("RT..min.", names(MSMS_tab))) == 1, 
                "Error in file format for plotting MS2 data.\nYou have to upload the following files:\n- A MSMSSpectrumInfo.txt file from BioPharma Finder (in the \"MS/MS File\" field).\n- The corresponding PSMs.txt or PrSMs.txt file (in the \"PSM File\" field).\n\nOne of the columns of the columns shoule be \"RT.in.min\" or \"RT..min.\" in the MSMS table.")
  )
  validate(need(sum(grepl("Spectrum.File", names(MSMS_tab))) == 1 | sum(grepl("First.Scan", names(PSM_tab))) == 1, 
                "Error in file format for plotting MS2 data.\nYou have to upload the following files:\n- A MSMSSpectrumInfo.txt file from BioPharma Finder (in the \"MS/MS File\" field).\n- The corresponding PSMs.txt or PrSMs.txt file (in the \"PSM File\" field).\n\nOne of the columns of the columns shoule be \"Spectrum.File\" in the MSMS table and \"First.Scan\" in the Prsm table.")
  )
}

