# library(data.table)

source("files/Rfiles/Keep_highest_signal.R")
source("files/Rfiles/rename_input_coulumns.R")

InputFileMS <- data.frame(
  name = "OFABR230324_17_Job_2_results.csv",
  datapath = "files/Signal_integration/test_files/20250729/OFABR230324_17_Job_2_results.csv"
)

input <- data.frame(
  IntensityThresh = 20
)

## Currently done in app.:

i = 1
lfiles <- vector(mode = "list")

lfiles[[i]] <- read.table(InputFileMS[i, 'datapath'], sep = "\t", header = F)
names(lfiles) <- InputFileMS$datapath
lfiles <- lapply(lfiles, RenameRoWinPro)
# Add placeholder columns for consistency with BioPharma format
lfiles[[i]] <- cbind(lfiles[[i]][,1:3], " Temp1" = rep(NA, nrow(lfiles[[i]])), "Temp2" = rep(NA, nrow(lfiles[[i]])))

lfiles <- ThresholdCleaning(lfiles, input$IntensityThresh)

##--------- NEW TO IMPLEMENT -------------------------------------------------##

## new variables:
path_to_input <- "files/Signal_integration/test_files/20250729/VisioProt_Signal_integration_input_template_Amélie.txt"

## Read metadata table:
input_summarization <- read.table(path_to_input, sep = "\t", header = T)

## Add MW-min and MW-max:
input_summarization$MW_min <- input_summarization$MW - input_summarization$delta_MW
input_summarization$MW_max <- input_summarization$MW + input_summarization$delta_MW

## get sum of intensities:
tab_signal = lfiles[[i]]

indexes <- lapply(1:nrow(input_summarization), function(r) {
  which(tab_signal$RT <= input_summarization$RT_stop[r] &
          tab_signal$RT >= input_summarization$RT_start[r] & 
          tab_signal$Mass <= input_summarization$MW_max[r] &
          tab_signal$Mass >= input_summarization$MW_min[r]
  )
})

input_summarization$Sum_intensities <- sapply(indexes, function(r) {
  sum(tab_signal$intensity[r], na.rm = T)
})
input_summarization$Prop_sum_intensities <- sapply(indexes, function(r) {
  sum(tab_signal$intensity[r], na.rm = T) / sum(tab_signal$intensity, na.rm = T)
})
input_summarization$Point_count <- sapply(indexes, function(r) {
  length(tab_signal$intensity[r][!is.na(tab_signal$intensity[r])])
})
input_summarization$Percentage_top_signal_considered <- input$IntensityThresh

## export:

write.table(input_summarization, sep = "\t", row.names = F, file = paste0("VisioProt-MS_", substring(InputFileMS$name, first = 1, last = (nchar(InputFileMS$name)-4)), "_", Sys.Date(), ".tsv"))
