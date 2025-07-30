## Test unidec output read
library(bit64)

# Advanced plotting libraries
# devtools::install_github('hadley/ggplot2')  # Use development version if needed
library(ggplot2)    # Grammar of graphics for static plots
library(plotly)     # Interactive web-based data visualization
library(dplyr)      # Data manipulation and transformation

# Visualization and UI enhancement packages
library(RColorBrewer)  # Color palettes for data visualization
library(shinyBS)       # Bootstrap components for Shiny (tooltips, modals, etc.)
library(data.table)    # High-performance data manipulation and reading

# Import file
library(rhdf5)  # For Unidec input compatibility


h5f = H5Fopen("files/test/test-unidec-visioprot/OFJMX190516_91.hdf5")



str(dat)   

h5readAttributes(file = "files/test/test-unidec-visioprot/OFJMX190516_91.hdf5", 
                 name = "/ms_dataset/30")
plot((1:length(H5Dread(h5f&'ms_dataset'&'30'&'mass_grid', bit64conversion = "bit64"))), H5Dread(h5f&'ms_dataset'&'30'&'mass_grid', bit64conversion = "bit64"), pch = ".")
H5Dread(h5f&'ms_dataset'&'mass_grid', bit64conversion = "bit64")
vec1 <- H5Dread(h5f&'ms_dataset'&'30'&'mass_data', bit64conversion = "bit64")[1,]
vec2 <- H5Dread(h5f&'ms_dataset'&'30'&'mass_data', bit64conversion = "bit64")[2,]
plot(vec1, vec2)
# plot(1:length(H5Dread(h5f&'ms_dataset'&'mass_grid', bit64conversion = "bit64")), H5Dread(h5f&'ms_dataset'&'mass_grid', bit64conversion = "bit64"), pch = ".")
H5Dread(h5f&'ms_dataset'&'3'&'processed_data', bit64conversion = "bit64")[1:10]
H5Dread(h5f&'ms_dataset'&'3'&'peaks', bit64conversion = "bit64")
H5Dread(h5f&'ms_dataset'&'mass_axis', bit64conversion = "bit64")
H5Dread(h5f&'ms_dataset'&'mass_data', bit64conversion = "bit64")

## action:

h5closeAll()
filepath <- "files/test/test-unidec-visioprot/OFJMX190516_91.hdf5"

### get all groups (spectra)

dat <- h5ls("files/test/test-unidec-visioprot/OFJMX190516_91.hdf5", recursive = 3)
all_md <- sort(unique(dat$group[dat$name == "mass_data"]))

# ### open object:
# h5f <- H5Fopen("files/test/test-unidec-visioprot/OFJMX190516_91.hdf5")

### get RT:

# h5readAttributes(file = filepath, 
# name = "/ms_dataset/30")$timemid

RT <- sapply(all_md, function(x) {
  h5readAttributes(file = filepath, 
                   name = x)$timemid
})

### get mass axis values:

Mass <- H5Dread(h5f&'ms_dataset'&'mass_axis', bit64conversion = "bit64")

### get intensities:

ltab <- lapply(seq_along(all_md), function(x) {
  tab <- h5read(filepath, name = paste0(all_md[x], '/mass_data'), bit64conversion = "bit64")
  tab <- as.data.frame(t(tab))
  names(tab) <- c("Mass", "Intensity")
  tab$RT <- RT[x]
  return(tab)
})
gtab <- rbindlist(ltab)

g <- ggplot(data = gtab[gtab$Intensity > 500000,], aes(x = RT, y = Mass, col = log10(Intensity))) +
  geom_point(alpha = 0.3) +
  scale_color_viridis_c(direction = -1) +
  theme_bw()
print(g)
