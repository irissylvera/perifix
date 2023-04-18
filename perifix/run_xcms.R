# Setup ----
library(tidyverse)
library(xcms)
options(pillar.sigfig=7)

# dataset_version <- "FT2040"
# dataset_version <- "MS3000"

# output_folder <- paste0("made_data_", dataset_version, "/")
# mzML_files <- list.files(paste0(output_folder, "mzMLs/"), full.names=TRUE)
prefilter_versioned <- c(3, 1e6)
# 
# file_data <- data.frame(filename=mzML_files) %>%
#   mutate(samp_type=str_extract(filename, "Blk|Smp|Std|Poo")) %>%
#   mutate(depth=str_extract(filename, "25m|DCM|175m|15m")) %>%
#   mutate(depth=ifelse(is.na(depth), "", depth)) %>%
#   mutate(colid=factor(paste0(depth, samp_type), levels=c("Blk", "25mSmp", "DCMSmp", "175mSmp", "15mSmp", "Std", "Poo"))) %>%
#   mutate(col=alpha(c("red", "blue", "green", "purple", "blue", "black", "#008080"), 0.8)[colid]) %>%
#   mutate(lwd=c(2, 1, 1, 1, 1, 2)[colid])
# if(!dir.exists(output_folder))dir.create(output_folder)

# raw_files <- list.files('../../mzMLs/pos')
# raw_files <- as.data.frame(raw_files)

fls <- list.files('../../mzMLs/pos', full.names = TRUE)
fls <- as.data.frame(fls) %>% 
  filter(!str_detect(fls, "DDA"))

# Creating metadata ----
metadataframe <- fls %>%
  # Grab just the unique filenames and rename the column
  distinct(filename=`fls`) %>%
  # Create a new column with sample type information
  mutate(samp_type=str_extract(filename, "Blk|Smp|Std|Poo")) %>%
  # Create a new column with timepoint information (either Tfinal or T0)
  mutate(timepoint=str_extract(filename, "Tfinal|T0")) %>%
  # Extract treatment information from the filename
  # Longest one needs to go first!
  mutate(treatment=str_extract(filename, "NPF|NF|PF|NP|N|P|F|C")) %>%
  # Replace accidental "P" treaments from "Pooled" with NAs
  mutate(treatment=ifelse(samp_type=="Poo", NA, treatment)) %>% 
  mutate(colid = str_extract(filename, "Std|Poo|Blk|NPF|NF|PF|NP|N|P|F|C|Tote"))
  # mutate(colid=(levels=c("Std", "Poo", "Blk", "NPF", "NF", "PF", "NP", "N", "P", "F", "C", "Tote")))

# XCMS things ----
register(BPPARAM = SerialParam(progressbar = TRUE))
msnexp <- readMSData(
  files = metadataframe$filename, 
  pdata = new("NAnnotatedDataFrame", metadataframe), 
  msLevel. = 1, 
  mode = "onDisk"
)


# saveRDS(msnexp, file = "msnexp.rds")
# msnexp <- readRDS("msnexp.rds")

register(BPPARAM = SnowParam(workers = 3, tasks = nrow(metadataframe), progressbar = TRUE))
cwp <- CentWaveParam(
  ppm = 5, 
  peakwidth = c(20, 80), 
  prefilter = prefilter_versioned, 
  snthresh = 0, 
  verboseColumns = TRUE, 
  extendLengthMSW = TRUE, 
  integrate = 2
)
msnexp_withpeaks <- findChromPeaks(msnexp, cwp)

register(BPPARAM = SerialParam(progressbar = TRUE))
obp <- ObiwarpParam(
  binSize = 0.1, 
  centerSample = round(nrow(metadataframe)/2), 
  response = 1, 
  distFun = "cor_opt"
)
msnexp_rtcor <- adjustRtime(msnexp_withpeaks, obp)

pdp <- PeakDensityParam(
  sampleGroups = metadataframe$colid, 
  bw = 12, 
  minFraction = 0.1, 
  binSize = 0.001, 
  minSamples = 2
)
msnexp_grouped <- groupChromPeaks(msnexp_rtcor, pdp)

fpp <- FillChromPeaksParam(ppm = 2.5)
msnexp_filled <- fillChromPeaks(msnexp_grouped, fpp)

saveRDS(msnexp_withpeaks, file = "msnexp_withpeaks.rds")
# msnexp <- readRDS("msnexp.rds")
saveRDS(msnexp_rtcor, file = "msnexp_rtcor.rds")
saveRDS(msnexp_grouped, file = "msnexp_grouped.rds")
saveRDS(msnexp_filled, file = "msnexp_filled.rds")
