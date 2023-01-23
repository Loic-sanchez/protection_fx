################################################################################
# script : 01_effect_sizes.R
#
# Commentary : This script creates both the "OS_occu.RData" output file and 
#              "df_abun.RData" file in /outputs. It contains the model estimates
#               for occurrence and abundance data.
#
#             CAUTION : these two functions take a long time to run. Occurrences
#             took two days with 50 cores, abundances took 9 days with 50 cores.
#
################################################################################

library(tidyverse)
load(here::here("data", "data", "melted_filt.RData"))
load(here::here("data", "data", "melted_filt_ab.RData"))

source(here::here("R", "boot_effect_sizes.R"))

OS_occu = compute_ES_occu()
df_abun = compute_ES_abun()