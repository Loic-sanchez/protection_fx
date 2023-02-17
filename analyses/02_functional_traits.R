################################################################################
# script : 02_functional_traits.R
#
# Commentary : This script uses both the "OS_occu.RData" output file and 
#              "df_abun.RData" file in /outputs. It computes mixed models to
#              check the effect of traits on the response to protection. 
#
#
################################################################################

load(here::here("outputs", "OS_occu.RData"))
load(here::here("outputs", "df_abun.RData"))
load(here::here("data", "data", "melted_filt.RData"))
load(here::here("data", "raw_data", "traits.RData"))

source(here::here("R", "traits_models.R"))
source(here::here("R", "plots.R"))
library(tidyverse)

rarity = compute_rarity()
taxonomy = extract_taxo()

occu_traits = join_occu_traits() # Join ES and traits 
abun_traits = join_abun_traits() # + Creates files necessary for plotting

occu_traits_both = join_rarity_taxo_occu() # Add rarity and taxonomy
abun_traits_both = join_rarity_taxo_abun() 

mod_occu_both = trait_model_occu() # Compute glmer
mod_abun_both = trait_model_abun() # + Creates files necessary for plotting

both_rsq() # [1] for occurrence, [2] for abundance
both_summaries() # [1] for occurrence, [2] for abundance
