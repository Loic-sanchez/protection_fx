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

occu_traits$tranche[occu_traits$Trophic.Level < 2.5] = "Low"
occu_traits$tranche[occu_traits$Trophic.Level > 3.5] = "High"
occu_traits$tranche[occu_traits$Trophic.Level >= 2.5 & occu_traits$Trophic.Level <= 3.5] = "Mid"
occu_traits$tranche = factor(occu_traits$tranche, levels = c("Low", "Mid", "High"))

occu_traits = occu_traits %>% drop_na(Trophic.Level)
abun_traits = abun_traits %>% drop_na(Trophic.Level)

occu_traits$sizeclass = cut(occu_traits$MaxLength, breaks=c(quantile(occu_traits$MaxLength, probs = seq(0, 1, by = 1/4))), 
    labels=c("low","midlow","midhigh", "high"))

occu_over = occu_traits[occu_traits$estimate_outside < 0.5,]
table(occu_traits$tranche, occu_traits$sizeclass)
ggplot(occu_traits, aes(x = tranche, y = MaxLength)) +
  geom_violin(aes(fill = tranche), show.legend = F) +
  theme_bw()
