################################################################################
# script : 00_load_and_tidy.R
#
# Commentary :This part can be jumped and you can start at 01; by loading 
#             the file in "data", "data" called "melted_filt.RData"
#
################################################################################

library(tidyverse)
source(here::here("R", "load_data.R"))
source(here::here("R", "tidy_data.R"))

load_rdata()

covariates_corrected = tidy_sites()

ssmatrix = make_ssmatrix(30)
ssmatrix_ab = make_ssmatrix_ab(30)

all_PC = make_PCA()

melted = melt_and_merge()
melted_ab = melt_and_merge_ab()

melted = recode_protection()
melted_ab = recode_protection_ab()

melted_filt = rid_rare_outside() #also saves "melted_filt.RData", in data/data
melted_filt_ab = rid_rare_outside_ab() #also saves "melted_filt_ab.RData"

##################################

#### Plot crastination ####

# To delete later

### Plots richness ###

melted_filt_traits$Trophic.group2 = as.factor(melted_filt_traits$Trophic.group2)
melted_filt_traits$SurveyID = as.factor(melted_filt_traits$SurveyID)
mft_higher = melted_filt_traits %>% filter(Trophic.group2 == "Higher carnivore") %>% filter(Presence == 1)
mft_bentinv = melted_filt_traits %>% filter(Trophic.group2 == "Benthic invertivore") %>% filter(Presence == 1)
mft_herb = melted_filt_traits %>% filter(Trophic.group2 == "Herbivore") %>% filter(Presence == 1)
mft_plank = melted_filt_traits %>% filter(Trophic.group2 == "Planktivore") %>% filter(Presence == 1)

richness_higher = mft_higher %>% group_by(SurveyID) %>% summarise(richness_high = length(unique(Species)))
richness_bent = mft_bentinv %>% group_by(SurveyID) %>% summarise(richness_bent = length(unique(Species)))
richness_herb = mft_herb %>% group_by(SurveyID) %>% summarise(richness_herb = length(unique(Species)))
richness_plank = mft_plank %>% group_by(SurveyID) %>% summarise(richness_plank = length(unique(Species)))

zbl = melted_filt_traits %>% group_by(SurveyID) %>% summarise(Protection = unique(Effectiveness))

rich_zbl = zbl %>% 
  left_join(richness_higher, by = "SurveyID") %>% 
  left_join(richness_bent, by = "SurveyID") %>% 
  left_join(richness_herb, by = "SurveyID") %>% 
  left_join(richness_plank, by = "SurveyID")
rich_zbl[is.na(rich_zbl)] = 0

rich_zbl_long = rich_zbl %>% pivot_longer(cols = c("richness_high", 
                                                   "richness_bent", 
                                                   "richness_plank", 
                                                   "richness_herb"))

ggplot(rich_zbl_long, aes(x = Protection, y = value, fill = Protection)) +
  geom_boxplot() +
  facet_wrap(~name) +
  theme_bw()

### Plot with trophic levels ###

melted_filt_traits$SurveyID = as.factor(melted_filt_traits$SurveyID)
mft_high = melted_filt_traits %>% filter(Trophic.Level > 3.5) %>% filter(Presence == 1)
mft_mid = melted_filt_traits %>% filter(Trophic.Level < 2.5) %>% filter(Presence == 1)
mft_low = melted_filt_traits %>% filter(Trophic.Level >= 2.5 & Trophic.Level <= 3.5) %>% filter(Presence == 1)

richness_high = mft_high %>% group_by(SurveyID) %>% summarise(richness_high = length(unique(Species)))
richness_mid = mft_mid %>% group_by(SurveyID) %>% summarise(richness_mid = length(unique(Species)))
richness_low = mft_low %>% group_by(SurveyID) %>% summarise(richness_low = length(unique(Species)))

zbl = melted_filt_traits %>% group_by(SurveyID) %>% summarise(Protection = unique(Effectiveness))

rich_zbl = zbl %>% 
  left_join(richness_high, by = "SurveyID") %>% 
  left_join(richness_mid, by = "SurveyID") %>% 
  left_join(richness_low, by = "SurveyID") 
rich_zbl[is.na(rich_zbl)] = 0

rich_zbl_long = rich_zbl %>% pivot_longer(cols = c("richness_high", 
                                                   "richness_mid", 
                                                   "richness_low"))
rich_zbl_long$name = factor(rich_zbl_long$name, levels = c("richness_low", "richness_mid", "richness_high"))
ggplot(rich_zbl_long, aes(x = Protection, y = value, fill = Protection)) +
  geom_violin(draw_quantiles = 0.5) +
  facet_wrap(~name) +
  theme_bw()

### Estimates in/out ###

occu_traits_long = occu_traits %>% 
  pivot_longer(cols= c(estimate_full, estimate_outside, estimate_part))
occu_traits_long$name = factor(occu_traits_long$name, levels = c("estimate_outside", "estimate_part", "estimate_full"))

ggplot(occu_traits_long, aes(x = name, y = value, fill = name)) +
  geom_violin(draw_quantiles = 0.5) +
  facet_wrap(~Trophic.group2) +
  theme_bw()

### Raw ratio/trophic level

ggplot(occu_traits_both, aes(x = Trophic.Level, y = RR)) +
  geom_point() +
  facet_wrap(~Protection) +
  theme_bw()
