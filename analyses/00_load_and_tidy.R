################################################################################
# script : 00_load_and_tidy.R
#
# Commentary :This part can be jumped and you can start at 01; by loading 
#             the file in "data", "data" called "melted_filt.RData"
#
################################################################################

devtools::load_all()

library(tidyverse)
source(here::here("R", "load_data.R"))
source(here::here("R", "tidy_data.R"))

load_rdata()

covariates_corrected = tidy_sites()

ssmatrix = make_ssmatrix(30)

all_PC = make_PCA()

melted = melt_and_merge()

melted = recode_protection()

melted_filt = rid_rare_outside()

##################################

#all species targeted in the family:

all_sp <- c("Acanthuridae", "Caesionidae", "Carangidae", "Ephippidae", "Haemulidae", "Kyphosidae",
            "Labridae", "Lethrinidae", "Lutjanidae", "Mullidae", "Nemipteridae", "Scaridae",
            "Scombridae", "Serranidae", "Siganidae", "Sparidae", "Sphyraenidae")

#families with species targeted if larger than 20cm:

target_larger_20cm <- c("Balistidae", "Holocentridae", "Pomacanthidae", "Priacanthidae")
traits = traits %>%
  mutate(target = ifelse(Family %in% all_sp | (Family %in% target_larger_20cm & MaxLength > 19.9), "Targeted", "Not Targeted"))
traits$target[is.na(traits$target)] = "Not Targeted"
traits$Species = traits$CURRENT_SPECIES_NAME
traits$target = as.factor(traits$target)

occu_traits = OS_occu %>% 
  filter(estimate_outside > 0.01) %>% 
  filter(estimate_full > 0.01) %>% 
  filter(estimate_part > 0.01) %>% 
  left_join(traits, by = "Species") %>% 
  mutate(RR_full = estimate_full/estimate_outside) %>% 
  mutate(RR_part = estimate_part/estimate_outside) %>% 
  drop_na(RR_full)

abun_traits = df_abun %>% 
  filter(estimate_outside > 0.01) %>% 
  filter(estimate_full > 0.01) %>% 
  filter(estimate_part > 0.01) %>% 
  left_join(traits, by = "Species") %>% 
  drop_na(IRR_full)

p = ggplot(occu_traits, aes(x = RR_part, y = RR_full, color = AUC, alpha = 0.5)) + 
  geom_point(show.legend = F) +
  theme_bw()

ggExtra::ggMarginal(p, type = "histogram", 
           fill = "lightblue", 
           xparams = list(bins = 50), 
           yparams =list(bins = 50))

p = ggplot(abun_traits, aes(x = IRR_part, y = IRR_full, color = RSQ, alpha = 0.5)) + 
  geom_point(show.legend = F) +
  theme_bw()

ggExtra::ggMarginal(p, type = "histogram", 
           fill = "lightblue", 
           xparams = list(bins = 50), 
           yparams =list(bins = 50))

abunoccu = occu_traits %>% 
  left_join(abun_traits %>% select(Species, IRR_full, IRR_part), by = "Species")

ggplot(abunoccu, aes(x = RR_full, y = IRR_full)) + 
  geom_point(colour = "#007FFF", show.legend = F, alpha = 0.5) +
  theme_bw()

ggplot(abunoccu, aes(x = RR_part, y = IRR_part)) + 
  geom_point(colour = "#007FFF", show.legend = F) +
  theme_bw()

occu_traits$Water.column = as.factor(fct_collapse(occu_traits$Water.column, "Pelagic" = c("pelagic non-site attached", "pelagic site attached")))
occu_traits$Trophic.group2 = as.factor(fct_collapse(occu_traits$Trophic.group2, "Benthic invertivore" = c("Benthic invertivore", "Corallivore")))
abun_traits$Trophic.group2 = as.factor(fct_collapse(abun_traits$Trophic.group2, "Benthic invertivore" = c("Benthic invertivore", "Corallivore")))
abun_traits$Water.column = as.factor(fct_collapse(abun_traits$Water.column, "Pelagic" = c("pelagic non-site attached", "pelagic site attached")))

taxonomy = rfishbase::load_taxa()
rarity$quantiles = as.character(ntile(rarity$mean_occu, 10))
rarity = rarity %>% 
  mutate(Rarity = ifelse(quantiles > 2, "Common", "Rare"))

occu_traits = occu_traits %>% 
  left_join(taxonomy %>% select(Species, Order), by = "Species") %>% 
  left_join(rarity, by = "Species") %>% 
  # mutate(extremes = ntile(RR_full, 100)) %>% 
  # filter(extremes < 100) %>% 
  drop_na(Order)

mod = lme4::glmer(RR_full ~ poly(Trophic.Level, 2) * Rarity + (1|Order),
          data = occu_traits,
          family = Gamma(link = "inverse"))

car::Anova(mod, type = 2, test.statistic = "F")
anova(mod)

plot(mod)
hist(residuals(mod), breaks = 50)
res = DHARMa::simulateResiduals(mod)
plot(res)

summary(mod)

p = predict(mod, type = "response")
cor(p, occu_traits$RR_full)^2

abun_traits = abun_traits %>% 
  left_join(taxonomy %>% select(Species, Order), by = "Species") %>% 
  left_join(rarity, by = "Species") %>% 
  # mutate(extremes = ntile(IRR_full, 100)) %>% 
  # filter(extremes < 100) %>% 
  drop_na(Order)

mod_abun = lme4::glmer(IRR_full ~ poly(Trophic.Level, 2) * Rarity + (1|Order),
               data = abun_traits,
               family = Gamma(link = "inverse"))

par(mfrow = c(1,2))
plot(mod_abun)

hist(residuals(mod_abun), breaks = 30)
summary(mod_abun)

p = predict(mod_abun, type = "response")
zbl = abun_traits$IRR_full[as.numeric(intersect(names(p), rownames(abun_traits)))]
cor(p, zbl)^2

par(mfrow = c(1,2))
p_occu = visreg::visreg(mod, 
               xvar = "Trophic.Level", 
               type = "conditional", 
               scale = "response", 
               by = "Rarity", 
               overlay = T, 
               rug = F,
               gg = T,
               line = list(lty = 1), 
               legend = F)

p_abun = visreg::visreg(mod_abun, 
               xvar = "Trophic.Level", 
               type = "conditional", 
               scale = "response", 
               by = "Rarity", 
               overlay = T,
               rug = F,
               gg = T, 
               line = list(lty =2), 
               legend = F)

theme_set(new = theme_bw())
p_occu + 
  p_abun$layers +
  xlab("Trophic level") +
  ylab("Effect Size of FPAs") +
  theme(legend.position = "none")

