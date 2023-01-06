load(here::here("outputs", "df_abund.RData"))
load(here::here("outputs", "df_vif.RData"))

library(tidyverse)

abund_vif = df_abund |>
  dplyr::left_join(df_vif |> dplyr::select(Species, med_maxVIF), by = "Species") |> 
  tidyr::drop_na(med_maxVIF) |> 
  dplyr::mutate(ESF = log(med_full/med_out), ESP = log(med_part/med_out))

abund_lowvif = abund_vif %>% 
  filter(med_maxVIF < 2.24)

plot(med_maxVIF ~ med_full, data = abund_lowvif)

abund_lowvif  %>% 
  drop_na(ESF) %>% 
  ggplot(aes(x = ESF)) +
  geom_histogram(color="black", fill="white", binwidth = 0.15, stat = "bin") +
  facet_wrap(vars(freq)) +
  theme_light()

abund_lowvif  %>% 
  drop_na(ESF) %>% 
  ggplot(aes(x = ESF)) +
  geom_histogram(color="black", fill="white", binwidth = 0.15, stat = "bin") +
  theme_light()

load(here::here("data", "raw_data", "traits.RData"))
traits$Species = traits$CURRENT_SPECIES_NAME

abund_n_traits = abund_lowvif |>
  left_join(traits, by = "Species")

ggplot(abund_n_traits, aes(x = log10(MaxLength), y = ESF)) +
  geom_point()

family_numb = abund_n_traits |>
  group_by(Family) |>
  summarize(number = n()) |>
  filter(number > 9)

ggplot(abund_n_traits[abund_n_traits$Family %in% family_numb$Family,], aes(x = Family, y = ESF)) +
  geom_boxplot() +
  theme_light()
