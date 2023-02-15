compute_rarity = function(){
  
  rarity = melted_filt |>
    group_by(Species) |> 
    summarise(mean_occu = mean(Presence))
  
  rarity$quantiles = as.character(ntile(rarity$mean_occu, 10))
  rarity = rarity %>% 
    mutate(Rarity = ifelse(quantiles > 2, "Common", "Rare"))
  
  save(rarity, file = here::here("outputs", "rarity.RData"))
  return(rarity)
}

extract_taxo = function(){
  
  taxonomy = rfishbase::load_taxa()
  return(taxonomy)
  
}

join_occu_traits = function(){

traits$Species = traits$CURRENT_TAXONOMIC_NAME
occu_traits = OS_occu %>% 
  filter(estimate_outside > 0.01) %>% 
  filter(estimate_full > 0.01) %>% 
  filter(estimate_part > 0.01) %>% 
  left_join(traits, by = "Species") %>% 
  mutate(RR_full = estimate_full/estimate_outside) %>% 
  mutate(RR_part = estimate_part/estimate_outside) 

save(occu_traits, file = here::here("outputs", "occu_traits.RData"))
return(occu_traits)

}

join_rarity_taxo_occu = function(){

traits$Species = traits$CURRENT_TAXONOMIC_NAME

occu_traits_both = occu_traits %>% 
  left_join(unique(taxonomy %>% select(Family, Order)), by = "Family") %>% 
  left_join(rarity, by = "Species") %>%
  pivot_longer(cols = c("RR_full", "RR_part")) %>% 
  rename(RR = "value", Protection = "name") %>% 
  # mutate(extremes = ntile(RR, 100)) %>%
  # filter(extremes < 99 & extremes > 2) %>%
  drop_na(Order)

return(occu_traits_both)

}

join_abun_traits = function(){

traits$Species = traits$CURRENT_TAXONOMIC_NAME

abun_traits = df_abun %>% 
  filter(estimate_outside > 0.01) %>%
  filter(estimate_full > 0.01) %>%
  filter(estimate_part > 0.01) %>%
  left_join(traits, by = "Species") 

save(abun_traits, file = here::here("outputs", "abun_traits.RData"))
return(abun_traits)

}

join_rarity_taxo_abun = function() {

traits$Species = traits$CURRENT_TAXONOMIC_NAME
  
abun_traits_both = abun_traits %>% 
  left_join(unique(taxonomy %>% select(Family, Order)), by = "Family") %>% 
  left_join(rarity, by = "Species") %>% 
  pivot_longer(cols = c("IRR_full", "IRR_part")) %>% 
  rename(IRR = "value", Protection = "name") %>% 
  # mutate(extremes = ntile(IRR, 100)) %>%
  # filter(extremes < 99 & extremes > 2) %>%
  drop_na(Order) %>% 
  drop_na(Rarity)

return(abun_traits_both)

}

trait_model_occu = function() {
  
  mod_occu_both = lme4::glmer(RR ~ poly(Trophic.Level, 2) * Protection + Rarity + (1|Order),
                              data = occu_traits_both,
                              family = Gamma(link = "inverse"))
  save(mod_occu_both, file = here::here("outputs", "mod_occu_both.RData"))
  return(mod_occu_both)
  
}

trait_model_abun = function() {
  
  mod_abun_both = lme4::glmer(IRR ~ poly(Trophic.Level, 2) * Protection + Rarity + (1|Order),
                              data = abun_traits_both,
                              family = Gamma(link = "inverse"))
  save(mod_abun_both, file = here::here("outputs", "mod_abun_both.RData"))
  return(mod_abun_both)
  
}

both_rsq = function(){ 
  
  p = predict(mod_occu_both, type = "response")
  cor_occu = cor(p, occu_traits_both$RR)^2
  
  p = predict(mod_abun_both, type = "response")
  cor_abun = cor(p, abun_traits_both$IRR)^2
  
  return(c(cor_occu, cor_abun))
  
}

both_summaries = function(){
  
  summary_occu = summary(mod_occu_both)
  summary_abun = summary(mod_abun_both)
  
  sum_list = list(summary_occu, summary_abun)
  return(sum_list)
  
}