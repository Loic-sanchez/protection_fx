tidy_sites = function() {
  
  sites_correct$Effectiveness_cor = "out"  #Set to outside all sites that were unprotected when surveyed
  
  covariates_corrected = dplyr::left_join(covariates, 
                                   sites_correct[,c(1,21)], 
                                   by = "SurveyID")
  
  covariates_corrected$Effectiveness = replace(covariates_corrected$Effectiveness, 
                                               covariates_corrected$Effectiveness_cor == "out",
                                               "out") # Correct it on the original covariates data
  
  return(covariates_corrected)
}

# Filter minimum occurrences per species + make a site-species matrix

make_ssmatrix = function(n) {
  
  fish_nospp <- fish[fish$TAXONOMIC_NAME %in% fish$TAXONOMIC_NAME,] |> 
    filter(!stringr::str_detect(TAXONOMIC_NAME,"spp.")) |> 
    filter(!stringr::str_detect(TAXONOMIC_NAME,"sp.")) |> 
    filter(!stringr::str_detect(TAXONOMIC_NAME,"(cf)")) # Keep only species levels
  
  tab <- table(fish_nospp$SPECIES_NAME)
  boxplot(tab, outline = F)
  tab <- tab[tab>n]
  fish_filt  <- fish_nospp[fish_nospp$SPECIES_NAME %in% names(tab),]
  fish_filt$SPECIES_NAME = as.factor(fish_filt$SPECIES_NAME)
  
  ssmatrix = t(fossil::create.matrix(fish_filt, 
                             tax.name = "SPECIES_NAME", 
                             locality = "SurveyID",
                             abund = F))
  
  return(ssmatrix)
}

# PCA

make_PCA = function() {
  
  all_cov = hab_filt |>
    left_join(socio_filt, by = "SurveyID") |>
    left_join(env_filt, by = "SurveyID") |>
    left_join(fine_habitat[,-c(2:4)], by = "SurveyID")
  
  ACP = FactoMineR::PCA(all_cov[,-1], ncp = 25, scale.unit = T)
  # factoextra::fviz_eig(ACP, ncp = 25)
  # factoextra::get_eig(ACP)
  PC_coords = ACP$ind$coord
  
  all_PC = cbind(all_cov[,1], PC_coords[,1:15])
  all_PC = as.data.frame(all_PC) |>
    rename(SurveyID = "V1")
  
  return(all_PC)
}

melt_and_merge = function() { # Melt matrix and merge with PC
  
  melted = ssmatrix |>
    reshape2::melt() |>
    rename(SurveyID = "Var1", Species = "Var2", Presence = "value") |>
    left_join(all_PC, by = "SurveyID") |>
    left_join(covariates_corrected[,-c(2:25,27)], by = "SurveyID") |>
    left_join(sites_info[,c(1:4,8,10)], by = "SurveyID") 

  
  melted = melted |>
    drop_na(Dim.1) 
    
  
  return(melted)
  
}

recode_protection = function() {
  
  # Recode protection levels as Grorud-Colvert
  
  melted$Effectiveness = fct_collapse(melted$Effectiveness, 
                                      "Fully Protected" = c("High No take", "High No take multizoned"),
                                      "Highly Protected" = c("Medium No take", "Medium No take multizoned"),
                                      "Lightly Protected" = c("High Restricted take", "High Restricted take multizoned", 
                                                              "Low No take", "Low No take multizoned", 
                                                              "Medium Restricted take", "Medium Restricted take multizoned"),
                                      "Minimally Protected" = c("Low Restricted take","Low Restricted take multizoned"),
                                      "Unprotected" = c("Medium Fishing", "Low Fishing", "out"))
  
  melted$Effectiveness = as.factor(as.character(
    plyr::revalue(melted$Effectiveness, c("Fully Protected" = "Full Protection",
                                          "Highly Protected" = "Partial Protection",
                                          "Lightly Protected" = "Partial Protection",
                                          "Minimally Protected" = "Partial Protection",
                                          "Unprotected" = "Unprotected"))))
  
  return(melted)
  
}

rid_rare_outside = function() { # Get rid of species with 0-4 occurrences in the reference level

  # Some values were dropped, make sure we still have at least N occurrences
  
  melted_sum = melted |>
    group_by(Species) |>
    summarise(Presence = sum(Presence)) |>
    dplyr::filter(Presence >= 30)
  
  spec_vec = as.character(melted_sum$Species)
  
  melted_filt = melted[melted$Species %in% spec_vec,]
  
  # Get rid of rares outside
  
  melted_fsum = melted_filt |>
    filter(Effectiveness == "Unprotected") |>
    group_by(Species) |>
    summarise(Presence = sum(Presence)) |>
    filter(Presence > 4) 

  spec_vec = as.character(melted_fsum$Species)

  melted_filt = melted_filt[melted_filt$Species %in% spec_vec,]
  melted_filt$Species = as.factor(as.character(melted_filt$Species))
  
  save(melted_filt, file = here::here("data", "data", "melted_filt.RData"))
  return(melted_filt)
}

