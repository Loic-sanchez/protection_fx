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

covariates_corrected = tidy_sites()

# Filter minimum occurrences per species + make a site-species matrix

make_ssmatrix = function(n) {
  
  fish_nospp <- fish[fish$TAXONOMIC_NAME %in% fish$TAXONOMIC_NAME,] |> 
    filter(!str_detect(TAXONOMIC_NAME,"spp.")) |> 
    filter(!str_detect(TAXONOMIC_NAME,"sp.")) |> 
    filter(!str_detect(TAXONOMIC_NAME,"(cf)")) # Keep only species levels
  
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

ssmatrix = make_ssmatrix(30)

# PCA

make_PCA = function() {
  
  all_cov = hab_filt %>% 
    left_join(socio_filt, by = "SurveyID") %>%
    left_join(env_filt, by = "SurveyID") %>%
    left_join(fine_habitat[,-c(2:4)], by = "SurveyID")
  
  ACP = FactoMineR::PCA(all_cov[,-1], ncp = 25, scale.unit = T)
  # factoextra::fviz_eig(ACP, ncp = 25)
  # factoextra::get_eig(ACP)
  PC_coords = ACP$ind$coord
  
  all_PC = cbind(all_cov[,1], PC_coords[,1:15])
  all_PC = as.data.frame(all_PC) %>%
    rename(SurveyID = "V1")
  
  return(all_PC)
}

all_PC = make_PCA()

# Melt matrix and merge with PC

melt_and_merge = function() {
  
  melted = ssmatrix %>%
    reshape2::melt() %>%
    rename(SurveyID = "Var1", Species = "Var2", Presence = "value") %>%
    left_join(all_PC, by = "SurveyID") %>%
    left_join(covariates_corrected[,-c(2:25,27)], by = "SurveyID") 
  
  melted = melted %>%
    drop_na(Dim.1)
  
  return(melted)
  
}

melted = melt_and_merge()

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

melted = recode_protection()

