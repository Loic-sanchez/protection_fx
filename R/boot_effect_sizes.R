compute_ES_occu = function() {
  
  all_sp_glm_qual = parallel::mclapply(1:length(unique(melted_filt$Species)), mc.cores = 50, function(i) {
    
    subdf = filter(melted_filt, Species == levels(melted_filt$Species)[i])
    
    pres_sp1 = subdf %>%
      filter(Presence > 0)
    
    pres_sp1$Ecoregion = as.factor(as.character(pres_sp1$Ecoregion))
    
    eco_table = table(pres_sp1$Ecoregion)
    eco_tab = names(eco_table[eco_table > 9])
    
    abs_sp1 = subdf %>%
      filter(Presence < 1)
    
    abs_sp1 = abs_sp1[abs_sp1$Ecoregion %in% eco_tab,]
    pres_sp1 = pres_sp1[pres_sp1$Ecoregion %in% eco_tab,]
    
    presabs_sp1 = rbind(pres_sp1, abs_sp1)
    
    istherefull = ifelse(sum(presabs_sp1$Presence[which(presabs_sp1$Effectiveness == "Full Protection")]) > 4, T, F)
    istherepart = ifelse(sum(presabs_sp1$Presence[which(presabs_sp1$Effectiveness == "Partial Protection")]) > 4, T, F)
    isthereunp = ifelse(sum(presabs_sp1$Presence[which(presabs_sp1$Effectiveness == "Unprotected")]) > 4, T, F)
    
    presabs_sp1_nofull = presabs_sp1 %>%
      filter(Effectiveness != "Full Protection")
    presabs_sp1_nopart = presabs_sp1 %>%
      filter(Effectiveness != "Partial Protection")
    
    if(istherefull == F){presabs_sp1 = presabs_sp1_nofull}
    if(istherepart == F){presabs_sp1 = presabs_sp1_nopart}
    
    if(isthereunp == F){return(NULL)} 
    if(istherefull == F & istherepart == F){return(NULL)}
    
    names_var = names(presabs_sp1[, c(4:12)])
    formula_var = paste0(names_var, collapse = "+")
    full_formula = paste("Presence ~ ", formula_var, sep = "")
    
    model = glm(formula = full_formula,
                family = "binomial",
                data = presabs_sp1)
    
    if(model$converged == F){return(NULL)}
    
    if(istherefull == F | istherepart == F){vif = sqrt((regclass::VIF(model)))} else{vif = regclass::VIF(model)[, 3]}
    
    while(max(vif) > 2.23){ 
      
      u = which.max(vif)
      names_var = names_var[-u]
      formula_var = paste0(names_var, collapse = "+")
      full_formula = paste("Presence ~ ", formula_var, sep = "")
      
      model = glm(formula = full_formula,
                  family = "binomial",
                  data = presabs_sp1)            
      
      if(istherefull == F | istherepart == F){vif = sqrt((regclass::VIF(model)))} else{vif = regclass::VIF(model)[, 3]}
      
      full_formula = full_formula
    }
    
    full_formula_sp = as.formula(paste(full_formula, "+ Matern(1 | SiteLongitude + SiteLatitude)", sep = ""))
    separ = spaMM::is_separated.formula(formula = full_formula_sp,
                                        data = presabs_sp1)
    
    if(separ == T){return(NULL)} 
    
    sp_model = spaMM::fitme(full_formula_sp,
                            family = "binomial",
                            method = "PQL/L",
                            data = presabs_sp1)
    
    mod_coefs = as.data.frame(sp_model$fixef)
    intercept = mod_coefs["(Intercept)",]
    
    estimate_outside = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"Intercept")) == 1,
                              boot::inv.logit(intercept),
                              NA)
    
    estimate_full = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessFull Protection")) == 1,
                           boot::inv.logit(intercept+mod_coefs["EffectivenessFull Protection",]),
                           NA)
    
    estimate_part = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessPartial Protection")) == 1,
                           boot::inv.logit(intercept+mod_coefs["EffectivenessPartial Protection",]),
                           NA)
    
    OR_full = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessFull Protection")) == 1,
                     exp(mod_coefs["EffectivenessFull Protection",]),
                     NA)
    
    OR_part = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessPartial Protection")) == 1,
                     exp(mod_coefs["EffectivenessPartial Protection",]),
                     NA)
    
    p = spaMM::predict.HLfit(sp_model, type = "response")
    roc.mod = pROC::roc(presabs_sp1$Presence, as.numeric(p))
    AUC = pROC::auc(roc.mod)
    
    Species = levels(melted_filt$Species)[i]
    
    coef = data.frame(Species,
                      estimate_outside,
                      estimate_full, 
                      estimate_part,
                      OR_full,
                      OR_part,
                      AUC)
    
    # cat(i, "\n")
    return(coef)
    
  }) 
  
  OS_occu = data.table::rbindlist(all_sp_glm_qual)
  save(OS_occu, file = here::here("outputs", "OS_occu.RData"))
  return(OS_occu)
  
}

compute_ES_abun = function() {
 
  melted_zero = melted_filt_ab |> 
    filter(Abundance < 1)
  melted_filt = melted_filt_ab |> 
    filter(Abundance > 0)
  
  melted_filt = melted_filt[, -c(12:18)]
  melted_zero = melted_zero[, -c(12:18)]
  
  all_sp_glm_abund = parallel::mclapply(1:length(unique(melted_filt$Species)), mc.cores = 50, function(i) {
    
    subdf = filter(melted_filt, Species == levels(melted_filt$Species)[i])
    subdf_zero = filter(melted_zero, Species == levels(melted_filt$Species)[i])
    
    pres_sp1 = subdf
    abs_sp1 = subdf_zero
    
    pres_sp1$Ecoregion = as.factor(as.character(pres_sp1$Ecoregion))
    
    eco_table = table(pres_sp1$Ecoregion)
    eco_tab = names(eco_table[eco_table > 9])
    
    abs_sp1 = abs_sp1[abs_sp1$Ecoregion %in% eco_tab,]
    pres_sp1 = pres_sp1[pres_sp1$Ecoregion %in% eco_tab,]
    
    presabs_sp1 = rbind(pres_sp1, abs_sp1)
    
    presabs_sp1_nofull = presabs_sp1 %>%
      filter(Effectiveness != "Full Protection")
    presabs_sp1_nopart = presabs_sp1 %>%
      filter(Effectiveness != "Partial Protection")
    
    istherefull = ifelse(length(presabs_sp1$Abundance[which(presabs_sp1$Effectiveness == "Full Protection")]) > 4, T, F) 
    istherepart = ifelse(length(presabs_sp1$Abundance[which(presabs_sp1$Effectiveness == "Partial Protection")]) > 4, T, F)
    isthereunp = ifelse(length(presabs_sp1$Abundance[which(presabs_sp1$Effectiveness == "Unprotected")]) > 4, T, F)
    
    if(istherefull == F){presabs_sp1 = presabs_sp1_nofull}
    if(istherepart == F){presabs_sp1 = presabs_sp1_nopart}
    
    if(isthereunp == F){return(NULL)} 
    if(istherefull == F & istherepart == F){return(NULL)}
    
    names_var = names(presabs_sp1[, c(4:12)])
    formula_var = paste0(names_var, collapse = "+")
    full_formula = paste("Abundance ~ ", formula_var, sep = "")
    
    tryCatch({model = MASS::glm.nb(formula = full_formula,
                                   data = presabs_sp1)}, error = function(e) model <<- NULL)
    
    if(is.null(model) == T){return(NULL)}
    if(istherefull == F | istherepart == F){vif = sqrt((regclass::VIF(model)))} else{vif = regclass::VIF(model)[, 3]}
    
    tryCatch({
      while(max(vif) > 2.23){ 
        
        u = which.max(vif)
        names_var = names_var[-u]
        formula_var = paste0(names_var, collapse = "+")
        full_formula = paste("Abundance ~ ", formula_var, sep = "")
        
        model = MASS::glm.nb(formula = full_formula,
                             data = presabs_sp1)
        
        if(istherefull == F | istherepart == F){vif = sqrt((regclass::VIF(model)))} else{vif = regclass::VIF(model)[, 3]}
        
        full_formula = full_formula
      }}, error = function(e) return(NULL))
    
    full_formula_sp = as.formula(paste(full_formula, "+ Matern(1 | SiteLongitude + SiteLatitude)", sep = ""))
    
    tryCatch({sp_model = spaMM::fitme(full_formula_sp,
                                      family = "negbin",
                                      method = c("ML", "exp"),
                                      data = presabs_sp1)}, error = function(e) sp_model <<- NULL)
    
    if(is.null(sp_model) == F){mod_coefs = as.data.frame(sp_model$fixef)}
    if(is.null(sp_model) == F){intercept = mod_coefs["(Intercept)",]}
    
    estimate_outside = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"Intercept")) == 1,
                              exp(intercept),
                              NA)
    
    estimate_full = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessFull Protection")) == 1,
                           exp(intercept+mod_coefs["EffectivenessFull Protection",]),
                           NA)
    
    estimate_part = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessPartial Protection")) == 1,
                           exp(intercept+mod_coefs["EffectivenessPartial Protection",]),
                           NA)
    
    IRR_full = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessFull Protection")) == 1,
                      exp(mod_coefs["EffectivenessFull Protection",]),
                      NA)
    
    IRR_part = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessPartial Protection")) == 1,
                      exp(mod_coefs["EffectivenessPartial Protection",]),
                      NA)
    
    tryCatch({p = spaMM::predict.HLfit(sp_model, type = "response")}, error = function(e) p <<- NA)
    tryCatch({RSQ = (cor(p, presabs_sp1$Abundance))^2}, error = function(e) RSQ <<- NA)
    
    Species = levels(melted_filt$Species)[i]
    
    tryCatch({coef = data.frame(Species,
                                estimate_outside,
                                estimate_full, 
                                estimate_part, 
                                IRR_full,
                                IRR_part,
                                RSQ)}, error = function(e) coef <<- NULL)
    
    # cat(i, "\n")
    return(coef)
    
  }) 
  
  df_abun = data.table::rbindlist(all_sp_glm_abund)
  save(df_abun, file = here::here("outputs", "df_abun.RData"))
  return(df_abun)
}