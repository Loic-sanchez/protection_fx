plot_cor_suppl1 = function(){
  
  theme_set(theme_bw())
  
  p1 = ggplot(occu_traits, aes(x = RR_part, y = RR_full, color = AUC)) + 
    geom_point(show.legend = T) +
    # geom_hline(yintercept = median(occu_traits$RR_full), 
    #            linetype = "dashed", 
    #            alpha = 0.5) +
    # geom_vline(xintercept = median(occu_traits$RR_part), 
    #            linetype = "dashed", 
    #            alpha = 0.5) +
    xlab("Effect size - partial protection") +
    ylab("Effect size - full protection") +
    scale_color_viridis_c(name = "AUC", direction = -1, option = "B") + 
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold"),
          panel.grid = element_blank()) 
  
  # p1_marg = ggExtra::ggMarginal(p1, type = "histogram", 
  #                               fill = "#A92E5EFF",
  #                               xparams = list(bins = 50),
  #                               yparams =list(bins = 50))
  
  p2 = ggplot(abun_traits, aes(x = IRR_part, y = IRR_full, colour = RSQ)) + 
    geom_point(show.legend = T) +
    # geom_hline(yintercept = median(abun_traits$IRR_full), 
    #            linetype = "dashed", 
    #            alpha = 0.5) +
    # geom_vline(xintercept = median(abun_traits$IRR_part), 
    #            linetype = "dashed", 
    #            alpha = 0.5) +
    xlab("Effect size - partial protection") +
    ylab("Effect size - full protection") +
    scale_color_viridis_c(name = "R-squared", direction = -1, option = "B") +
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold"),
          panel.grid = element_blank()) 
  
  
  # p2_marg = ggExtra::ggMarginal(p2, type = "histogram", 
  #                               fill = "#A92E5EFF", 
  #                               xparams = list(bins = 50), 
  #                               yparams =list(bins = 50))
  # 
  fig1 = cowplot::plot_grid(p1, p2, 
                     ncol = 2, 
                     labels = c("(a)", "(b)"),
                     label_x = 0.01,
                     label_y = 0.975)
  
  ggsave(here::here("figures", "Figure_1.pdf"),
         width = 11.4, 
         height = 9.5,
         unit = "cm",
         fig1)
  
}

plot_cor_suppl2 = function(){
  
  theme_set(theme_bw())
  
  abun_occu = abun_traits %>% 
    left_join(occu_traits %>% select(Species, RR_full, RR_part), by = "Species") %>% 
    drop_na(RR_full)
  
  p1 = ggplot(abun_occu, aes(x = RR_full, y = IRR_full)) + 
    geom_point(show.legend = T) +
    xlab("Effect size - Occurrence") +
    ylab("Effect size - Abundance") +
    # ylim(c(0,10)) +
    theme(panel.grid = element_blank()) 
  
  p2 = ggplot(abun_occu, aes(x = RR_part, y = IRR_part)) + 
    geom_point(show.legend = T) +
    xlab("Effect size - Occurrence") +
    ylab("Effect size - Abundance") +
    # ylim(c(0,10)) +
    theme(panel.grid = element_blank()) 
  
  fig1 = cowplot::plot_grid(p1, p2, 
                            ncol = 2, 
                            labels = c("(a)", "(b)"),
                            label_x = 0.01,
                            label_y = 0.975)
  
  ggsave(here::here("figures", "Figure_1.pdf"),
         width = 11.4, 
         height = 9.5,
         unit = "cm",
         fig1)

}

plot_histograms_panel = function() {
  
  theme_set(theme_bw())
  
  abun_traits_both$RR = abun_traits_both$IRR
  occu_abun = rbind(occu_traits_both %>% select(Species, RR, Protection),
                    abun_traits_both %>% select(Species, RR, Protection))
  
  occu_abun$Protection = plyr::revalue(occu_abun$Protection,
                                       c("RR_full" = "Occurrence - full protection",
                                         "RR_part" = "Occurrence - partial protection",
                                         "IRR_full" = "Abundance - full protection", 
                                         "IRR_part" = "Abundance - partial protection"))
  occu_abun$Protection = factor(occu_abun$Protection, levels = c("Occurrence - full protection", 
                                                                 "Occurrence - partial protection",
                                                                 "Abundance - full protection",
                                                                 "Abundance - partial protection"))
  
  medians <- occu_abun %>% 
    group_by(Protection) %>% 
    summarise(med_RR = median(RR))
  quartiles <- occu_abun %>%
    group_by(Protection) %>% 
    summarise(qart3_RR = quantile(RR, 0.75),
              qart1_RR = quantile(RR, 0.25))
  
  fig1 = ggplot(occu_abun, 
         aes(x = RR, 
             fill = Protection)) + 
    geom_histogram(alpha = 0.2,
                   color = "black",
                   binwidth = 0.1,
                   aes(y=after_stat(density)),
                   show.legend = F) +
    geom_density(aes(color = Protection), alpha = 0.5, show.legend = F, n = 1024) +
    geom_vline(data = medians, 
               aes(xintercept = med_RR, 
                   color = Protection), 
               # color = "black",
               linetype = 1,
               size = 1,
               show.legend = F) +
    geom_vline(data = quartiles, 
               aes(xintercept = qart1_RR, color = Protection), 
               # color = "black",
               linetype = 4,
               size = 1,
               show.legend = F) +
    geom_vline(data = quartiles, 
               aes(xintercept = qart3_RR, color = Protection), 
               # color = "black",
               linetype = 4,
               size = 1,
               show.legend = F) +
    scale_x_continuous(breaks = seq(0,10, 1), 
                       limits = c(0, 10)) +
    facet_wrap(~Protection) + 
    xlab("Effect size of protection") + 
    ylab("Density") +
    # xlim(c(0,10)) + 
    scale_fill_manual(values = c("#842681FF", "#AE347BFF", "#F56B5CFF", "#FEAC76FF")) +
    scale_color_manual(values = c("#842681FF", "#AE347BFF", "#F56B5CFF", "#FEAC76FF")) +
    # scale_color_gradient2(midpoint = 1, low = "#4D117BFF", high = "#762181FF", mid = "grey") +
    # scale_fill_gradient2(midpoint = 1, low = "#4D117BFF", high = "#762181FF", mid = "grey") +
    theme(strip.background = element_rect(colour="black", fill="white", 
                                          linewidth = 1.25),
          strip.text = element_text(color = "black",
                                    face = "bold",
                                    size = 10),
          panel.grid = element_blank())
  
  ggsave(here::here("figures", "Figure_1.pdf"),
         fig1)
  
}

plot_fx = function() { 
  
  theme_set(theme_bw())
  
  fx_occu_both = ggeffects::ggeffect(model = mod_occu_both, 
                                     terms = c("Trophic.Level",
                                               "Rarity", 
                                               "Protection"),
                                     type = "re")
  fx_abun_both = ggeffects::ggeffect(model = mod_abun_both, 
                                     terms = c("Trophic.Level", 
                                               "Rarity", 
                                               "Protection"),
                                     type = "re")
  fx_abun_both$type = "Abundance"
  fx_occu_both$type = "Occurrence"
  
  fx_both = rbind(fx_occu_both, fx_abun_both)
  fx_both = fx_both %>% 
    rename(Protection = "facet")
  fx_both$Protection = fct_collapse(fx_both$Protection, "Full" = c("RR_full", "IRR_full"), "Partial" = c("RR_part", "IRR_part"))
  fx_both$type = factor(fx_both$type, levels = c("Occurrence", "Abundance"))
  
  labels = c("Occurrence" = "(a) Occurrence","Abundance" = "(b) Abundance")

  fig2 = ggplot(fx_both) +
    geom_line(aes(x = x,
                  y = predicted,
                  colour = group,
                  linetype = Protection),
              linewidth = 1.25, 
              alpha = 0.75) +
    scale_colour_manual(values = c("#F78311FF", "#140B35FF")) +
    facet_wrap(~type, labeller = labeller(type = labels)) +
    xlab("Trophic level") + 
    ylab("Effect size of protection") +
    guides(colour = guide_legend(title = "Rarity"),
           linetype = guide_legend(title = "Protection")) +
    theme(strip.background = element_rect(colour="black", fill="white", 
                                          linewidth = 1.25),
          strip.text = element_text(color = "black",
                                    face = "bold",
                                    size = 10),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 1.25),
          legend.title = element_text(face = "bold"))
  
  ggsave(here::here("figures", "Figure_2.pdf"),
         unit = "in",
         fig2)
  
}

# plot_map = function() {
#   
#   world = rnaturalearth::ne_countries(scale = "small", 
#                                       returnclass = "sf", 
#                                       continent = c("africa", 
#                                                     "north america", 
#                                                     "south america", 
#                                                     "oceania", 
#                                                     "europe", 
#                                                     "asia"))
#   
#   
#   coords_sites = covariates %>% select(SurveyID, Effectiveness) %>%
#     left_join(sites_info %>% select(SurveyID, SiteLongitude, SiteLatitude), by = "SurveyID")
#   
#   ggplot(world) + 
#     geom_sf() +
#     geom_point(data = coords_sites, aes(x = SiteLongitude, 
#                                         y = SiteLatitude,
#                                         color = Effectiveness),
#                size = 0.2) +
#     xlab("Longitude") + 
#     ylab("Latitude") +
#     ggspatial::annotation_north_arrow(location = "bl",
#                            pad_x = unit(0.4, "in"), pad_y = unit(0.25, "in"),
#                            style = north_arrow_fancy_orienteering)
# }

# world = rnaturalearth::ne_countries(scale = "small", 
#                                     returnclass = "sf", 
#                                     continent = c("africa", 
#                                                   "north america", 
#                                                   "south america", 
#                                                   "oceania", 
#                                                   "europe", 
#                                                   "asia"))
# 

# 
# coords_sites$Effectiveness = fct_collapse(coords_sites$Effectiveness,
#                                     "Fully Protected" = c("High No take", "High No take multizoned"),
#                                     "Highly Protected" = c("Medium No take", "Medium No take multizoned"),
#                                     "Lightly Protected" = c("High Restricted take", "High Restricted take multizoned",
#                                                             "Low No take", "Low No take multizoned",
#                                                             "Medium Restricted take", "Medium Restricted take multizoned"),
#                                     "Minimally Protected" = c("Low Restricted take","Low Restricted take multizoned"),
#                                     "Unprotected" = c("Medium Fishing", "Low Fishing", "out"))
# # 
# coords_sites$Effectiveness = as.factor(as.character(
#   plyr::revalue(coords_sites$Effectiveness, c("Fully Protected" = "Full Protection",
#                                         "Highly Protected" = "Partial Protection",
#                                         "Lightly Protected" = "Partial Protection",
#                                         "Minimally Protected" = "Partial Protection",
#                                         "Unprotected" = "Unprotected"))))
# 
# sp_coords = sf::st_as_sf(coords_sites, coords = c("SiteLongitude", "SiteLatitude"))
# sp_coords = sf::st_set_crs(sp_coords, 4269)


# theme_set(theme_bw())
# ggplot(world) + 
#   geom_sf(fill= "#252525") +
#   coord_sf(xlim = c(-180, 170),
#            expand = FALSE) +
#   geom_point(data = coords_sites, 
#              aes(x = SiteLongitude, 
#                  y = SiteLatitude, 
#                  color = Effectiveness, 
#                  group = Effectiveness), 
#              position = position_dodge(width = 0.01), 
#              show.legend = F) +
#   xlab("Longitude") + 
#   ylab("Latitude") +
#   ggspatial::annotation_north_arrow(location = "bl",
#                                     pad_x = unit(0.4, "in"), pad_y = unit(0.25, "in"),
#                                     style = ggspatial::north_arrow_fancy_orienteering) +
#   scale_color_manual(values = c("#0F092CFF", "#8D2369FF", "#F67E14FF")) +
#   theme(panel.grid = element_blank(), 
#         # panel.grid.major = element_line(color = gray(.9), 
#         #                                                               linetype = "dashed", 
#         #                                                               size = 0.5), 
#         panel.background = element_rect(fill = "#969696"))
#   
# scales::show_col(RColorBrewer::brewer.pal(10, "Greys"))
