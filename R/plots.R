plot_histograms_margins = function(){
  
  theme_set(theme_bw())
  
  p1 = ggplot(occu_traits, aes(x = RR_part, y = RR_full, color = AUC)) + 
    geom_point(show.legend = T) +
    geom_hline(yintercept = median(occu_traits$RR_full), 
               linetype = "dashed", 
               alpha = 0.5) +
    geom_vline(xintercept = median(occu_traits$RR_part), 
               linetype = "dashed", 
               alpha = 0.5) +
    xlab("Risk Ratio - Partial Protection") +
    ylab("Risk Ratio - Full Protection") +
    scale_color_viridis_c(name = "AUC", direction = -1, option = "B") + 
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold")) 
  
  p1_marg = ggExtra::ggMarginal(p1, type = "histogram", 
                                fill = "#A92E5EFF",
                                xparams = list(bins = 50),
                                yparams =list(bins = 50))
  
  p2 = ggplot(abun_traits, aes(x = IRR_part, y = IRR_full, colour = RSQ)) + 
    geom_point(show.legend = T) +
    geom_hline(yintercept = median(abun_traits$IRR_full), 
               linetype = "dashed", 
               alpha = 0.5) +
    geom_vline(xintercept = median(abun_traits$IRR_part), 
               linetype = "dashed", 
               alpha = 0.5) +
    xlab("Incidence Rate Ratio - Partial Protection") +
    ylab("Incidence Rate Ratio - Full Protection") +
    scale_color_viridis_c(name = "R-squared", direction = -1, option = "B") +
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold")) 
  
  
  p2_marg = ggExtra::ggMarginal(p2, type = "histogram", 
                                fill = "#A92E5EFF", 
                                xparams = list(bins = 50), 
                                yparams =list(bins = 50))
  
  fig1 = cowplot::plot_grid(p1_marg, p2_marg, 
                     ncol = 2, 
                     labels = c("a", "b"),
                     label_x = 0.05,
                     label_y = 0.95)
  
  ggsave(here::here("figures", "Figure_1.pdf"),
         width = 180, 
         height = 150,
         unit = "mm",
         fig1)
  
}

plot_effects = function() { 
  
  fx_occu_both = ggeffects::ggeffect(model = mod_occu_both, 
                                     terms = c("Trophic.Level", "Rarity", "Protection"),
                                     type = "re")
  fx_abun_both = ggeffects::ggeffect(model = mod_abun_both, 
                                     terms = c("Trophic.Level", "Rarity", "Protection"),
                                     type = "re")
  fx_abun_both$type = "Abundance"
  fx_occu_both$type = "Occurrence"
  
  fx_both = rbind(fx_occu_both, fx_abun_both)
  fx_both = fx_both %>% 
    rename(Protection = "facet")
  fx_both$Protection = fct_collapse(fx_both$Protection, "Full" = c("RR_full", "IRR_full"), "Partial" = c("RR_part", "IRR_part"))
  fx_both$type = factor(fx_both$type, levels = c("Occurrence", "Abundance"))
  
  fig2 = ggplot(fx_both) +
    geom_line(aes(x = x,
                  y = predicted,
                  colour = group,
                  linetype = Protection),
              linewidth = 1.25) +
    scale_colour_manual(values = c("#F78311FF", "#140B35FF")) +
    facet_wrap(~type) +
    xlab("Trophic Level") + 
    ylab("Effect Size of Protection") +
    guides(colour = guide_legend(title = "Rarity"),
           linetype = guide_legend(title = "Protection")) +
    theme(strip.background = element_rect(colour="white", fill="white", 
                                          linewidth = 1.25),
          strip.text = element_text(color = "black",
                                    face = "bold.italic",
                                    size = 10),
          panel.border = element_rect(colour = "black", linewidth = 1.25),
          legend.title = element_text(face = "bold"))
  
  ggsave(here::here("figures", "Figure_2.pdf"),
         width = 180, 
         height = 150,
         unit = "mm",
         fig2)
  
}