################################################################
# IPTp and child growth
# Figures showing intervention-mediator and mediator-outcome results
# for z-scores and incidence 
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))


################################################################
#INTERVENTION-MEDIATOR
################################################################

IM_zscore_3mo = readRDS(paste0(here::here(), "/4-result-data/intervention_mediator_main_results_stratified.RDS")) %>%
  filter(gravidity_strata %in% c("multi", "single")) %>%
  mutate(gravidae = ifelse(gravidity_strata=="single", "Primigravidae", "Multigravidae")) %>% 
  mutate(gravidae = factor(gravidae, levels = c("Primigravidae", "Multigravidae"))) %>%
  filter(!mediator %in% c("anemia_36binary", "antiviral_binary", "antiparasitic_binary",
                          "antibacterial_other_binary", "antimalarial_binary", "sulfonamide_binary", 
                          "antibacterial_binary", "betalactam_binary", "birthweight",
                          "birthweight_01kg", "TNF"))%>% 
  mutate(mediator_label = case_when(
    mediator == "gestational_weightchange" ~ "Gestational weight\nchange (kg)",
    mediator == "LBW" ~ "Low birth weight",
    mediator == "placentalmal" ~ "Placental malaria",
    mediator == "preterm" ~ "Preterm birth",
    mediator == "SCF" ~ "Stem cell factor (SCF)",
    mediator == "IL4" ~ "Interleukin-4 (IL4)",
    mediator == "IL10" ~ "Interleukin-10 (IL10)",
    mediator == "TRANCE" ~ "TNF-related activation-\ninduced cytokine\n(TRANCE)",
    mediator == "IFN_gamma" ~ "Interferon gamma\n(IFN-gamma)",
    mediator == "anemia_28binary" ~ "Maternal anemia",
    mediator == "birthlength" ~ "Birth length (cm)",
    mediator == "birthweight_kg" ~ "Birth weight (kg)"
  )) %>% 
  mutate(mediator_label = factor(mediator_label, levels = rev(c( "Maternal anemia",
                                                                                     "Gestational weight\nchange (kg)",
                                                                                     "Placental malaria",
                                                                                     "Stem cell factor (SCF)",
                                                                                     "Interleukin-4 (IL4)",
                                                                                     "Interleukin-10 (IL10)",
                                                                                     "TNF-related activation-\ninduced cytokine\n(TRANCE)",
                                                                                     "Interferon gamma\n(IFN-gamma)",
                                                                                     "Preterm birth",
                                                                                     "Low birth weight",
                                                                                     "Birth length (cm)",
                                                                                     "Birth weight (kg)"
  )))) 
  
binary_mediators_main = c("anemia_28binary", "placentalmal", "preterm", "LBW")

continuous_mediator = c("gestational_weightchange", "IL4", "IL10", "SCF",
                        "TRANCE", "IFN_gamma", "birthweight_kg", "birthlength")

# make plots -------------------------------------------------------------------

IM_plot_binary <- ggplot(IM_zscore_3mo %>% filter(mediator %in% binary_mediators_main),
                         aes(x = mediator_label, y = adjusted_point_estimate, col = gravidae)) +
  geom_point(size = 1.5, position = position_dodge(width = 0.4)) +
  geom_linerange(aes(ymin = adjusted_lower_95CI, ymax = adjusted_upper_95CI, col = gravidae),
                 position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c("#999999", "#6600CC"), drop = FALSE) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position="none",
        axis.title.y = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 7)) +
  ylab("Risk Ratio (95% CI)") +
  ggtitle("Dichotomous Mediators") +
  theme(strip.text.y = element_text(angle = 0),
        plot.title =element_text(hjust = -0.1, size = 11))

IM_plot_continuous <- ggplot(IM_zscore_3mo %>% filter(mediator %in% continuous_mediator),
                         aes(x = mediator_label, y = adjusted_point_estimate, col = gravidae)) +
  geom_point(size = 1.5, position = position_dodge(width = 0.4)) +
  geom_linerange(aes(ymin = adjusted_lower_95CI, ymax = adjusted_upper_95CI, col = gravidae),
                 position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#999999", "#6600CC"), drop = FALSE) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text.y = element_text(size= 9),
        axis.text.x = element_text(size= 7))+
  ylab("Mean difference (95% CI)") +
  ggtitle("Continuous Mediators") +
  theme(strip.text.y = element_text(angle = 0),
        plot.title =element_text(hjust = -0.2, size = 11))

IM_plot <- grid.arrange(
  IM_plot_binary, 
  IM_plot_continuous, 
  heights = c(2, 4.3), 
  ncol = 1,
  top = textGrob("Intervention-Mediator Effects by Gravidity", gp = gpar(fontsize = 13), hjust = 1)
)

IM_plot

ggsave(IM_plot, filename = paste0(figure_path, "plot-tx-med.png"),
       width=7, height=7)


################################################################
#MEDIATOR-OUTCOME
################################################################

MO_zscore_3mo = readRDS(paste0(here::here(), "/4-result-data/mediator_outcome_zscore_results_3mo_stratified.RDS")) %>%
  mutate(age_group = factor(age_group, levels = c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months")))



# define data processing function -------------------------------------------------------------------
process_data <- function(data, outcome_type){
  output_data <- data %>%  filter(!independent_variable %in% c("anemia_36binary",
                                                               #"anemia_28binary",
                                                               "antiviral_binary",
                                                               "antiparasitic_binary",
                                                               "antibacterial_other_binary",
                                                               "antimalarial_binary",
                                                               "sulfonamide_binary",
                                                               "fluoroquinolone_binary",
                                                               "antibacterial_binary",
                                                               "betalactam_binary",
                                                               "birthweight",
                                                               "birthweight_01kg",
                                                               #"LBW",
                                                               "TNF"))%>% 
    
    mutate(independent_variable_label = case_when(
      independent_variable == "gestational_weightchange" ~ "Gestational weight\nchange (kg)",
      independent_variable == "LBW" ~ "Low birth weight",
      independent_variable == "placentalmal" ~ "Placental malaria",
      independent_variable == "preterm" ~ "Preterm birth",
      independent_variable == "SCF" ~ "Stem cell factor (SCF)",
      independent_variable == "IL4" ~ "Interleukin-4 (IL4)",
      independent_variable == "IL10" ~ "Interleukin-10 (IL10)",
      independent_variable == "TRANCE" ~ "TNF-related activation-\ninduced cytokine (TRANCE)",
      #independent_variable == "TNF" ~ "Tumor necrosis factor (TNF)",
      independent_variable == "IFN_gamma" ~ "Interferon gamma\n(IFN-gamma)",
      # independent_variable == "antibacterial_binary" ~ "Antibacterial use",
      # independent_variable == "betalactam_binary" ~ "Beta-lactam use",
      independent_variable == "anemia_28binary" ~ "Maternal anemia",
      independent_variable == "birthlength" ~ "Birth length (cm)",
      independent_variable == "birthweight_kg" ~ "Birth weight (kg)"
    )) %>% 
    mutate(independent_variable_label = factor(independent_variable_label, levels = c( "Maternal anemia",
                                                                                       "Gestational weight\nchange (kg)",
                                                                                       #"Antibacterial use",
                                                                                       # "Beta-lactam\nuse",
                                                                                       "Placental malaria",
                                                                                       "Stem cell factor (SCF)",
                                                                                       "Interleukin-4 (IL4)",
                                                                                       "Interleukin-10 (IL10)",
                                                                                       #"Tumor necrosis factor (TNF)",
                                                                                       "TNF-related activation-\ninduced cytokine (TRANCE)",
                                                                                       "Interferon gamma\n(IFN-gamma)",
                                                                                       "Preterm birth",
                                                                                       "Low birth weight",
                                                                                       "Birth length (cm)",
                                                                                       "Birth weight (kg)"
    ))) %>% 
    mutate(age_group = fct_rev(age_group))  %>% 
    filter(gravidae!="all") %>% 
    mutate(gravidae = ifelse(gravidae=="single", "Primigravidae", "Multigravidae")) %>% 
    mutate(gravidae = factor(gravidae, levels = c("Primigravidae", "Multigravidae")))
  
  return(output_data)
}


# process Z-score data -------------------------------------------------------------------
zscore_plotdata = process_data(MO_zscore_3mo, outcome_type="continuous") %>% 
  mutate(dependent_variable_label = case_when(
    dependent_variable == "haz_quarter" ~ "LAZ",
    dependent_variable == "whz_quarter" ~ "WLZ",
    dependent_variable == "waz_quarter" ~ "WAZ"
  ))



# make plots -------------------------------------------------------------------

Zscore_plot <- ggplot(zscore_plotdata %>% 
                        filter(dependent_variable_label != "WAZ") %>%
                        filter(age_group %in% c("Birth", "1 day-3 months", ">3-6 months")),
       aes(x = age_group, y = point_estimate, col = gravidae)) +
  geom_point(size = 1.2, position = position_dodge(width = 0.6)) +
  geom_linerange(aes(ymin = lower_95CI, ymax = upper_95CI, col = gravidae),
                 position = position_dodge(width = 0.6)) +
  facet_grid(independent_variable_label~dependent_variable_label,
             scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#999999", "#6600CC"), drop = FALSE) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text.y = element_text(size= 7),
        axis.text.x = element_text(size= 6))+
  xlab("Age Group") +
  ylab("Mean difference in Z-score (95% CI)") +
  ggtitle("Mediator-Outcome Effects by Infant Age and Gravidity") +
  theme(strip.text.y = element_text(angle = 0),
  plot.title =element_text(hjust = -0.25, size = 13))

ggsave(Zscore_plot, filename = paste0(figure_path, "plot-med-outcome-zscore.png"),
       width=8, height=10)


