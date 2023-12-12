################################################################
# IPTp and child growth
# Figures showing single mediator results 
# for z-scores and incidence 
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

# load results -------------------------------------------------------------------
result_single_mediator_zscore_3mo = readRDS(paste0(results_path,"aim2_single_mediator_zscore_results_3mo.RDS")) %>%
  filter(interaction == 0) %>%
  dplyr::select(mediator, outcome, age_group, gravidae,
                ACME_average, ACME_average_lower_CI, ACME_average_upper_CI, 
                ADE_average, ADE_average_lower_CI, ADE_average_upper_CI) %>%
  mutate(age_group = factor(age_group, levels = c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))) %>%
  filter(age_group %in% c("Birth", "1 day-3 months", ">3-6 months")) %>% 
  filter(!mediator %in% c(
      "antibacterial_binary",
      "betalactam_binary",
      "birthweight",
      "birthweight_01kg",
      "TNF"
    )
  )

results_single_mediator_laz = subset(result_single_mediator_zscore_3mo, outcome == "haz_quarter")
results_single_mediator_waz = subset(result_single_mediator_zscore_3mo, outcome == "waz_quarter")
results_single_mediator_wlz = subset(result_single_mediator_zscore_3mo, outcome == "whz_quarter")


# data processing functions -------------------------------------------------------------------

process_data <- function(input_data, outcome_type){
  output_data <- input_data  %>% 
    pivot_longer(
      cols = c("ACME_average", "ACME_average_lower_CI", "ACME_average_upper_CI",
               "ADE_average", "ADE_average_lower_CI", "ADE_average_upper_CI"),
      names_to = c("measure", ".value"),
      names_pattern = "^(ACME|ADE)_(\\w+)"
    ) %>% 
    mutate(mediator_label = case_when(
      mediator == "anemia_28binary" ~ "Maternal anemia",
      mediator == "gestational_weightchange" ~ "Gestational\nweight change (kg)",
      #mediator == "antibacterial_binary" ~ "Antibacterial use",
      #mediator == "betalactam_binary" ~ "Beta−lactam use",
      mediator == "LBW" ~ "Low birth weight",
      mediator == "placentalmal" ~ "Placental malaria",
      mediator == "preterm" ~ "Preterm birth",
      mediator == "SCF" ~ "Stem cell factor (SCF)",
      mediator == "IL4" ~ "Interleukin-4 (IL4)",
      mediator == "IL10" ~ "Interleukin-10 (IL10)",
      mediator == "TRANCE" ~ "TNF-related activation-\ninduced cytokine (TRANCE)",
      #mediator == "TNF" ~ "Tumor necrosis\nfactor (TNF)",
      mediator == "IFN_gamma" ~ "Interferon gamma\n(IFN-gamma)",
      mediator == "birthlength" ~ "Birth length (cm)",
      mediator == "birthweight_kg" ~ "Birth weight (kg)"
    )) %>% 
    mutate(age_group = fct_rev(age_group)) %>% 
    mutate(mediator_label = factor(mediator_label, levels = c("Maternal anemia", 
                                                              "Gestational\nweight change (kg)",
                                                              #"Antibacterial use",
                                                              #"Beta−lactam use",
                                                              "Placental malaria",
                                                              "Stem cell factor",
                                                              "Stem cell factor (SCF)",
                                                              "Interleukin-4 (IL4)",
                                                              "Interleukin-10 (IL10)",
                                                              #"Tumor necrosis\nfactor (TNF)",
                                                              "TNF-related activation-\ninduced cytokine (TRANCE)",
                                                              "Interferon gamma\n(IFN-gamma)",
                                                              "Preterm birth",
                                                              "Low birth weight",
                                                              "Birth length (cm)",
                                                              "Birth weight (kg)"))) %>%
    mutate(gravidity = case_when(gravidae == "all" ~ "All",
                                 gravidae == "multi" ~ "Multigradividy",
                                 gravidae == "single" ~ "Primigravidity"))
  
  if(outcome_type=="binary"){
    output_data <- output_data %>% 
      mutate(favors = case_when(
        average < 1 & average_upper_CI <1 ~ "DP promotes\ngrowth",
        average_lower_CI <=1 & average_upper_CI >= 1 ~ "Null",
        average > 1 & average_lower_CI >1 ~ "SP promotes\ngrowth" 
      )) %>% 
      mutate(favors = factor(favors, levels = c("SP promotes\ngrowth", "Null", "DP promotes\ngrowth"))) 
  }
  
  
  if(outcome_type=="continuous"){
    output_data <- output_data %>% 
      mutate(favors = case_when(
        average < 0 & average_upper_CI < 0 ~ "SP promotes\ngrowth",
        average_lower_CI <= 0 & average_upper_CI >= 0 ~ "Null",
        average > 0 & average_lower_CI > 0 ~ "DP promotes\ngrowth" 
      )) %>% 
      mutate(favors = factor(favors, levels = c("SP promotes\ngrowth", "Null", "DP promotes\ngrowth"))) 
  }
  
  
  return(output_data)
  
}

iptp_color = c("#9e161d", "#909190","#164c9e")



# LAZ -------------------------------------------------------------------
laz_l <- process_data(input_data = results_single_mediator_laz,
                      outcome_type = "continuous")

ACME_plot = ggplot(laz_l %>% filter(measure=="ACME"), aes(x = age_group, y = average)) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(aes(col = favors), position = position_dodge(width=0.4), size=1.5) + 
  geom_linerange(aes(ymin = average_lower_CI, ymax = average_upper_CI,
                     col = favors), position = position_dodge(width=0.2)) + 
  scale_color_manual(values =iptp_color, drop = FALSE) +
  scale_y_continuous(breaks = c(-0.5, -0.25, -0.1,0, 0.1, 0.25, 0.5),
                     labels = c(-0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5)) +
  coord_flip() +
  theme_minimal() + 
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size= 6),
    axis.text.y = element_text(size= 7), 
    legend.position = "none"
  ) +
  facet_grid(mediator_label~gravidity) +
  theme(strip.text.y = element_text(angle = 0))+
  ggtitle("Average causal mediated effect (NIE)") + 
  theme(plot.title= element_text(hjust = -0.1))

# Print the plot
print(ACME_plot)

ADE_plot = ggplot(laz_l %>% filter(measure=="ADE"), aes(x = age_group, y = average)) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(aes(col = favors), position = position_dodge(width=0.4), size= 1.5) + 
  geom_linerange(aes(ymin = average_lower_CI, ymax = average_upper_CI,
                     col = favors), position = position_dodge(width=0.2)) + 
  scale_color_manual(values = iptp_color, drop = FALSE) +
  scale_y_continuous(breaks = c(-0.5, -0.25, -0.1,0, 0.1, 0.25, 0.5),
                     labels = c(-0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5)) +
  coord_flip() +
  theme_minimal() + 
  ylab("Mean difference in Z-score (95% CI)") + 
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size= 7), 
        legend.position = "bottom",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())  +
  facet_grid(mediator_label~gravidity) +
  theme(strip.text.y = element_text(angle = 0))+
  ggtitle("Average causal direct effect (NDE)") + 
  theme(plot.title= element_text(hjust = -0.1))

combined_plot <- grid.arrange(ACME_plot, ADE_plot, nrow=2, heights = c(2, 2.2))

ggsave(combined_plot, filename = paste0(figure_path, "plot-mediation-laz_stratified_full.png"),
       width=9, height=12)



# WLZ -------------------------------------------------------------------
wlz_l <- process_data(input_data = results_single_mediator_wlz,
                      outcome_type = "continuous")

ACME_plot = ggplot(wlz_l %>% filter(measure=="ACME"), aes(x = age_group, y = average)) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(aes(col = favors), position = position_dodge(width=0.2), size= 1.5) + 
  geom_linerange(aes(ymin = average_lower_CI, ymax = average_upper_CI,
                     col = favors), position = position_dodge(width=0.2)) + 
  scale_color_manual(values = iptp_color, drop = FALSE) +
  scale_y_continuous(breaks = c(-0.5, -0.25, -0.1,0, 0.1, 0.25, 0.5),
                     labels = c(-0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5)) +
  coord_flip() +
  theme_minimal() + 
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size= 7), 
    legend.position = "none"
  ) +
  facet_grid(mediator_label~gravidity) +
  theme(strip.text.y = element_text(angle = 0))+
  ggtitle("Average causal mediated effect (NIE)") + 
  theme(plot.title= element_text(hjust = -0.1))


ADE_plot = ggplot(wlz_l %>% filter(measure=="ADE"), aes(x = age_group, y = average)) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(aes(col = favors), position = position_dodge(width=0.2), size= 1.5) + 
  geom_linerange(aes(ymin = average_lower_CI, ymax = average_upper_CI,
                     col = favors), position = position_dodge(width=0.2)) + 
  scale_color_manual(values = iptp_color, drop = FALSE) +
  scale_y_continuous(breaks = c(-0.5, -0.25, -0.1,0, 0.1, 0.25, 0.5),
                     labels = c(-0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5)) +
  coord_flip() +
  theme_minimal() + 
  ylab("Mean difference in Z-score (95% CI)") + 
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size= 7), 
        legend.position = "bottom",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())  +
  facet_grid(mediator_label~gravidity) +
  theme(strip.text.y = element_text(angle = 0))+
  ggtitle("Average causal direct effect (NDE)") + 
  theme(plot.title= element_text(hjust = -0.1))

combined_plot <- grid.arrange(ACME_plot, ADE_plot, nrow=2, heights = c(2, 2.2))



ggsave(combined_plot, filename = paste0(figure_path, "plot-mediation-wlz-stratified_full.png"),
       width=8, height=12)


