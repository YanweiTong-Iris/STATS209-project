################################################################
# IPTp and child growth
# Script for single-mediator analysis (Txarm as the exposure)
# Last updated: Sep 20, 2023
################################################################

rm(list = ls())
source(paste0(here::here(), "/0-config.R"))

#--------------------------------------------------------
# load data
#--------------------------------------------------------

data_continuous = readRDS(paste0(data_path, "analysis_data_continuous.RDS"))
data_monthly_round = readRDS(paste0(data_path, "analysis_data_monthly_round.RDS"))
data_monthly_ceiling = readRDS(paste0(data_path, "analysis_data_monthly_ceiling.RDS"))
data_zscore_quarterly = readRDS(paste0(data_path, "analysis_data_zscore_quarterly.RDS"))
data_prevalence_quarterly = readRDS(paste0(data_path, "analysis_data_prevalence_quarterly.RDS"))
data_incidence_3month = readRDS(paste0(data_path,"analysis_data_incidence_quarterly.RDS"))
data_incidence_6month = readRDS(paste0(data_path,"analysis_data_incidence_biannual.RDS"))
data_incidence_12month = readRDS(paste0(data_path,"analysis_data_incidence_annual.RDS"))
data_velocity_1month = readRDS(paste0(data_path,"analysis_data_velocity_1month.RDS"))
data_velocity_2month = readRDS(paste0(data_path,"analysis_data_velocity_2month.RDS")) 
data_velocity_3month = readRDS(paste0(data_path,"analysis_data_velocity_3month.RDS")) 


#--------------------------------------------------------
# create outcome, age, and mediator list 
#--------------------------------------------------------

outcome_monthly_round_continuous = c("haz", "whz", "waz")
outcome_monthly_round_binary = c("haz_ms_stunt", "haz_s_stunt",
                          "whz_ms_waste", "whz_s_waste", "waz_underwt")

outcome_monthly_ceiling_continuous = c("wgv1", "wlz_gv1") 
outcome_monthly_ceiling_binary = c("incident_haz_ms_stunt_agemonthcat",
                                   "incident_haz_s_stunt_agemonthcat",
                                   "incident_whz_ms_waste_agemonthcat",
                                   "incident_whz_s_waste_agemonthcat",
                                   "incident_waz_underwt_agemonthcat")

outcome_z_score_quarter = c("haz_quarter", "whz_quarter", "waz_quarter")
outcome_prevalence_quarterly = c("haz_ms_stunt_quarter", "haz_s_stunt_quarter", "whz_ms_waste_quarter", "whz_s_waste_quarter", "waz_underwt_quarter")

outcome_velocity_1mo = c("wgv1", "wlz_gv1")
outcome_velocity_2mo = c("lgv2", "laz_gv2")
outcome_velocity_3mo = c("wgv3", "wlz_gv3", "lgv3", "laz_gv3")

#"SGA" only eligible at birth, will be added separately
outcome_incidence_3mo = c(
  "incident_haz_ms_stunt_agecat_birth",
  "incident_haz_s_stunt_agecat_birth",
  "incident_whz_ms_waste_agecat_birth",
  "incident_whz_s_waste_agecat_birth",
  "incident_waz_underwt_agecat_birth"
)
outcome_incidence_6mo = c(
  "incident_haz_ms_stunt_age_6_12_month",
  "incident_haz_s_stunt_age_6_12_month",
  "incident_whz_ms_waste_age_6_12_month",
  "incident_whz_s_waste_age_6_12_month",
  "incident_waz_underwt_age_6_12_month"
)
outcome_incidence_12mo = c(
  "incident_haz_ms_stunt_age_1_12",
  "incident_haz_s_stunt_age_1_12",
  "incident_whz_ms_waste_age_1_12",
  "incident_whz_s_waste_age_1_12",
  "incident_waz_underwt_age_1_12"
)

age_list_1mo_ceiling = as.numeric(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
#age_list_1mo_round = as.numeric(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
age_list_2mo = factor(c("<2 months", "2- 4 months", "4- 6 months", "6- 8 months", "8- 10 months", "10- 12 months"), levels = c("<2 months", "2- 4 months", "4- 6 months", "6- 8 months", "8- 10 months", "10- 12 months"))
age_list_3mo = factor(c("0-3 months", ">3-6 months", ">6-9 months", ">9-12 months"), levels = c("0-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))
age_list_3mo_birth = factor(c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"), levels= c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))
age_list_6mo_birth = factor(c("birth", "1 day- <6 months", "6-12 months"), levels = c("birth", "1 day- <6 months", "6-12 months"))
age_list_12mo_birth = factor(c("birth", "1 day- 12 months"), levels =c("birth", "1 day- 12 months"))

# can include "anyhb_28", "anyhb_36" as continuous measures if needed in the future 
binary_maternal_mediators = c("LBW", "preterm", "placentalmal", "anemia_28binary", "anemia_36binary")
#maternal continuous: "gestational_weightchange", "SCF_NPX"
#new_anymalaria for monthly analysis

# missing: non-malarial prenatal infection and maternal inflammation 
# Medications: "antibacterial", "betalactam", "fluoroquinolone", "sulfonamide", "antibacterial_other", "antimalarial", "antiparasitic", "antiviral", "anymed"


#--------------------------------------------------------
# single-mediator analysis wrapper function
#--------------------------------------------------------

# Documentation: single_mediator_analysis
# Description: to identify potential mediators
# Intervention: Txarm OR rand_Txarm (for blind testing)
# Mediators of interest: infant malaria (anymalaria; new_anymalaria), 
#                        low birth weight (LBW), preterm birth (preterm), 
#                        placental malaria (placentalmal), non-malarial prenatal infection, 
#                        anti-microbial use (?), maternal anemia (anemia_28binary, anemia_36binary),
#                        gestational weight gain (gestational_weightchange), 
#                        and maternal inflammation (?). 
# Use GLM with a Gaussian family with identity link for continuous outcomes
# and a binomial family with logit link for dichotomous mediators. 
# Args/Options:
# time_unit: "1 month", "2 month", "quarter", "biannual", "annual"
# age_category: to indicate which age column to use (as a string)
#     - data_continuous & data_prevalence: agemonth_round
#     - data_zscore_quarterly & data_prevalence_quarterly: agecat_birth
#     - data_incidence_1month & data_velocity_1month: agemonth_ceiling
#     - data_velocity_2month: age_2_month
#     - data_velocity_3month: agecat
#     - data_incidence_3month: agecat_birth
#     - data_incidence_6month: age_6_12_month
#     - data_incidence_12month: age_1_12
# outcome (detailed names see above):
# - continuous: "haz", "whz", "waz", "wgv1", "wlz_gv1", "wgv3", "wlz_gv3",
#             "lgv2", "laz_gv2"， "lgv3", "laz_gv3"
# - binary: "SGA", haz_ms_stunt", "haz_s_stunt", "whz_ms_waste", "whz_s_waste", 
#         "waz_underweight", "incident_haz_ms_stunt", "incident_haz_s_stunt", 
#         "incident_whz_ms_waste", "incident_whz_s_waste", "incident_waz_underwt"
# outcome_type: "continuous" OR "binary" (as a string)

single_mediator_analysis <-
  function(data,
           time_unit,
           age_category,
           age_group,
           mediator,
           mediator_type,
           outcome,
           outcome_type
           ) {
    
    print(paste0("mediator: ", mediator, ", outcome: ", outcome, ", age: ", age_group))
    data_age <- data %>% filter(data[[age_category]] == age_group) 
    data_age = data_age %>% filter(complete.cases(data_age[[mediator]], data_age[[outcome]]))
    data_age = data_age%>% mutate(!!mediator := as.numeric(as.character(.[[mediator]])))
    
    if (mediator %in% c("gestational_weightchange", "anemia_28binary", "anemia_36binary")){
      mediator_formula = as.formula(paste(mediator, "~", "Txarm"))
      outcome_formula = as.formula(paste(outcome, "~", "Txarm+", mediator))
      outcome_full_formula = as.formula(paste(outcome, "~", "Txarm +", mediator, "+ Txarm:", mediator))
      adjustment_factor = "Unadjusted"
    }
    else if (mediator %in% c("preterm","LBW")){
      mediator_formula = as.formula(paste(mediator, "~", "Txarm + Txarm:gravidity_cat+ gravidity_cat + GAenroll + enrollage + APdichenroll + educdich + wealthcat + sex"))
      outcome_formula = as.formula(paste(outcome, "~", "Txarm + Txarm:gravidity_cat+ gravidity_cat + GAenroll + enrollage + APdichenroll + educdich + wealthcat+ sex +", mediator))
      outcome_full_formula = as.formula(paste(outcome, "~", "Txarm + Txarm:gravidity_cat+ gravidity_cat + GAenroll + enrollage + APdichenroll + educdich + wealthcat+ sex + ", mediator, "+ Txarm:", mediator))
      adjustment_factor = "Full"
    }
    else if(mediator == "new_anymalaria"){
      mediator_formula = as.formula(paste(mediator, "~", "Txarm + gravidity_cat"))
      outcome_formula = as.formula(paste(outcome, "~", "Txarm +  gravidity_cat +", mediator))
      outcome_full_formula = as.formula(paste(outcome, "~", "Txarm + gravidity_cat + ", mediator, "+ Txarm:", mediator))
      adjustment_factor = "Gravidity only"
    }
    else {
      mediator_formula = as.formula(paste(mediator, "~", "Txarm + Txarm:gravidity_cat+ gravidity_cat + GAenroll + enrollage + APdichenroll + educdich + wealthcat"))
      outcome_formula = as.formula(paste(outcome, "~", "Txarm + Txarm:gravidity_cat+ gravidity_cat + GAenroll + enrollage + APdichenroll + educdich + wealthcat+", mediator))
      outcome_full_formula = as.formula(paste(outcome, "~", "Txarm + Txarm:gravidity_cat+ gravidity_cat + GAenroll + enrollage + APdichenroll + educdich + wealthcat+", mediator, "+ Txarm:", mediator))
      adjustment_factor = "Full- sex"
    }
    
    # Check if any of the variables in the model formula have enough variation
    # if (mediator_type == "binary" && outcome_type == "binary") {
    #   counts = table(data_age[[mediator]], data_age[[outcome]])
    #   print(counts)
    #   if (!all(counts >= 15)) {
    #     message("One of the levels doesn't have enough variability. Skipping analysis.")
    #     return(NULL)
    #   }
    # } else if (mediator_type == "binary") {
    #   counts = table(data_age[[mediator]])
    #   print(counts)
    #   if (!(length(unique(data_age[[mediator]])) >= 2 && all(counts >= 10))) {
    #     message("Variable '", mediator, "' doesn't have enough variability. Skipping analysis.")
    #     return(NULL)
    #   } 
    # } else if (outcome_type == "binary") {
    #   counts = table(data_age[[outcome]])
    #   print(counts)
    #   if (!(length(unique(data_age[[outcome]])) >= 2 && all(counts >= 15))) {
    #     message("Variable '", outcome, "' doesn't have enough variability. Skipping analysis.")
    #     return(NULL)
    #   }
    # }
    
      
    if (mediator_type == "binary") {
      intervention_mediator_fit = glm(
        mediator_formula,
        data = data_age,
        family = binomial(link = "logit"),
        na.action = na.omit
      )
    } else {
      intervention_mediator_fit = glm(mediator_formula,
                                      data = data_age,
                                      family = gaussian(link = "identity"))
    }
    
    if (outcome_type == "binary") {
      outcome_fit = glm(
        outcome_formula,
        data = data_age,
        family = binomial(link = "logit"),
        na.action = na.omit
      )
      outcome_fit_interaction = glm(
        outcome_full_formula,
        data = data_age,
        family = binomial(link = "logit"),
        na.action = na.omit
      )
      
      set.seed(2023)
      mediate.fit = mediate.RR(
        model.m = intervention_mediator_fit,
        model.y = outcome_fit,
        sims = 1000,
        treat = "Txarm",
        control.value = "SP",
        treat.value = "DP",
        mediator = mediator,
        conf.level = 0.95,
        boot = FALSE,
        robustSE = TRUE
      )
      summary_mediate = summary(mediate.fit)
      print(summary_mediate)
      
      set.seed(2023)
      mediate.fit.full = mediate.RR(
        model.m = intervention_mediator_fit,
        model.y = outcome_fit_interaction,
        sims = 1000,
        treat = "Txarm",
        control.value = "SP",
        treat.value = "DP",
        mediator = mediator,
        conf.level = 0.95,
        boot = FALSE,
        robustSE = TRUE
      )
      
      summary_mediate_full = summary(mediate.fit.full)
      print(summary_mediate_full)
      
    } else{
      outcome_fit = glm(outcome_formula,
                        data = data_age,
                        family = gaussian(link = "identity"),
                        na.action = na.omit)
      outcome_fit_interaction = glm(outcome_full_formula,
                                    data = data_age,
                                    family = gaussian(link = "identity"),
                                    na.action = na.omit)
      
      set.seed(2023)
      mediate.fit = mediation::mediate(
        model.m = intervention_mediator_fit,
        model.y = outcome_fit,
        sims = 1000,
        treat = "Txarm",
        control.value = "SP",
        treat.value = "DP",
        mediator = mediator,
        conf.level = 0.95,
        boot = FALSE,
        robustSE = TRUE
      )
      summary_mediate = summary(mediate.fit)
      print(summary_mediate)
      
      set.seed(2023)
      mediate.fit.full = mediation::mediate(
        model.m = intervention_mediator_fit,
        model.y = outcome_fit_interaction,
        sims = 1000,
        treat = "Txarm",
        control.value = "SP",
        treat.value = "DP",
        mediator = mediator,
        conf.level = 0.95,
        boot = FALSE,
        robustSE = TRUE
      )
      
      summary_mediate_full = summary(mediate.fit.full)
      print(summary_mediate_full)
    }
    
    N <- data_age %>% distinct(id) %>% count()
    
    outcome_remark = case_when(
      outcome == "haz_quarter" | outcome == "haz" ~ "height-for-age z score",
      outcome == "whz_quarter" | outcome == "whz"~ "weight-for-height z score",
      outcome == "waz_quarter" | outcome == "waz" ~ "weight-for-age z score",
      outcome == "wgv1" ~ "absolute weight growth velocity (1 month increment)",
      outcome == "wgv3" ~ "absolute weight growth velocity (3 month increment)",
      outcome == "wlz_gv1" ~ "weight velocity based on weight-for-length Z-score (1 month increment)",
      outcome == "wlz_gv3" ~ "weight velocity based on weight-for-length Z-score (3 month increment)",
      outcome == "lgv2" ~ "absolute length growth velocity (2 month increment)",
      outcome == "lgv3" ~ "absolute length growth velocity (3 month increment)",
      outcome == "laz_gv2" ~ "length velocity based on length-for-age Z-score (2 month increment)",
      outcome == "laz_gv3" ~ "length velocity based on length-for-age Z-score (3 month increment)",
      grepl("SGA", outcome) ~ "small for gestational weight",
      outcome == "haz_ms_stunt_quarter" | outcome == "haz_ms_stunt" ~ "prevalence: moderate to severe stunting",
      outcome == "haz_s_stunt_quarter"| outcome == "haz_s_stunt" ~ "prevalence: severe stunting",
      outcome == "whz_ms_waste_quarter" | outcome == "whz_ms_waste" ~ "prevalence: moderate to severe wasting",
      outcome == "whz_s_waste_quarter" | outcome == "whz_s_waste" ~ "prevalence: severe wasting",
      outcome == "waz_underwt_quarter" | outcome == "waz_underwt" ~ "prevalence: underweight",
      grepl("incident_haz_ms_stunt", outcome) ~ "incidence: moderate to severe stunting",
      grepl("incident_haz_s_stunt", outcome) ~ "incidence: severe stunting",
      grepl("incident_whz_ms_waste", outcome) ~ "incidence: moderate to severe wasting",
      grepl("incident_whz_s_waste", outcome) ~ "incidence: severe wasting",
      grepl("incident_waz_underwt", outcome) ~ "incidence: underweight"
    )
    
    mediator_remark = case_when(
      mediator == "anemia_28binary" ~ "anemia at gestational week 28",
      mediator == "anemia_36binary" ~ "anemia at gestational week 36",
      mediator == "placentalmal" ~ "placental malaria",
      mediator == "incidentmalaria"  ~ "infant malaria incidence",
      mediator == "new_anymalaria"  ~ "infant malaria occurance",
      mediator == "LBW"  ~ "low birth weight",
      mediator == "preterm"  ~ "pre-term birth",
      mediator == "gestational_weightchange" ~ "gestational weight change between wk 20 and 36",
      mediator == "placentalLAMPdich" ~ "placental malaria (blood)",
      mediator == "SCF_NPX" ~ "Olink NPX value of Assay SCF"
    )
    
    #d= ACME, z= ADE, tau= total effect, n = prop.mediated
    output_full <- data.frame(
      mediator = mediator,
      mediator_remark = mediator_remark,
      outcome = outcome,
      outcome_remark = outcome_remark,
      time_unit = time_unit,
      age_group = age_group,
      N_from_analysis = N,
      interaction = 1,
      adjustment = adjustment_factor,
      #prop_mediated_avg = round(summary_mediate_full$n.avg, 4),
      #prop_mediated_avg_CI = paste0("[", round(summary_mediate_full$n.avg.ci[1], 4), ", ", round(summary_mediate_full$n.avg.ci[2],4), "]"),
      #prop_mediated_avg_lower_CI = round(summary_mediate_full$n.avg.ci[1], 4),
      #prop_mediated_avg_upper_CI = round(summary_mediate_full$n.avg.ci[2], 4), 
      
      total_effect = round(summary_mediate_full$tau.coef, 4),
      total_effect_CI = paste0("[", round(summary_mediate_full$tau.ci[1],4), ", ", round(summary_mediate_full$tau.ci[2],4), "]"),
      total_effect_lower_CI = round(summary_mediate_full$tau.ci[1],4),
      total_effect_upper_CI = round(summary_mediate_full$tau.ci[2],4),
      
      ACME_average = summary_mediate_full$d.avg,
      ACME_average_CI = paste0("[", round(summary_mediate_full$d.avg.ci[1],4), ", ", round(summary_mediate_full$d.avg.ci[2],4), "]"),
      ACME_average_lower_CI = round(summary_mediate_full$d.avg.ci[1],4),
      ACME_average_upper_CI = round(summary_mediate_full$d.avg.ci[2],4),
      
      ADE_average = summary_mediate_full$z.avg,
      ADE_average_CI = paste0("[", round(summary_mediate_full$z.avg.ci[1],4), ", ", round(summary_mediate_full$z.avg.ci[2],4), "]"),
      ADE_average_lower_CI = round(summary_mediate_full$z.avg.ci[1],4),
      ADE_average_upper_CI = round(summary_mediate_full$z.avg.ci[2],4),
      
      ACME_control = summary_mediate_full$d0,
      ACME_control_CI = paste0("[", round(summary_mediate_full$d0.ci[1],4), ", ", round(summary_mediate_full$d0.ci[2],4), "]"),
      ACME_control_lower_CI = round(summary_mediate_full$d0.ci[1],4),
      ACME_control_upper_CI = round(summary_mediate_full$d0.ci[2],4),
      
      ACME_treated = summary_mediate_full$d1,
      ACME_treated_CI = paste0("[", round(summary_mediate_full$d1.ci[1],4), ", ", round(summary_mediate_full$d1.ci[2],4), "]"),
      ACME_treated_lower_CI = round(summary_mediate_full$d1.ci[1],4), 
      ACME_treated_upper_CI = round(summary_mediate_full$d1.ci[2],4),
      
      ADE_control = summary_mediate_full$z0,
      ADE_control_CI = paste0("[", round(summary_mediate_full$z0.ci[1],4), ", ", round(summary_mediate_full$z0.ci[2],4), "]"),
      ADE_control_lower_CI = round(summary_mediate_full$z0.ci[1],4), 
      ADE_control_upper_CI = round(summary_mediate_full$z0.ci[2],4),
      
      ADE_treated = summary_mediate_full$z1,
      ADE_treated_CI = paste0("[", round(summary_mediate_full$z1.ci[1],4), ", ", round(summary_mediate_full$z1.ci[2],4), "]"),
      ADE_treated_lower_CI = round(summary_mediate_full$z1.ci[1],4), 
      ADE_treated_upper_CI = round(summary_mediate_full$z1.ci[2],4)
      
      #prop_control_avg = summary_mediate_full$n0,
      #prop_control_avg_CI = paste0("[", round(summary_mediate_full$n0.ci[1],4), ", ", round(summary_mediate_full$n0.ci[2],4), "]"),
      #prop_control_avg_lower_CI = round(summary_mediate_full$n0.ci[1],4), 
      #prop_control_avg_upper_CI = round(summary_mediate_full$n0.ci[2],4),
      
      #prop_treated_avg = summary_mediate_full$n1,
      #prop_treated_avg_CI = paste0("[", round(summary_mediate_full$n1.ci[1],4), ", ", round(summary_mediate_full$n1.ci[2],4), "]"),
      #prop_treated_avg_lower_CI =round(summary_mediate_full$n1.ci[1],4), 
      #prop_treated_avg_upper_CI = round(summary_mediate_full$n1.ci[2],4)
    )
    
    output_simple <- data.frame(
      mediator = mediator,
      mediator_remark = mediator_remark,
      outcome = outcome,
      outcome_remark = outcome_remark,
      time_unit = time_unit,
      age_group = age_group,
      N_from_analysis = N,
      interaction = 0,
      adjustment = adjustment_factor,
      #prop_mediated_avg = summary_mediate$n.avg,
      #prop_mediated_avg_CI = paste0("[", round(summary_mediate$n.avg.ci[1],4), ", ", round(summary_mediate$n.avg.ci[2],4), "]"),
      #prop_mediated_avg_lower_CI = round(summary_mediate$n.avg.ci[1],4),
      #prop_mediated_avg_upper_CI = round(summary_mediate$n.avg.ci[2],4),
      
      total_effect = summary_mediate$tau.coef,
      total_effect_CI = paste0("[", round(summary_mediate$tau.ci[1],4), ", ", round(summary_mediate$tau.ci[2],4), "]"),
      total_effect_lower_CI = round(summary_mediate$tau.ci[1],4), 
      total_effect_upper_CI = round(summary_mediate$tau.ci[2],4), 
      
      ACME_average = summary_mediate$d.avg,
      ACME_average_CI = paste0("[", round(summary_mediate$d.avg.ci[1],4), ", ", round(summary_mediate$d.avg.ci[2],4), "]"),
      ACME_average_lower_CI = round(summary_mediate$d.avg.ci[1],4), 
      ACME_average_upper_CI = round(summary_mediate$d.avg.ci[2],4),
      
      ADE_average = summary_mediate$z.avg,
      ADE_average_CI = paste0("[", round(summary_mediate$z.avg.ci[1],4), ", ", round(summary_mediate$z.avg.ci[2],4), "]"),
      ADE_average_lower_CI = round(summary_mediate$z.avg.ci[1],4),
      ADE_average_upper_CI = round(summary_mediate$z.avg.ci[2],4)
      
      )
    
    output = dplyr::bind_rows(output_full, output_simple)
    rownames(output) = NULL
    
    return(output %>% mutate_if(is.numeric, round, digits=4))
  }



#--------------------------------------------------------
# function application
#--------------------------------------------------------

single_mediator_analysis_application <-
  function(data_set,
           crossing_set) {
    full_results_list = list()
    
    for (i in 1:nrow(crossing_set)) {
      try(full_results_list[[i]] <- single_mediator_analysis(
        data = data_set,
        time_unit = as.character(crossing_set[i, "time_unit"]),
        age_category = as.character(crossing_set[i, "age_category"]),
        age_group = as.character(crossing_set[i, "age_group"]),
        mediator = as.character(crossing_set[i, "mediator"]),
        mediator_type = as.character(crossing_set[i, "mediator_type"]),
        outcome = as.character(crossing_set[i, "outcome"]),
        outcome_type = as.character(crossing_set[i, "outcome_type"])))
    }
    
    combined_results_df <- do.call(rbind, full_results_list)
    View(combined_results_df)
    return(combined_results_df)
  }


#--------------------------------------------------------
# create crossing sets 
#--------------------------------------------------------

crossing_malaria_outcome_round = rbind(
  crossing(
    mediator = "new_anymalaria",
    outcome = outcome_monthly_round_continuous,
    time_unit = "1 month",
    age_group = c(1,2,3,4,5,6,7,8,9,10,11,12),
    age_category = "agemonth_round",
    mediator_type = "binary",
    outcome_type = "continuous"
  ), crossing(
    mediator = "new_anymalaria",
    outcome = outcome_monthly_round_binary,
    time_unit = "1 month",
    age_group = c(1,2,3,4,5,6,7,8,9,10,11,12),
    age_category = "agemonth_round",
    mediator_type = "binary",
    outcome_type = "binary"
  ))

crossing_malaria_outcome_ceiling = rbind(crossing(
  mediator = "new_anymalaria",
  outcome = outcome_monthly_ceiling_continuous,
  time_unit = "1 month",
  age_group = factor(c(1,2,3,4,5,6,7,8,9,10,11,12)),
  age_category = "agemonth_ceiling",
  mediator_type = "binary",
  outcome_type = "continuous"
),
crossing(
  mediator = "new_anymalaria",
  outcome = outcome_monthly_ceiling_binary,
  time_unit = "1 month",
  age_group = factor(c(1,2,3,4,5,6,7,8,9,10,11,12)),
  age_category = "agemonth_ceiling",
  mediator_type = "binary",
  outcome_type = "binary"
))

crossing_zscore_3mo =rbind(crossing(
  mediator = binary_maternal_mediators,
  outcome = outcome_z_score_quarter,
  time_unit = "3 month",
  age_group = age_list_3mo_birth,
  age_category = "agecat_birth",
  mediator_type = "binary",
  outcome_type = "continuous"
),crossing(
  mediator = c("gestational_weightchange", "SCF"),
  outcome = outcome_z_score_quarter,
  time_unit = "3 month",
  age_group = age_list_3mo_birth,
  age_category = "agecat_birth",
  mediator_type = "continuous",
  outcome_type = "continuous"
)) %>% mutate(age_group = as.character(age_group))

crossing_prevalence_3mo = rbind(crossing(
  mediator = binary_maternal_mediators,
  outcome = outcome_prevalence_quarterly,
  time_unit = "3 month",
  age_group = age_list_3mo_birth,
  age_category = "agecat_birth",
  mediator_type = "binary",
  outcome_type = "binary"
),crossing(
  mediator = c("gestational_weightchange", "SCF"),
  outcome = outcome_prevalence_quarterly,
  time_unit = "3 month",
  age_group = age_list_3mo_birth,
  age_category = "agecat_birth",
  mediator_type = "continuous",
  outcome_type = "binary"
)) %>% mutate(age_group = as.character(age_group))

crossing_velocity_1mo = rbind(crossing(
  mediator = binary_maternal_mediators,
  outcome = outcome_velocity_1mo,
  time_unit = "1 month",
  age_group = age_list_1mo_ceiling,
  age_category = "agemonth_ceiling",
  mediator_type = "binary",
  outcome_type = "continuous"
),crossing(
  mediator = c("gestational_weightchange", "SCF"),
  outcome = outcome_velocity_1mo,
  time_unit = "1 month",
  age_group = age_list_1mo_ceiling,
  age_category = "agemonth_ceiling",
  mediator_type = "continuous",
  outcome_type = "continuous"
)) %>% mutate(age_group = as.character(age_group))

crossing_velocity_2mo = rbind(crossing(
  mediator = binary_maternal_mediators,
  outcome = outcome_velocity_2mo,
  time_unit = "2 month",
  age_group = age_list_2mo,
  age_category = "age_2_month",
  mediator_type = "binary",
  outcome_type = "continuous"
), crossing(
  mediator = c("gestational_weightchange", "SCF"),
  outcome = outcome_velocity_2mo,
  time_unit = "2 month",
  age_group = age_list_2mo,
  age_category = "age_2_month",
  mediator_type = "continuous",
  outcome_type = "continuous"
)) %>% mutate(age_group = as.character(age_group))

crossing_velocity_3mo = rbind(crossing(
  mediator = binary_maternal_mediators,
  outcome = outcome_velocity_3mo,
  time_unit = "3 month",
  age_group = age_list_3mo,
  age_category = "agecat",
  mediator_type = "binary",
  outcome_type = "continuous"
), crossing(
  mediator = c("gestational_weightchange", "SCF"),
  outcome = outcome_velocity_3mo,
  time_unit = "3 month",
  age_group = age_list_3mo,
  age_category = "agecat",
  mediator_type = "continuous",
  outcome_type = "continuous"
)) %>% mutate(age_group = as.character(age_group))

crossing_incidence_3mo = rbind(
  crossing(
    mediator = binary_maternal_mediators,
    outcome = outcome_incidence_3mo,
    time_unit = "3 month",
    age_group = age_list_3mo_birth,
    age_category = "agecat_birth",
    mediator_type = "binary",
    outcome_type = "binary"
  ),
  crossing(
    mediator = c("gestational_weightchange", "SCF"),
    outcome = outcome_incidence_3mo,
    time_unit = "3 month",
    age_group = age_list_3mo_birth,
    age_category = "agecat_birth",
    mediator_type = "continuous",
    outcome_type = "binary"
  )
) %>% mutate(age_group = as.character(age_group))

crossing_incidence_6mo = rbind(
  crossing(
    mediator = binary_maternal_mediators,
    outcome = outcome_incidence_6mo,
    time_unit = "6 month",
    age_group = age_list_6mo_birth,
    age_category = "age_6_12_month",
    mediator_type = "binary",
    outcome_type = "binary"
  ),
  crossing(
    mediator = c("gestational_weightchange", "SCF"),
    outcome = outcome_incidence_6mo,
    time_unit = "6 month",
    age_group = age_list_6mo_birth,
    age_category = "age_6_12_month",
    mediator_type = "continuous",
    outcome_type = "binary"
  )
) %>% mutate(age_group = as.character(age_group))

crossing_incidence_12mo = rbind(
  crossing(
    mediator = binary_maternal_mediators,
    outcome = outcome_incidence_12mo,
    time_unit = "12 month",
    age_group = age_list_12mo_birth,
    age_category = "age_1_12",
    mediator_type = "binary",
    outcome_type = "binary"
  ),
  crossing(
    mediator = binary_maternal_mediators,
    outcome = "SGA_annual",
    time_unit = "12 month",
    age_group = "birth",
    age_category = "age_1_12",
    mediator_type = "binary",
    outcome_type = "binary"
  ),
  crossing(
    mediator = "gestational_weightchange",
    outcome = outcome_incidence_12mo,
    time_unit = "12 month",
    age_group = age_list_12mo_birth,
    age_category = "age_1_12",
    mediator_type = "continuous",
    outcome_type = "binary"
  ),
  crossing(
    mediator = c("gestational_weightchange", "SCF"),
    outcome = "SGA_annual",
    time_unit = "12 month",
    age_group = "birth",
    age_category = "age_1_12",
    mediator_type = "continuous",
    outcome_type = "binary"
  )
) %>%
  mutate(age_group = as.character(age_group))


# crossing_incidence_3mo = crossing(
#     mediator = binary_maternal_mediators,
#     outcome = outcome_incidence_3mo,
#     time_unit = "3 month",
#     age_group = age_list_3mo_birth,
#     age_category = "agecat_birth",
#     mediator_type = "binary",
#     outcome_type = "binary"
#   ) %>% mutate(age_group = as.character(age_group))
# 
# crossing_incidence_6mo = crossing(
#     mediator = binary_maternal_mediators,
#     outcome = outcome_incidence_6mo,
#     time_unit = "6 month",
#     age_group = age_list_6mo_birth,
#     age_category = "age_6_12_month",
#     mediator_type = "binary",
#     outcome_type = "binary"
#   ) %>% mutate(age_group = as.character(age_group))
# 
# crossing_incidence_12mo = rbind(crossing(
#     mediator = binary_maternal_mediators,
#     outcome = outcome_incidence_12mo,
#     time_unit = "12 month",
#     age_group = age_list_12mo_birth,
#     age_category = "age_1_12",
#     mediator_type = "binary",
#     outcome_type = "binary"
#   ),
#   crossing(
#     mediator = binary_maternal_mediators,
#     outcome = "SGA_annual",
#     time_unit = "12 month",
#     age_group = "birth",
#     age_category = "age_1_12",
#     mediator_type = "binary",
#     outcome_type = "binary"
#    )) %>% mutate(age_group = as.character(age_group))




#--------------------------------------------------------
# apply the single-mediator analysis wrapper function
#--------------------------------------------------------

single_mediator_zscore_3mo = single_mediator_analysis_application(data_set = data_zscore_quarterly, crossing_set = crossing_zscore_3mo)
View(single_mediator_zscore_3mo)
saveRDS(single_mediator_zscore_3mo, paste0(results_path,"aim2_single_mediator_zscore_results_3mo.RDS"))

single_mediator_velocity_1mo = single_mediator_analysis_application(data_set = data_velocity_1month, crossing_set = crossing_velocity_1mo)
View(single_mediator_velocity_1mo)
saveRDS(single_mediator_velocity_1mo, paste0(results_path,"aim2_single_mediator_velocity_results_1mo.RDS"))

single_mediator_velocity_2mo = single_mediator_analysis_application(data_set = data_velocity_2month, crossing_set = crossing_velocity_2mo)
View(single_mediator_velocity_2mo)
saveRDS(single_mediator_velocity_2mo, paste0(results_path,"aim2_single_mediator_velocity_results_2mo.RDS"))

single_mediator_velocity_3mo = single_mediator_analysis_application(data_set = data_velocity_3month, crossing_set = crossing_velocity_3mo)
View(single_mediator_velocity_3mo)
saveRDS(single_mediator_velocity_3mo, paste0(results_path,"aim2_single_mediator_velocity_results_3mo.RDS"))

single_mediator_prevalence_3mo = single_mediator_analysis_application(
  data_set = data_prevalence_quarterly,
  crossing_set = crossing_prevalence_3mo)
View(single_mediator_prevalence_3mo)
saveRDS(single_mediator_prevalence_3mo, paste0(results_path,"aim2_single_mediator_prevalence_results_3mo.RDS"))

single_mediator_incidence_3mo = single_mediator_analysis_application(data_set = data_incidence_3month, crossing_set = crossing_incidence_3mo)
View(single_mediator_incidence_3mo)
saveRDS(single_mediator_incidence_3mo, paste0(results_path,"aim2_single_mediator_incidence_results_3mo.RDS"))

single_mediator_incidence_6mo = single_mediator_analysis_application(data_set = data_incidence_6month, crossing_set = crossing_incidence_6mo)
View(single_mediator_incidence_6mo)
saveRDS(single_mediator_incidence_6mo, paste0(results_path,"aim2_single_mediator_incidence_results_6mo.RDS"))

single_mediator_incidence_12mo = single_mediator_analysis_application(data_set = data_incidence_12month, crossing_set = crossing_incidence_12mo)
View(single_mediator_incidence_12mo)
saveRDS(single_mediator_incidence_12mo, paste0(results_path,"aim2_single_mediator_incidence_results_12mo.RDS"))

single_mediator_malaria_outcome_1mo_ceiling = single_mediator_analysis_application(data_set = data_monthly_ceiling, crossing_set = crossing_malaria_outcome_ceiling)
single_mediator_malaria_outcome_1mo_round = single_mediator_analysis_application(data_set = data_monthly_round, crossing_set = crossing_malaria_outcome_round)
single_mediator_malaria_outcome_1mo = rbind(single_mediator_malaria_outcome_1mo_ceiling, single_mediator_malaria_outcome_1mo_round)
View(single_mediator_malaria_outcome_1mo)
saveRDS(single_mediator_malaria_outcome_1mo, paste0(results_path,"aim2_single_mediator_malaria_outcome_1mo.RDS"))



#--------------------------------------------------------
# !CHECKING!
#--------------------------------------------------------
# prevalence vs. incidence
#View(data_prevalence_quarterly)
for (age in unique(data_prevalence_quarterly$agecat_birth)){
  print(age)
  data_age = data_prevalence_quarterly[data_prevalence_quarterly$agecat_birth == age,]
  print(table(data_age$incident_waz_underwt_agecat_birth))
}
data_age = data_prevalence_quarterly[data_prevalence_quarterly$agecat_birth == "1 day-3 months",]
View(data_age)

#View(data_incidence_3month)
for (age in unique(data_incidence_3month$agecat_birth)){
  print(age)
  data_age = data_incidence_3month[data_incidence_3month$agecat_birth == age,]
  print(table(data_age$incident_waz_underwt_agecat_birth))
}
data_age2 = data_incidence_3month[data_incidence_3month$agecat_birth == "1 day-3 months",]
View(data_age2)


table(data_prevalence_quarterly[data_prevalence_quarterly$agecat_birth == "1 day-3 months",]$waz_underwt_quarter)
View(data_prevalence_quarterly[data_prevalence_quarterly$agecat_birth == "1 day-3 months",])
length(unique(data_prevalence_quarterly[data_prevalence_quarterly$agecat_birth == "1 day-3 months",]$"waz_underwt_quarter"))
single_mediator_analysis(data = data_prevalence_quarterly, time_unit = "3 month",
                         age_category = "agecat_birth",
                         age_group= "1 day-3 months",
                         mediator ="anemia_36binary",
                         mediator_type = "binary",
                         outcome = "waz_underwt_quarter",
                         outcome_type = "binary")

data_age = data_prevalence_quarterly[data_prevalence_quarterly$agecat_birth == ">3-6 months",]
table(data_age$anemia_36binary, data_age$whz_ms_waste_quarter)
glm.fit = glm(
  whz_ms_waste_quarter~ Txarm + anemia_36binary+ Txarm*anemia_36binary,
  data = data_age,
  family = binomial(link = "logit"),
  na.action = na.omit
)
vcov(glm.fit)




#--------------------------------------------------------
# *mediate package test
#--------------------------------------------------------
data(jobs)

b <- lm(job_seek ~ treat + econ_hard + sex + age, data=jobs)
c <- lm(depress2 ~ treat + job_seek + econ_hard + sex + age, data=jobs)
contcont <- mediation::mediate(b, c, sims=50, treat="treat", mediator="job_seek")

summary(contcont)
plot(contcont)

direct.fit <- lm(depress2 ~ treat, data=jobs)
summary(direct.fit)
mediator.fit <- lm(job_seek ~ treat, data=jobs)
summary(mediator.fit)
full.fit <- lm(depress2 ~ treat + job_seek, data=jobs)
summary(full.fit)
full.fit2 <- lm(depress2 ~ treat + job_seek + treat*job_seek, data=jobs)
summary(full.fit2)

fit <- mediate(mediator.fit, full.fit, sim = 1000, treat="treat", mediator="job_seek", robustSE= TRUE)
summary(fit)
fit2 <- mediate(mediator.fit, full.fit2, sim = 1000, treat="treat", mediator="job_seek", robustSE= TRUE)
summary = summary(fit2)
class(summary)

b <- glm(job_seek ~ educ + sex, data=jobs)
c <- glm(depress2 ~ educ + job_seek + job_seek*educ + sex, data=jobs)
d <- glm(depress2 ~ educ + job_seek + sex, data=jobs)
model.cat <- mediate(b, c, treat="educ", mediator="job_seek", sims=50,
                     control.value = "gradwk", treat.value = "somcol")
summary(model.cat)

model.cat <- mediate(b, d, treat="educ", mediator="job_seek", sims=50,
                     control.value = "gradwk", treat.value = "somcol")
summary(model.cat)


