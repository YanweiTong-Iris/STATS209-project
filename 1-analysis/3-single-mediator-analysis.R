################################################################
# IPTp and child growth
# Script for single-mediator analysis
# (Txarm as the exposure; sratified by gravidity)
# Last updated: Dec 4, 2023
################################################################

rm(list = ls())
source(paste0(here::here(), "/0-config.R"))

#--------------------------------------------------------
# load data
#--------------------------------------------------------

data_zscore_quarterly = readRDS(paste0(here::here(), "/0-input-data/analysis_data_zscore_quarterly.RDS"))

#--------------------------------------------------------
# create outcome, age, and mediator list 
#--------------------------------------------------------

outcome_z_score_quarter = c("haz_quarter", "whz_quarter", "waz_quarter")

age_list_3mo_birth = factor(c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"), levels= c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))

binary_maternal_mediators = c("LBW", "preterm", "placentalmal", "anemia_28binary")

continuous_mediator = c("gestational_weightchange", "SCF", "IL4", "IL10",
                        "TRANCE", "IFN_gamma","birthweight_kg", "birthlength")

#--------------------------------------------------------
# single-mediator analysis wrapper function
#--------------------------------------------------------

# Documentation: single_mediator_analysis (statified by gravidity)
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
           gravidity_strata,
           mediator,
           mediator_type,
           outcome,
           outcome_type
           ) {
    
    print(paste0("mediator: ", mediator, ", outcome: ", outcome, ", age: ", age_group))
    data_age <- data %>% filter(data[[age_category]] == age_group) 
    data_age = data_age %>% filter(complete.cases(data_age[[mediator]], data_age[[outcome]]))
    data_age = data_age%>% mutate(!!mediator := as.numeric(as.character(.[[mediator]])))
    
    if(gravidity_strata == "all"){
      data_stratified = data_age
    }else {
      data_stratified = data_age %>% filter(gravidity_cat == gravidity_strata)
    }
      
      if (gravidity_strata == "all") {
          mediator_formula = as.formula(paste(mediator, "~", "Txarm + sex+ gravidity_cat + GAenroll + enrollage + APdichenroll + educdich + wealthcat"))
          outcome_formula = as.formula(paste(outcome, "~", "Txarm + sex+ gravidity_cat + GAenroll + enrollage + APdichenroll + educdich + wealthcat +", mediator))
          outcome_full_formula = as.formula(paste(outcome, "~", "Txarm + sex+ gravidity_cat + GAenroll + enrollage + APdichenroll + educdich + wealthcat + ", mediator, "+ Txarm:", mediator))
          adjustment_factor = "Full"
    } else{
          mediator_formula = as.formula(paste(mediator, "~", "Txarm + sex+ GAenroll + enrollage + APdichenroll + educdich + wealthcat"))
          outcome_formula = as.formula(paste(outcome, "~", "Txarm + sex+ GAenroll + enrollage + APdichenroll + educdich + wealthcat +", mediator))
          outcome_full_formula = as.formula(paste(outcome, "~", "Txarm + sex+ GAenroll + enrollage + APdichenroll + educdich + wealthcat + ", mediator, "+ Txarm:", mediator))
          adjustment_factor = "Full"
    }
      

    # Check if any of the variables in the model formula have enough variation
    if (mediator_type == "binary" && outcome_type == "binary") {
      counts = table(data_stratified[[mediator]], data_stratified[[outcome]])
      print(counts)
      if (!all(counts >= 5)) {
        message("One of the levels doesn't have enough variability. Skipping analysis.")
        return(NULL)
      }
    } else if (mediator_type == "binary") {
      counts = table(data_stratified[[mediator]])
      print(counts)
      if (!(length(unique(data_stratified[[mediator]])) >= 2 && all(counts >= 5))) {
        message("Variable '", mediator, "' doesn't have enough variability. Skipping analysis.")
        return(NULL)
      }
    } else if (outcome_type == "binary") {
      counts = table(data_stratified[[outcome]])
      print(counts)
      if (!(length(unique(data_stratified[[outcome]])) >= 2 && all(counts >= 5))) {
        message("Variable '", outcome, "' doesn't have enough variability. Skipping analysis.")
        return(NULL)
      }
    }
    
      
    if (mediator_type == "binary") {
      intervention_mediator_fit = glm(
        mediator_formula,
        data = data_stratified,
        family = binomial(link = "logit"),
        na.action = na.omit
      )
    } else {
      intervention_mediator_fit = glm(mediator_formula,
                                      data = data_stratified,
                                      family = gaussian(link = "identity"))
    }
    
    if (outcome_type == "binary") {
      outcome_fit = glm(
        outcome_formula,
        data = data_stratified,
        family = binomial(link = "logit"),
        na.action = na.omit
      )
      outcome_fit_interaction = glm(
        outcome_full_formula,
        data = data_stratified,
        family = binomial(link = "logit"),
        na.action = na.omit
      )
      
      outcome_fit_RR = glm(
        outcome_formula,
        data = data_stratified,
        family = poisson(link = "log"),
        na.action = na.omit
      )
      outcome_fit_interaction_RR = glm(
        outcome_full_formula,
        data = data_stratified,
        family = poisson(link = "log"),
        na.action = na.omit
      )
      
      # aim1 effect
      estimates <- summary(outcome_fit_RR)$coefficients
      coefs = names(coef(outcome_fit_RR))
      treatment_row = as.numeric(grep("armDP$", coefs))
      aim1_point_estimate <- exp(estimates[treatment_row, 1])
      
      inter_estimates <- summary(outcome_fit_interaction_RR)$coefficients
      inter_coefs = names(coef(outcome_fit_interaction_RR))
      inter_treatment_row = as.numeric(grep("armDP$", coefs))
      aim1_point_estimate_inter <- exp(inter_estimates[inter_treatment_row, 1])
      
      
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
                        data = data_stratified,
                        family = gaussian(link = "identity"),
                        na.action = na.omit)
      outcome_fit_interaction = glm(outcome_full_formula,
                                    data = data_stratified,
                                    family = gaussian(link = "identity"),
                                    na.action = na.omit)
      
      # aim1 effect
      estimates <- summary(outcome_fit)$coefficients
      coefs = names(coef(outcome_fit))
      treatment_row = as.numeric(grep("armDP$", coefs))
      aim1_point_estimate <- estimates[treatment_row, 1]
      
      inter_estimates <- summary(outcome_fit_interaction)$coefficients
      inter_coefs = names(coef(outcome_fit_interaction))
      inter_treatment_row = as.numeric(grep("armDP$", coefs))
      aim1_point_estimate_inter <- inter_estimates[inter_treatment_row, 1]
      
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
    
    N <- data_stratified %>% distinct(id) %>% count()
    
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
      mediator == "SCF" ~ "Olink NPX value of Assay SCF"
    )
    
    #d= ACME, z= ADE, tau= total effect, n = prop.mediated
    output_full <- data.frame(
      mediator = mediator,
      mediator_remark = mediator_remark,
      outcome = outcome,
      outcome_remark = outcome_remark,
      gravidae= gravidity_strata,
      time_unit = time_unit,
      age_group = age_group,
      N_from_analysis = N,
      interaction = 1,
      adjustment = adjustment_factor,
      
      aim1_effect = aim1_point_estimate_inter,
     
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
    )
    
    output_simple <- data.frame(
      mediator = mediator,
      mediator_remark = mediator_remark,
      outcome = outcome,
      outcome_remark = outcome_remark,
      gravidae= gravidity_strata,
      time_unit = time_unit,
      age_group = age_group,
      N_from_analysis = N,
      interaction = 0,
      adjustment = adjustment_factor,
      
      aim1_effect = aim1_point_estimate,
      
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
    # output_list[[i]] = output
    # 
    # }
    # output_df = bind_rows(output_list)
    return(output %>% mutate_if(is.numeric, round, digits=2))
  }
    

#--------------------------------------------------------
# function application
#--------------------------------------------------------

single_mediator_analysis_application <-
  function(data_set,
           crossing_set) {
    full_results_list = list()
    
    for (i in 1:nrow(crossing_set)) {
      try(full_results_list[[i]] <-
            single_mediator_analysis(
                    data = data_set,
                    time_unit = as.character(crossing_set[i, "time_unit"]),
                    age_category = as.character(crossing_set[i, "age_category"]),
                    age_group = as.character(crossing_set[i, "age_group"]),
                    mediator = as.character(crossing_set[i, "mediator"]),
                    mediator_type = as.character(crossing_set[i, "mediator_type"]),
                    outcome = as.character(crossing_set[i, "outcome"]),
                    outcome_type = as.character(crossing_set[i, "outcome_type"]),
                    gravidity_strata = as.character(crossing_set[i, "gravidity_strata"])
                  ))}
  
    combined_results_df <- do.call(rbind, full_results_list)
    View(combined_results_df)
    return(combined_results_df)
  }


#--------------------------------------------------------
# create crossing sets 
#--------------------------------------------------------

crossing_zscore_3mo =rbind(crossing(
  mediator = binary_maternal_mediators,
  outcome = outcome_z_score_quarter,
  time_unit = "3 month",
  age_group = age_list_3mo_birth,
  age_category = "agecat_birth",
  mediator_type = "binary",
  outcome_type = "continuous",
  gravidity_strata = c("all", "single", "multi")
),crossing(
  mediator = continuous_mediator,
  outcome = outcome_z_score_quarter,
  time_unit = "3 month",
  age_group = age_list_3mo_birth,
  age_category = "agecat_birth",
  mediator_type = "continuous",
  outcome_type = "continuous",
  gravidity_strata = c("all", "single", "multi")
)) %>% mutate(age_group = as.character(age_group))



#--------------------------------------------------------
# apply the single-mediator analysis wrapper function
#--------------------------------------------------------

#prevent using scientific notations
options(scipen = 999)

single_mediator_zscore_3mo_stratified = single_mediator_analysis_application(data_set = data_zscore_quarterly, crossing_set = crossing_zscore_3mo)
View(single_mediator_zscore_3mo_stratified)
saveRDS(single_mediator_zscore_3mo_stratified, paste0(results_path,"aim2_single_mediator_zscore_results_3mo.RDS"))



