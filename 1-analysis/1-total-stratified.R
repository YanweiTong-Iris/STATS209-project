################################################################
# IPTp and child growth
# Script for Aim1: estimate the impact of IPTp on child growth
# Last updated: Oct 14, 2023
################################################################

rm(list = ls())
source(paste0(here::here(), "/0-config.R"))

#prevent using scientific notations
options(scipen = 999)

#--------------------------------------------------------
# load data
#--------------------------------------------------------

data_zscore_monthly = readRDS(paste0(data_path, "analysis_data_zscore_monthly.RDS"))
data_zscore_quarterly = readRDS(paste0(data_path, "analysis_data_zscore_quarterly.RDS"))

#--------------------------------------------------------
# create outcome, age, and modifier list 
#--------------------------------------------------------

outcome_zscore = c("haz", "whz", "waz")
outcome_zscore_quarter = c("haz_quarter", "whz_quarter", "waz_quarter")

age_list_1mo_round = as.numeric(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
age_list_1mo_ceiling = as.numeric(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
age_list_3mo_birth = factor(c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"), levels= c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))

#--------------------------------------------------------
# aim1 analysis wrapper function
#--------------------------------------------------------

# Documentation: aim1_analysis
# Usage: aim1_analysis(data, time_unit, age_group, outcome, outcome_type, treatment, modifier)
# Description: Estimate the impact of IPTp on child growth by each age category (stratified by gravidity: both, single, multi)
# For continuous outcomes, use GLM with a Gaussian family and identity link.
# For binary outcomes, use GLM with a Poisson family and log link.
# If the sample size is sufficient, investigate whether effects are modified by
# the modifiers below in both a relative and additive scale.
# Return an outcome data frame with outcome name, age group, N from analysis,
# modifier name (or “none” if no modifier), point estimate, SE, lower bound of 95% CI of the point estimate,
# upper bound of 95% CI of the point estimate, additive modifier level (or NA if no modifier),
# lower bound of 95% CI of the additive modifier level, upper bound of 95% CI of the additive modifier level,
# relative modifier level (or NA if no modifier), lower bound of 95% CI of the relative modifier level,
# upper bound of 95% CI of the relative modifier level

# Args/Options:
# time_unit: "1 month", "2 month", "quarter", "biannual", "annual"
# age_group: 
# data_continuous & data_prevalence: agemonth_round
# data_incidence_1month & data_velocity_1month: agemonth_ceiling
# data_velocity_2month: age_2_month
# data_velocity_3month: agecat
# data_incidence_3month: agecat_birth
# data_incidence_6month: age_6_12_month
# data_incidence_12month: age_1_12
# outcome: outcome name, as a string (w/ different time units at the end)
# continuous: "haz", "whz", "waz", "wgv", "wlz_gv", "lgv", "laz_gv"
# binary: "haz_ms_stunt", "haz_s_stunt", "whz_ms_waste", "whz_s_waste", "waz_underweight"
#         "incident_haz_ms_stunt", "incident_haz_s_stunt", "incident_whz_ms_waste", 
#         "incident_whz_s_waste", "incident_waz_underwt","SGA"
# outcome_type: type of outcome, as a string (either "continuous" or "binary")
# # treatment: either "Txarm" or "rand_Txarm")
# modifier_name: name of effect modifier, as a string (default = NULL)
#   "enrollage" (maternal age), "gravidity_cat" (first pregnancy or not), 
#   "gestational_weightchange_cat" (weight difference at gestaional week 36 and 20), 
#   "mother_heightcat" (the mother's height measured during pregnancy),
#   "preterm" (pre-term delivery), "sex" (child sex), "anymalaria" (child malaria infection),
#   and "wet_season" (wet season or not, not available for biannual and annual datasets),

# Returns: the data frame as described above
# Output: prints the data frame described above

aim1_analysis <-
  function(data,
           time_unit,
           age_group,
           age_category,
           outcome,
           outcome_type,
           treatment = "Txarm",
           modifier = "NA") {
    
    print(paste0("outcome: ", outcome, ", age: ", age_group, ", modifier: ", modifier))
    
    data_age = data[data[[age_category]] == age_group, ]
    
    output_list = list()
    strata_list = list(
      "all" = c("single", "multi"),
      "single" = c("single"),
      "multi" = c("multi")
    )
    
    for (i in 1:3) {
      strata_name = names(strata_list)[[i]]
      strata = strata_list[[strata_name]]
      data_stratified = data_age %>% filter(gravidity_cat %in% strata)
      
      modifier_counts <- table(data_stratified[[modifier]])
      
      if (modifier == "NA" && strata_name != "all") {
        formula <- as.formula(paste(outcome, "~", treatment,
                                    "+sex +GAenroll + enrollage + APdichenroll + educdich + wealthcat"))
      } else if(modifier == "NA" && strata_name == "all"){
        formula <- as.formula(paste(outcome, "~", treatment,
                                    "+sex + gravidity_cat + GAenroll + enrollage + APdichenroll + educdich + wealthcat"))
      } else if (strata_name != "all" && length(modifier_counts) > 1 && all(modifier_counts > 3)){
        data_stratified <- data_stratified[!is.na(data_stratified[[modifier]]), ]
        formula <-
          as.formula(paste(outcome, "~", treatment, "+", modifier, "+", treatment, ":", modifier,
                           "+sex + GAenroll + enrollage + APdichenroll + educdich + wealthcat"))
      } else if (strata_name == "all" && length(modifier_counts) > 1 && all(modifier_counts > 3)) {
        data_stratified <- data_stratified[!is.na(data_stratified[[modifier]]), ]
        formula <-
          as.formula(paste(outcome, "~", treatment, "+", modifier, "+", treatment, ":", modifier,
                           "+sex + gravidity_cat + GAenroll + enrollage + APdichenroll + educdich + wealthcat"))
      } else {
        message("One of the levels doesn't have enough variability. Skipping analysis.")
        return(NULL)
      }
      print(formula)
      
      if (outcome_type == "continuous") {
        model <- glm(formula, data = data_stratified, family = gaussian(link = "identity"), na.action = na.omit)
        
      } else if (outcome_type == "binary") {
        model <- glm(formula, data = data_stratified, family = poisson(link = "log"), na.action = na.omit)
      }
      
      N = data_stratified %>%
        filter(!is.na(outcome)) %>% 
        distinct(id) %>% count()
      
      # model without interaction tx by gravidity
      estimates <- summary(model)$coefficients
      coefs = names(coef(model))
      treatment_row = as.numeric(grep("armDP$", coefs))
      modifier_row = as.numeric(grep(paste0("^", modifier), coefs))
      modifier_levels = grep(paste0("^", modifier), coefs, value = TRUE)
      interaction_row = as.numeric(grep(paste0("^", treatment, "DP:", modifier), coefs))
      
      point_estimate <- estimates[treatment_row, 1]
      point_estimate_modified <- ifelse(outcome_type == "continuous", point_estimate, exp(point_estimate))
      SE <- estimates[treatment_row, 2]
      lower_CI <- ifelse(outcome_type == "continuous", point_estimate - qnorm(0.975) * SE, exp(point_estimate - qnorm(0.975) * SE))
      upper_CI <- ifelse(outcome_type == "continuous", point_estimate + qnorm(0.975) * SE, exp(point_estimate + qnorm(0.975) * SE))
      
      # Assessment of effect modification
      if (modifier != "NA" && length(modifier_counts) > 1 && all(modifier_counts > 3)) {
        beta2 = coef(model)[modifier_row]
        beta3 = coef(model)[interaction_row]
        # Additive interaction level
        RERI = 1 + exp(beta3 + beta2 + point_estimate) - exp(point_estimate) - exp(beta2)
        addi_g <-
          as.formula(~ exp(x3 + x2 + x1) + 1 - exp(x1) - exp(x2))
        RERI_SE = deltamethod(g = addi_g,
                              mean = coef(model)[c(treatment_row, modifier_row, interaction_row)],
                              cov = vcov(model)[c(treatment_row, modifier_row, interaction_row), c(treatment_row, modifier_row, interaction_row)])
        RERI_CI <-
          cbind(RERI - qnorm(0.975) * RERI_SE, RERI + qnorm(0.975) * RERI_SE)
        # Relative (multiplicative) interaction level
        rel_int <- exp(beta3)
        rel_g = as.formula(~ exp(x1))
        rel_SE = deltamethod(g = rel_g,
                             mean = coef(model)[interaction_row],
                             cov = vcov(model)[interaction_row, interaction_row])
        rel_CI <-
          cbind(rel_int - qnorm(0.975) * rel_SE,
                rel_int + qnorm(0.975) * rel_SE)
        modifier_levels = levels(data_stratified[[modifier]])[c(2:length(levels(data_stratified[[modifier]])))]
      } else {
        RERI <- NA
        rel_int <- NA
        RERI_CI <- matrix(NA, nrow = 1, ncol = 2)
        rel_CI <- matrix(NA, nrow = 1, ncol = 2)
        modifier_levels = NA
      }
      
      outcome_remark = case_when(
        outcome == "haz"| outcome == "haz_quarter" ~ "height-for-age z score",
        outcome == "whz"| outcome == "whz_quarter" ~ "weight-for-height z score",
        outcome == "waz"| outcome == "waz_quarter" ~ "weight-for-age z score",
        outcome == "wgv1" ~ "absolute weight growth velocity (1 month increment)",
        outcome == "wgv3" ~ "absolute weight growth velocity (3 month increment)",
        outcome == "wlz_gv1" ~ "weight velocity based on weight-for-length Z-score (1 month increment)",
        outcome == "wlz_gv3" ~ "weight velocity based on weight-for-length Z-score (3 month increment)",
        outcome == "lgv2" ~ "absolute length growth velocity (2 month increment)",
        outcome == "lgv3" ~ "absolute length growth velocity (3 month increment)",
        outcome == "laz_gv2" ~ "length velocity based on length-for-age Z-score (2 month increment)",
        outcome == "laz_gv3" ~ "length velocity based on length-for-age Z-score (3 month increment)", 
        outcome== "haz_ms_stunt" | outcome == "haz_ms_stunt_quarter" ~ "prevalence: moderate to severe stunting",
        outcome == "haz_s_stunt" | outcome== "haz_s_stunt_quarter" ~ "prevalence: severe stunting",
        outcome == "whz_ms_waste" | outcome == "whz_ms_waste_quarter" ~ "prevalence: moderate to severe wasting",
        outcome == "whz_s_waste"| outcome == "whz_s_waste_quarter" ~ "prevalence: severe wasting",
        outcome == "waz_underwt"| outcome== "waz_underwt_quarter" ~ "prevalence: underweight",
        grepl("SGA", outcome) ~ "small for gestational weight",
        grepl("incident_haz_ms_stunt", outcome) ~ "incidence: moderate to severe stunting",
        grepl("incident_haz_s_stunt", outcome) ~ "incidence: severe stunting",
        grepl("incident_whz_ms_waste", outcome) ~ "incidence: moderate to severe wasting",
        grepl("incident_whz_s_waste", outcome) ~ "incidence: severe wasting",
        grepl("incident_waz_underwt", outcome) ~ "incidence: underweight",
        TRUE ~ NA
      )
      
      # Build output data frame
      output <- data.frame(
        outcome = outcome,
        outcome_remark = outcome_remark,
        time_unit = time_unit,
        age_group = age_group,
        gravidity_strata = strata_name,
        N_from_analysis = N,
        point_estimate = point_estimate_modified,
        SE = SE,
        lower_95CI = lower_CI,
        upper_95CI = upper_CI,
        modifier_name = ifelse(modifier=="NA", NA, modifier),
        modifier_level = modifier_levels,
        additive_effect = RERI,
        additive_CI = RERI_CI,
        relative_effect = rel_int,
        relative_CI = rel_CI
      )
      
      rownames(output) <- NULL
      output_list[[i]] = output
      
    }
    
    output_df = bind_rows(output_list)
    return(output_df)
  }


#--------------------------------------------------------
# apply the function and generate results
#--------------------------------------------------------

full_modifier_list = c("NA", "anymalaria", "maternal_agecat", "mother_heightcat", "wet_season")

crossing_zscore_1mo = crossing(
  outcome = outcome_zscore,
  time_unit = "1 month",
  age_group = age_list_1mo_round,
  age_category = "agemonth_round") 


crossing_zscore_3mo = crossing(
  outcome = outcome_zscore_quarter,
  time_unit = "3 month",
  age_group = age_list_3mo_birth,
  age_category = "agecat_birth"
) %>%
  mutate(age_group = as.character(age_group))


# --------------------------------------------------------

aim1_application <-
  function(data_set,
           crossing_set,
           modifier_list,
           outcome_type) {
    full_results_list = list()
    
    for (i in 1:nrow(crossing_set)) {
      results_list = lapply(modifier_list, function(x)
        aim1_analysis(
          outcome = as.character(crossing_set[i, "outcome"]),
          data = data_set,
          time_unit = as.character(crossing_set[i, "time_unit"]),
          age_category = as.character(crossing_set[i, "age_category"]),
          age_group = as.character(crossing_set[i, "age_group"]),
          outcome_type = outcome_type,
          treatment = "Txarm",
          modifier = x
        ))
      full_results_list = append(full_results_list, results_list)
    }
    
    combined_results_df <- bind_rows(full_results_list)
    combined_results_df = data.frame(lapply(combined_results_df, function(x) if(is.numeric(x)) round(x, 4) else x))
    
    View(combined_results_df)
    return(combined_results_df)
  }


#--------------------------------------------------------
# apply the wrapper function and save the results
#--------------------------------------------------------

z_score_monthly_results_stratified = aim1_application(data_set = data_zscore_monthly, crossing_set = crossing_zscore_1mo, modifier_list = full_modifier_list, outcome_type = "continuous")
View(z_score_monthly_results_stratified)
saveRDS(z_score_monthly_results_stratified, paste0(results_path,"aim1-stratified/aim1_zscore_monthly_results_stratified.RDS"))

z_score_quarterly_results_stratified = aim1_application(data_set = data_zscore_quarterly, crossing_set = crossing_zscore_3mo, modifier_list = full_modifier_list, outcome_type = "continuous")
View(z_score_quarterly_results_stratified)
saveRDS(z_score_quarterly_results_stratified, paste0(results_path,"aim1-stratified/aim1_zscore_quarterly_results_stratified.RDS"))

