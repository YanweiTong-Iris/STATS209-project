################################################################
# IPTp and child growth
# Script for intervention-mediator and mediator-outcome analysis
# Last updated: June 4, 2023
################################################################

rm(list = ls())
source(paste0(here::here(), "/0-config.R"))

#--------------------------------------------------------
# load data
#--------------------------------------------------------

data_continuous = readRDS(paste0(here::here(), "/0-input-data/analysis_data_continuous.RDS"))
data_monthly_ceiling = readRDS(paste0(here::here(), "/0-input-data/analysis_data_monthly_ceiling.RDS"))
data_zscore_monthly = readRDS(paste0(here::here(), "/0-input-data/analysis_data_zscore_monthly.RDS"))
data_zscore_quarterly = readRDS(paste0(here::here(), "/0-input-data/analysis_data_zscore_quarterly.RDS")) %>% 
  mutate(birthweight_01kg = birthweight_kg * 10)
#data_prevalence_quarterly = readRDS(paste0(here::here(), "/0-input-data/analysis_data_prevalence_quarterly.RDS"))
data_incidence_3month = readRDS(paste0(here::here(), "/0-input-data/analysis_data_incidence_quarterly.RDS"))
data_incidence_6month = readRDS(paste0(here::here(), "/0-input-data/analysis_data_incidence_biannual.RDS"))
data_incidence_12month = readRDS(paste0(here::here(), "/0-input-data/analysis_data_incidence_annual.RDS"))


#--------------------------------------------------------
# create outcome, age, and mediator list 
#--------------------------------------------------------

outcome_monthly_round_continuous = c("haz", "whz", "waz")
outcome_monthly_round_binary = c("haz_ms_stunt", "haz_s_stunt", 
                          "whz_ms_waste", "whz_s_waste", "waz_underwt")

outcome_monthly_ceiling_continuous = c("haz", "whz", "waz") 
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
age_list_2mo = factor(c("<2 months", "2- 4 months", "4- 6 months", "6- 8 months", "8- 10 months", "10- 12 months"), levels = c("<2 months", "2- 4 months", "4- 6 months", "6- 8 months", "8- 10 months", "10- 12 months"))
age_list_3mo = factor(c("0-3 months", ">3-6 months", ">6-9 months", ">9-12 months"), levels = c("0-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))
age_list_3mo_birth = factor(c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"), levels= c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))
age_list_6mo_birth = factor(c("birth", "1 day- <6 months", "6-12 months"), levels = c("birth", "1 day- <6 months", "6-12 months"))
age_list_12mo_birth = factor(c("birth", "1 day- 12 months"), levels =c("birth", "1 day- 12 months"))

# can include "anyhb_28", "anyhb_36" as continuous measures
binary_mediators = c("LBW", "preterm", "anyHP", "placentalLAMPdich", 
                     "placentalmal", "anemia_28binary", "anemia_36binary",
                     "antibacterial_binary", "betalactam_binary",
                     "fluoroquinolone_binary", "sulfonamide_binary",
                     "antibacterial_other_binary", "antimalarial_binary", 
                     "antiparasitic_binary", "antiviral_binary")
# anymalaria varies with infant age, will be addressed individually (monthly only)

# anymed_binary, any_nonmalaria_med_binary lack of variability
maternal_mediators = c("LBW", "preterm", "placentalmal", "anemia_28binary", 
                       "anemia_36binary", "gestational_weightchange", "SCF",
                       "antibacterial_binary", "betalactam_binary",
                       "fluoroquinolone_binary", "sulfonamide_binary",
                       "antibacterial_other_binary", "antimalarial_binary", 
                       "antiparasitic_binary", "antiviral_binary")

maternal_mediators_main = c("anemia_28binary", "gestational_weightchange",
                            "antibacterial_binary", "betalactam_binary",
                            "placentalmal", "IL4", "IL10", "SCF",
                             "TRANCE", "TNF", "IFN_gamma",
                            "preterm", "LBW","birthweight", "birthlength", 
                            "birthweight_kg", "birthweight_01kg")

maternal_mediators_2.1 = c("anemia_28binary", "gestational_weightchange",
                           "antibacterial_binary", "betalactam_binary",
                           "placentalmal", "SCF","IL4", "IL10")

preceding_malaria_mediator = c("prec_malaria_stunt", "prec_malaria_waste")

#--------------------------------------------------------
# intervention-mediator & mediator-outcome analysis wrapper function
#--------------------------------------------------------

# Documentation: mediator_analysis
# Usage: mediator_analysis(data, time_unit, age_group, model, independent_var, dependent_var, dependent_type)
# Description: intervention-mediator & mediator-outcome analysis
# Intervention: Txarm OR rand_Txarm (for blind testing)
# Mediators of interest: infant malaria (incidentmalaria), 
#                        low birth weight (LBW), preterm birth (preterm), 
#                        placental malaria (anyHP), non-malarial prenatal infection (not included), 
#                        anti-microbial use (see above), maternal anemia (anemia_28binary, anemia_36binary),
#                        gestational weight gain (gestational_weightchange), 
#                        and maternal inflammation (SCF only). 
# In intervention-mediator models, use GLM with a Gaussian family and identity link for continuous mediators 
# and a Poisson family and log link for dichotomous mediators. All intervention-mediator models will be unadjusted. 
# For mediator-outcome models, use GEE with a Gaussian family, identity link, and exchangeable correlation matrix 
# for continuous outcomes since there were repeated measures within each age category. 
# For dichotomous outcomes, use a Poisson family and log link since there are no repeated measures. 
# Mediator-outcome models will be adjusted for confounders of mediators and outcomes by choosing nodes 
# sufficient to block backdoor pathways

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
# age_group: to subset the full dataset to a single age range
# model: "intervention-mediator" OR "mediator-outcome"
# independent_var: x of the formula (as a string)
# dependent_var: y of the formula (as a string)
# variable to choose: 
# intervention: "Txarm", "rand_Txarm"
# mediator: "anymalaria", "LBW", "preterm", "anyHP", "anemia_28", “anemia_36", 
#             "gestational_weightchange", "placentalLAMPdich"
#           "antibacterial", "betalactam", "fluoroquinolone", "sulfonamide", "antibacterial_other",
#           "antimalarial", "antiparasitic", "antiviral", "anymed"
# confounders: "sex", "gravidity_cat", "APdichenroll", "enrollage_binary", "educdich", "wealth_binary"
# outcome:
# - continuous: "haz", "whz", "waz", "wgv1", "wlz_gv1", "wgv3", "wlz_gv3",
#             "lgv2", "laz_gv2"， "lgv3", "laz_gv3"
# - binary: "SGA", haz_ms_stunt", "haz_s_stunt", "whz_ms_waste", "whz_s_waste", 
#         "waz_underweight", "incident_haz_ms_stunt", "incident_haz_s_stunt", 
#         "incident_whz_ms_waste", "incident_whz_s_waste", "incident_waz_underwt"
# dependent_type: "continuous" OR "binary" (as a string)

# Returns: the dataframe as described above
# Output: prints the data frame described above

mediator_analysis <- function(data, time_unit, age_category, age_group, model, independent_var, dependent_var, dependent_type, gravidity_strata) {
  print(paste0("model: ", model, ", x: ", independent_var, ", y: ", dependent_var, ", age: ", age_group))
  
  comment({#when mediator = anymalaria, for infant ever infected, only look at data post first infection
  if (independent_var == "new_anymalaria"){
    data = data %>% 
      group_by(id) %>%
      mutate(post_malaria = cumsum(as.numeric(as.character(new_anymalaria))) > 0 | max(as.numeric(as.character(new_anymalaria))) == 0) %>%
      ungroup()
    data = data[data$post_malaria,]
  }})
  
  data_age <- data[data[[age_category]] == age_group, ]
  data_age <- data_age[!is.na(data_age[[dependent_var]]), ]
  print(nrow(data_age))
  
  confounder <- NA
  adjusted <- 0
  
  # output_list = list()
  # strata_list = list(
  #   "all" = c("single", "multi"),
  #   "single" = c("single"),
  #   "multi" = c("multi")
  # )
  # 
  # for (i in 1:3) {
  #   strata_name = names(strata_list)[[i]]
  #   strata = strata_list[[strata_name]]
  if(gravidity_strata == "all"){
    data_stratified = data_age
  }else {
    data_stratified = data_age %>% filter(gravidity_cat == gravidity_strata)
    }
    
   
    # IM analysis now in file 7a
    # else if (model == "intervention-mediator") {
    #   data_age[[dependent_var]] <-
    #     as.numeric(as.character(data_age[[dependent_var]]))
    #   formula <-
    #     as.formula(paste(dependent_var, "~", independent_var))
    # } 
    

    # Compose the formula
    # if (model == "mediator-outcome" &
    #     independent_var %in% c("preterm", "LBW", "birthweight", "birthlength", "gestational_weightchange")) {
       adjusted <- 1
       confounder <- "full"
      if (gravidity_strata == "all") {
        formula <-
          as.formula(paste(dependent_var, "~",independent_var,
             "+ Txarm + sex + gravidity_cat+ GAenroll + enrollage + APdichenroll + educdich + wealthcat"))
      } else{
        formula <-
          as.formula(paste(dependent_var,"~",independent_var,
              "+ Txarm + sex + GAenroll + enrollage + APdichenroll + educdich + wealthcat"))
      }
    # } else {
    #   adjusted <- 1
    #   confounder <- "full- sex"
    #   if (gravidity_strata == "all") {
    #     formula <-as.formula(paste(dependent_var,"~", independent_var,
    #           "+ Txarm + gravidity_cat+ GAenroll + enrollage + APdichenroll + educdich + wealthcat"))
    #   } else{
    #     formula <-as.formula(paste(dependent_var,"~",independent_var,
    #                                "+ Txarm + GAenroll + enrollage + APdichenroll + educdich + wealthcat"))
    #   }
    # }
    print(formula)
  
    # Check if the outcome has enough variation
    if (dependent_type == "binary") {
      counts = table(data_stratified[[dependent_var]])
      print(counts)
      if (!(length(unique(data_stratified[[dependent_var]])) >= 2 && all(counts >= 10))) {
        message("Variable '", dependent_var, "' doesn't have enough variability. Skipping analysis.")
        return(NULL)
      }
    }

    if(independent_var %in% c("anemia_28binary", "antibacterial_binary",
                             "betalactam_binary", "placentalmal",
                             "preterm", "LBW", "anymalaria", "new_anymalaria")){
      counts = table(data_stratified[[independent_var]])
      print(counts)
      if (!(length(unique(data_stratified[[independent_var]])) >= 2 && all(counts >= 10))) {
        message("Variable '", independent_var, "' doesn't have enough variability. Skipping analysis.")
        return(NULL)
      }
    }
  
  
  # Run GLM model
  if (dependent_type == "continuous") {
    model.fit <- glm(formula, data = data_stratified, family = gaussian(link = "identity"), na.action = na.omit)
  } else if (dependent_type == "binary") {
    data_stratified <- data_stratified[!is.na(data_stratified[[independent_var]]), ]
    model.fit <- glm(formula, data = data_stratified, family = poisson(link = "log"), na.action = na.omit)
  }
  #print(summary(model.fit))
  
  N <- data_stratified %>% filter(!is.na(dependent_var)) %>% distinct(id) %>% count()
  N_tx <- data_stratified %>% 
    group_by(Txarm) %>% 
    summarise(N = n())
  
  N_DP <- N_tx$N[N_tx$Txarm == "DP"]
  N_SP <- N_tx$N[N_tx$Txarm == "SP"]
  
  estimates <- summary(model.fit)$coefficients
  coefs <- names(coef(model.fit))
  independent_var_row <- as.numeric(grep(paste0("^", independent_var), coefs))
  
  point_estimate <- estimates[independent_var_row, 1]
  point_estimate_modified <- ifelse(dependent_type == "continuous", point_estimate, exp(point_estimate))
  SE <- estimates[independent_var_row, 2]
  lower_CI <- ifelse(dependent_type == "continuous", point_estimate - qnorm(0.975) * SE, exp(point_estimate - qnorm(0.975) * SE))
  upper_CI <- ifelse(dependent_type == "continuous", point_estimate + qnorm(0.975) * SE, exp(point_estimate + qnorm(0.975) * SE))
  
  outcome_remark <- case_when(
    dependent_var == "haz_quarter" | dependent_var == "haz" ~ "height-for-age z score",
    dependent_var == "whz_quarter" | dependent_var == "whz"~ "weight-for-height z score",
    dependent_var == "waz_quarter" | dependent_var == "waz" ~ "weight-for-age z score",
    dependent_var == "wgv1" ~ "absolute weight growth velocity (1 month increment)",
    dependent_var == "wgv3" ~ "absolute weight growth velocity (3 month increment)",
    dependent_var == "wlz_gv1" ~ "weight velocity based on weight-for-length Z-score (1 month increment)",
    dependent_var == "wlz_gv3" ~ "weight velocity based on weight-for-length Z-score (3 month increment)",
    dependent_var == "lgv2" ~ "absolute length growth velocity (2 month increment)",
    dependent_var == "lgv3" ~ "absolute length growth velocity (3 month increment)",
    dependent_var == "laz_gv2" ~ "length velocity based on length-for-age Z-score (2 month increment)",
    dependent_var == "laz_gv3" ~ "length velocity based on length-for-age Z-score (3 month increment)",
    grepl("SGA", dependent_var) ~ "small for gestational weight",
    dependent_var == "haz_ms_stunt_quarter" | dependent_var == "haz_ms_stunt" ~ "prevalence: moderate to severe stunting",
    dependent_var == "haz_s_stunt_quarter"| dependent_var == "haz_s_stunt" ~ "prevalence: severe stunting",
    dependent_var == "whz_ms_waste_quarter" | dependent_var == "whz_ms_waste" ~ "prevalence: moderate to severe wasting",
    dependent_var == "whz_s_waste_quarter" | dependent_var == "whz_s_waste" ~ "prevalence: severe wasting",
    dependent_var == "waz_underwt_quarter" | dependent_var == "waz_underwt" ~ "prevalence: underweight",
    grepl("incident_haz_ms_stunt", dependent_var) ~ "incidence: moderate to severe stunting",
    grepl("incident_haz_s_stunt", dependent_var) ~ "incidence: severe stunting",
    grepl("incident_whz_ms_waste", dependent_var) ~ "incidence: moderate to severe wasting",
    grepl("incident_whz_s_waste", dependent_var) ~ "incidence: severe wasting",
    grepl("incident_waz_underwt", dependent_var) ~ "incidence: underweight",
    TRUE ~ NA
  )
  
  mediator_remark <- case_when(
    independent_var == "anemia_28binary" | dependent_var == "anemia_28binary" ~ "Anemia at gestational week 28",
    independent_var == "anemia_36binary" | dependent_var == "anemia_36binary" ~ "Anemia at gestational week 36",
    independent_var == "anyHP" | dependent_var == "anyHP" ~ "Placental malaria",
    grepl("anymalaria", independent_var) | grepl("anymalaria", dependent_var) ~ "Infant malaria",
    independent_var == "LBW" | dependent_var == "LBW" ~ "Low birth weight",
    independent_var == "preterm" | dependent_var == "preterm" ~ "Pre-term birth",
    independent_var == "gestational_weightchange" | dependent_var == "gestational_weightchange" ~ "Gestational weight change (kg)",
    independent_var == "placentalLAMPdich" | dependent_var == "placentalLAMPdich" ~ "Placental malaria (blood)",
    independent_var == "placentalmal" | dependent_var == "placentalmal" ~ "Placental malaria (blood or histopath)",
    independent_var == "antibacterial_binary" | dependent_var == "antibacterial_binary" ~ "Any antibacterial use",
    independent_var == "betalactam_binary" | dependent_var == "betalactam_binary" ~ "Beta-lactam use",
    independent_var == "SCF" | dependent_var == "betalactam_binary" ~ "Stem cell factor",
    independent_var == "birthweight" | dependent_var == "birthweight" ~ "Birth weight",
    independent_var == "birthweight_kg" | dependent_var == "birthweight_kg" ~ "Birth weight (kg)",
    independent_var == "birthlength" | dependent_var == "birthlength" ~ "Birth length",
    TRUE ~ NA
  )
  
  # Build output data frame
  if (model == "intervention-mediator") {
    output <- data.frame(
      model = model,
      gravidae = gravidity_strata,
      adjusted = adjusted,
      confounder = confounder,
      independent_variable = independent_var,
      dependent_variable = dependent_var,
      dependent_type = dependent_type,
      mediator_remark = mediator_remark,
      time_unit = time_unit,
      age_group = age_group,
      N_from_analysis = N,
      N_DP = N_DP,
      N_SP = N_SP,
      point_estimate = point_estimate_modified,
      SE = SE,
      lower_95CI = lower_CI,
      upper_95CI = upper_CI
    )
  } else {
    output <- data.frame(
      model = model,
      gravidae = gravidity_strata,
      adjusted = adjusted,
      confounder = confounder,
      independent_variable = independent_var,
      dependent_variable = dependent_var,
      dependent_type = dependent_type,
      mediator_remark = mediator_remark,
      outcome_remark = outcome_remark,
      time_unit = time_unit,
      age_group = age_group,
      N_from_analysis = N,
      N_DP = N_DP,
      N_SP = N_SP,
      point_estimate = point_estimate_modified,
      SE = SE,
      lower_95CI = lower_CI,
      upper_95CI = upper_CI
    )
  }
  rownames(output) <- NULL
  #output_list[[i]] = output
  #output_df = bind_rows(output_list)
  #print(output_df)
  #return(output_df)
  return(output)
}


#--------------------------------------------------------
# function application
#--------------------------------------------------------

mediator_analysis_application <-
  function(data_set,
           crossing_set,
           model) {
    full_results_list = list()
    
    for (i in 1:nrow(crossing_set)) {
      try(full_results_list[[i]] <- mediator_analysis(
        data = data_set,
        time_unit = as.character(crossing_set[i, "time_unit"]),
        age_category = as.character(crossing_set[i, "age_category"]),
        age_group =  as.character(crossing_set[i, "age_group"]),
        gravidity_strata = as.character(crossing_set[i, "gravidity_strata"]),
        model = model,
        independent_var = as.character(crossing_set[i, "independent_var"]),
        dependent_var = as.character(crossing_set[i, "dependent_var"]),
        dependent_type = as.character(crossing_set[i, "dependent_type"])))
    }
    
    combined_results_df <- do.call(rbind, full_results_list)
    combined_results_df = data.frame(lapply(combined_results_df, function(x) if(is.numeric(x)) round(x, 4) else x))
    View(combined_results_df)
    return(combined_results_df)
  }


#--------------------------------------------------------
# create crossing sets 
#--------------------------------------------------------

# IM = intervention-mediator
# MO = mediator-outcome

# anymalaria at birth removed due to low variability

# crossing_malaria_outcome_MO_round = rbind(
#   crossing(
#     independent_var = "new_anymalaria",
#     dependent_var = outcome_monthly_round_continuous,
#     time_unit = "1 month",
#     age_group = c(0,1,2,3,4,5,6,7,8,9,10,11,12),
#     age_category = "agemonth_round",
#     dependent_type = "continuous"
#   ), crossing(
#     independent_var = "new_anymalaria",
#     dependent_var = outcome_monthly_round_binary,
#     time_unit = "1 month",
#     age_group = c(0,1,2,3,4,5,6,7,8,9,10,11,12),
#     age_category = "agemonth_round",
#     dependent_type = "binary"
#   ))
# 
crossing_malaria_outcome_MO_ceiling = rbind(crossing(
    independent_var = "new_anymalaria",
    dependent_var = outcome_monthly_ceiling_continuous,
    time_unit = "1 month",
    age_group = factor(c(1,2,3,4,5,6,7,8,9,10,11,12)),
    age_category = "agemonth_ceiling",
    dependent_type = "continuous",
    gravidity_strata = c("all", "single", "multi")),
    crossing(
      independent_var = "anymalaria",
      dependent_var = outcome_monthly_ceiling_continuous,
      time_unit = "1 month",
      age_group = factor(c(1,2,3,4,5,6,7,8,9,10,11,12)),
      age_category = "agemonth_ceiling",
      dependent_type = "continuous",
      gravidity_strata = c("all", "single", "multi"))
    )
  
  # ),
  # crossing(
  #   independent_var = "anymalaria",
  #   dependent_var = outcome_monthly_ceiling_binary,
  #   time_unit = "1 month",
  #   age_group = factor(c(1,2,3,4,5,6,7,8,9,10,11,12)),
  #   age_category = "agemonth_ceiling",
  #   dependent_type = "binary",
  #   gravidity_strata = c("all", "single", "multi"),
  # ))

crossing_preceding_malaria = crossing(
  independent_var = preceding_malaria_mediator,
  dependent_var = c("incident_haz_ms_stunt_agecat_birth",
              "incident_whz_ms_waste_agecat_birth"),
  time_unit = "3 month",
  age_group = age_list_3mo_birth,
  age_category = "agecat_birth",
  dependent_type = "binary",
  gravidity_strata = c("all", "single", "multi")
) %>% mutate(age_group = as.character(age_group))

crossing_zscore_MO_3mo =crossing(
    independent_var = maternal_mediators_main,
    dependent_var = outcome_z_score_quarter,
    time_unit = "3 month",
    age_group = age_list_3mo_birth,
    age_category = "agecat_birth",
    gravidity_strata = c("all", "single", "multi"),
    dependent_type = "continuous"
  ) %>% mutate(age_group = as.character(age_group))

crossing_prevalence_MO_3mo = crossing(
    independent_var = maternal_mediators_main,
    dependent_var = outcome_prevalence_quarterly,
    time_unit = "3 month",
    age_group = age_list_3mo_birth,
    age_category = "agecat_birth",
    gravidity_strata = c("all", "single", "multi"),
    dependent_type = "binary"
  ) %>% mutate(age_group = as.character(age_group))

crossing_velocity_1mo = crossing(
  independent_var = maternal_mediators_main,
  dependent_var = outcome_velocity_1mo,
  time_unit = "1 month",
  age_group = age_list_1mo_ceiling,
  age_category = "agemonth_ceiling",
  gravidity_strata = c("all", "single", "multi"),
  dependent_type = "continuous"
) %>% mutate(age_group = as.character(age_group))

crossing_velocity_2mo = crossing(
  independent_var = maternal_mediators_main,
  dependent_var = outcome_velocity_2mo,
  time_unit = "2 month",
  age_group = age_list_2mo,
  age_category = "age_2_month",
  gravidity_strata = c("all", "single", "multi"),
  dependent_type = "continuous"
) %>% mutate(age_group = as.character(age_group))

crossing_velocity_3mo = crossing(
  independent_var = maternal_mediators_main,
  dependent_var = outcome_velocity_3mo,
  time_unit = "3 month",
  age_group = age_list_3mo,
  age_category = "agecat",
  gravidity_strata = c("all", "single", "multi"),
  dependent_type = "continuous"
) %>% mutate(age_group = as.character(age_group))

crossing_incidence_3mo = rbind(
  crossing(independent_var = maternal_mediators_main,
           dependent_var = outcome_incidence_3mo,
           time_unit = "3 month",
           age_group = age_list_3mo_birth,
           age_category = "agecat_birth",
           gravidity_strata = c("all", "single", "multi"),
           dependent_type = "binary"),
  crossing(independent_var = maternal_mediators_main,
           dependent_var = "SGA_quarterly",
           time_unit = "3 month",
           age_group = "Birth",
           age_category = "agecat_birth",
           gravidity_strata = c("all", "single", "multi"),
           dependent_type = "binary")) %>%
  mutate(age_group = as.character(age_group))

crossing_incidence_6mo = rbind(
  crossing(independent_var = maternal_mediators_main,
           dependent_var = outcome_incidence_6mo,
           time_unit = "6 month",
           age_group = age_list_6mo_birth,
           age_category = "age_6_12_month",
           gravidity_strata = c("all", "single", "multi"),
           dependent_type = "binary"),
  crossing(independent_var = maternal_mediators_main,
           dependent_var = "SGA_biannual",
           time_unit = "6 month",
           age_group = "birth",
           age_category = "age_6_12_month",
           gravidity_strata = c("all", "single", "multi"),
           dependent_type = "binary")) %>%
  mutate(age_group = as.character(age_group))

crossing_incidence_12mo = rbind(
  crossing(independent_var = maternal_mediators_main,
           dependent_var = outcome_incidence_12mo,
           time_unit = "12 month",
           age_group = age_list_12mo_birth,
           age_category = "age_1_12",
           gravidity_strata = c("all", "single", "multi"),
           dependent_type = "binary"),
  crossing(independent_var = maternal_mediators_main,
           dependent_var = "SGA_annual",
           time_unit = "12 month",
           age_group = "birth",
           age_category = "age_1_12",
           gravidity_strata = c("all", "single", "multi"),
           dependent_type = "binary")) %>%
  mutate(age_group = as.character(age_group))

crossing_2.1 = rbind(
  crossing(independent_var = maternal_mediators_2.1,
           dependent_var = c("LBW", "preterm"),
           time_unit = "3 month",
           age_group = "Birth",
           age_category = "agecat_birth",
           gravidity_strata = c("all", "single", "multi"),
           dependent_type = "binary"),
  crossing(independent_var = maternal_mediators_2.1,
           dependent_var = c("birthlength", "birthweight_kg"),
           time_unit = "3 month",
           age_group = "Birth",
           age_category = "agecat_birth",
           gravidity_strata = c("all", "single", "multi"),
           dependent_type = "continuous")) %>%
  mutate(age_group = as.character(age_group))

#--------------------------------------------------------
# save the results
#--------------------------------------------------------

#prevent using scientific notations
options(scipen = 999)

# ! IM results are now analyzed in 7a-intervention-mediator-analysis
# IM_results = mediator_analysis_application(data_set = data_zscore_quarterly, crossing_set = crossing_IM, model = "intervention-mediator")
# View(IM_results)
# saveRDS(IM_results, paste0(results_path,"intervention_mediator_results.RDS"))

MO_malaria_outcome_1mo_ceiling = mediator_analysis_application(data_set = data_monthly_ceiling, crossing_set = crossing_malaria_outcome_MO_ceiling, model = "mediator-outcome")
# MO_malaria_outcome_1mo_round = mediator_analysis_application(data_set = data_monthly_round, crossing_set = crossing_malaria_outcome_MO_round, model = "mediator-outcome")
# MO_malaria_outcome_1mo = rbind(MO_malaria_outcome_1mo_ceiling, MO_malaria_outcome_1mo_round)
# View(MO_malaria_outcome_1mo)
# saveRDS(MO_malaria_outcome_1mo, paste0(results_path,"mediator_outcome_malaria_outcome_1mo.RDS"))

MO_zscore_3mo_stratified = mediator_analysis_application(data_set = data_zscore_quarterly, crossing_set = crossing_zscore_MO_3mo, model = "mediator-outcome")
View(MO_zscore_3mo_stratified)
saveRDS(MO_zscore_3mo_stratified, paste0(here::here(), "/4-result-data/mediator_outcome_zscore_results_3mo_stratified.RDS"))

MO_prevalence_3mo_stratified = mediator_analysis_application(data_set = data_prevalence_quarterly, crossing_set = crossing_prevalence_MO_3mo, model = "mediator-outcome")
View(MO_prevalence_3mo_stratified)
saveRDS(MO_prevalence_3mo_stratified, paste0(results_path,"IM-MO-stratified/mediator_outcome_prevalence_results_3mo_stratified.RDS"))

MO_velocity_1mo_stratified = mediator_analysis_application(data_set = data_velocity_1month, crossing_set = crossing_velocity_1mo, model = "mediator-outcome")
View(MO_velocity_1mo_stratified)
saveRDS(MO_velocity_1mo_stratified, paste0(results_path,"IM-MO-stratified/mediator_outcome_velocity_results_1mo_stratified.RDS"))

MO_velocity_2mo_stratified = mediator_analysis_application(data_set = data_velocity_2month, crossing_set = crossing_velocity_2mo, model = "mediator-outcome")
View(MO_velocity_2mo_stratified)
saveRDS(MO_velocity_2mo_stratified, paste0(results_path,"IM-MO-stratified/mediator_outcome_velocity_results_2mo_stratified.RDS"))

MO_velocity_3mo_stratified = mediator_analysis_application(data_set = data_velocity_3month, crossing_set = crossing_velocity_3mo, model = "mediator-outcome")
View(MO_velocity_3mo_stratified)
saveRDS(MO_velocity_3mo_stratified, paste0(results_path,"IM-MO-stratified/mediator_outcome_velocity_results_3mo_stratified.RDS"))

MO_incidence_3mo_stratified = mediator_analysis_application(data_set = data_incidence_3month, crossing_set = crossing_incidence_3mo, model = "mediator-outcome")
View(MO_incidence_3mo_stratified)
saveRDS(MO_incidence_3mo_stratified, paste0(results_path,"IM-MO-stratified/mediator_outcome_incidence_results_3mo_stratified.RDS"))

MO_incidence_6mo_stratified = mediator_analysis_application(data_set = data_incidence_6month, crossing_set = crossing_incidence_6mo, model = "mediator-outcome")
View(MO_incidence_6mo_stratified)
saveRDS(MO_incidence_6mo_stratified, paste0(results_path,"IM-MO-stratified/mediator_outcome_incidence_results_6mo_stratified.RDS"))

MO_incidence_12mo_stratified = mediator_analysis_application(data_set = data_incidence_12month, crossing_set = crossing_incidence_12mo, model = "mediator-outcome")
View(MO_incidence_12mo_stratified)
saveRDS(MO_incidence_12mo_stratified, paste0(results_path,"IM-MO-stratified/mediator_outcome_incidence_results_12mo_stratified.RDS"))

MO_aim2_1_stratified = mediator_analysis_application(data_set = data_incidence_3month %>% mutate(preterm = as.numeric(as.character(preterm)),
                                                                                                 LBW = as.numeric(as.character(LBW))), 
                                                     crossing_set = crossing_2.1, model = "mediator-outcome")
View(MO_aim2_1_stratified)
saveRDS(MO_aim2_1_stratified, paste0(results_path,"IM-MO-stratified/mediator_outcome_aim2_1_stratified.RDS"))

MO_prec_malaria_3mo_stratified = mediator_analysis_application(data_set = data_incidence_3month, crossing_set = crossing_preceding_malaria, model = "mediator-outcome")
View(MO_prec_malaria_3mo_stratified)
saveRDS(MO_prec_malaria_3mo_stratified, paste0(results_path,"IM-MO-stratified/mediator_outcome_prec_malaria_3mo_stratified.RDS"))


# RESULT CHECK
check_data_binary = data_incidence_12month[data_incidence_12month$age_1_12 == "1 day- 12 months", ]
check_data_binary = check_data_binary[!is.na(check_data_binary$incident_whz_s_waste_age_1_12),]
model.fit1 = glm(
  incident_whz_s_waste_age_1_12 ~ incidentmalaria + sex+ gravidity_cat,
  data = check_data_binary,
  family = poisson(link = "log"),
  na.action = na.omit
)
summary(model.fit1)
  
check_data_continuous = data_zscore_quarterly[data_zscore_quarterly$agecat_birth == "1 day-3 months",]
check_data_continuous <- check_data_continuous[!is.na(check_data_continuous[["gestational_weightchange"]]),]
model.fit2 = glm(	
  haz_quarter~ gestational_weightchange,
  data = check_data_continuous,
  family = gaussian(link = "identity"))
summary(model.fit2)


ggplot(data_monthly_ceiling, aes(x = agemonth_ceiling, y = haz, color = as.factor(anymalaria), group = as.factor(anymalaria))) +
  stat_summary(fun = mean, geom = "point") +  
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "line", aes(group = as.factor(anymalaria))) +  
  labs(color = "Any Malaria Infection") +
  ggtitle("LAZ")+
  theme_minimal()

ggplot(data_monthly_ceiling, aes(x = agemonth_ceiling, y = haz, color = as.factor(new_anymalaria), group = as.factor(new_anymalaria))) +
  stat_summary(fun = mean, geom = "point") +  
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "line", aes(group = as.factor(new_anymalaria))) +  
  labs(color = "Imputed Malaria Infection") +
  ggtitle("LAZ")+
  theme_minimal()
