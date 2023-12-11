################################################################
# IPTp and child growth
# Script for intervention-mediator analysis
# Last updated: Oct 3, 2023
################################################################

rm(list = ls())
source(paste0(here::here(), "/0-config.R"))

#--------------------------------------------------------
# load data
#--------------------------------------------------------
data_zscore_quarterly = readRDS(paste0(data_path, "analysis_data_zscore_quarterly.RDS"))


#--------------------------------------------------------
# create outcome, age, and mediator list 
#--------------------------------------------------------

age_list_3mo_birth = factor(c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"), levels= c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))

binary_mediators = c("LBW", "preterm", "anyHP", "placentalLAMPdich", 
                     "placentalmal", "anemia_28binary", "anemia_36binary",
                     "antibacterial_binary", "betalactam_binary",
                     "fluoroquinolone_binary", "sulfonamide_binary",
                     "antibacterial_other_binary", "antimalarial_binary", 
                     "antiparasitic_binary", "antiviral_binary")

binary_mediators_main = c("anemia_28binary", "placentalmal", "preterm", "LBW")

continuous_mediator = c("gestational_weightchange", "IL4", "IL10", "SCF",
                        "TRANCE", "IFN_gamma", "birthweight_kg", "birthlength")

#--------------------------------------------------------
# intervention-mediator analysis wrapper function
#--------------------------------------------------------

intervention_mediator_analysis <- function(data, mediator, mediator_type) {
  print(paste0("mediator: ", mediator))
  
  data_age <- data[data[["agecat_birth"]] == "Birth", ]
  data_age <- data_age[!is.na(data_age[[mediator]]), ]
  data_age[[mediator]] <- as.numeric(as.character(data_age[[mediator]]))
  print(nrow(data_age))
  
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
    
    mean.risk <- data_stratified %>% group_by(Txarm) %>%
      summarise(value = if_else(
        mediator_type == "binary",
        sum(.data[[mediator]] == 1, na.rm = TRUE) / n(),
        mean(.data[[mediator]], na.rm = TRUE)
      ))
    
    mean.risk.intervention <- mean.risk$value[mean.risk$Txarm == "DP"]
    mean.risk.control <- mean.risk$value[mean.risk$Txarm == "SP"]
    
    
    unadjusted_formula <- as.formula(paste(mediator, "~ Txarm"))
    
    if (mediator %in% c("preterm", "LBW", "birthweight", "birthweight_kg", "birthlength", "gestational_weightchange")) {
      confounder <- "full"
      if (strata_name == "all"){
        adjusted_formula <- as.formula(paste(mediator,
              "~ Txarm + sex + gravidity_cat + GAenroll + enrollage + APdichenroll + educdich + wealthcat"))
      } else{
        adjusted_formula <- as.formula(paste(mediator,
                                             "~ Txarm + sex + GAenroll + enrollage + APdichenroll + educdich + wealthcat"))
      }
    } else{
      confounder <- "full- sex"
      if (strata_name == "all") {
        adjusted_formula <- as.formula(
          paste(
            mediator,
            "~ Txarm + gravidity_cat + GAenroll + APdichenroll + educdich + wealthcat"
          )
        )
      } else{
        adjusted_formula <- as.formula(paste(
          mediator,
          "~ Txarm + GAenroll + APdichenroll + educdich + wealthcat"
        ))
      }
    }
    
    # Run GLM model
    if (mediator_type == "continuous") {
      unadjusted.model.fit <-
        glm(
          unadjusted_formula,
          data = data_stratified,
          family = gaussian(link = "identity"),
          na.action = na.omit
        )
      adjusted.model.fit <-
        glm(
          adjusted_formula,
          data = data_stratified,
          family = gaussian(link = "identity"),
          na.action = na.omit
        )
    } else if (mediator_type == "binary") {
      unadjusted.model.fit <-
        glm(
          unadjusted_formula,
          data = data_stratified,
          family = poisson(link = "log"),
          na.action = na.omit
        )
      adjusted.model.fit <-
        glm(
          adjusted_formula,
          data = data_stratified,
          family = poisson(link = "log"),
          na.action = na.omit
        )
    }
    
    N <- data_stratified %>% filter(!is.na(mediator)) %>% distinct(id) %>% count()
    
    print(summary(unadjusted.model.fit))
    unadjusted_estimates <-
      summary(unadjusted.model.fit)$coefficients
    unadjusted_coefs <- names(coef(unadjusted.model.fit))
    unadjusted_Tx_row <-
      as.numeric(grep(paste0("^", "TxarmDP", "$"), unadjusted_coefs))
    unadjusted_point_estimate <-
      unadjusted_estimates[unadjusted_Tx_row, 1]
    unadjusted_point_estimate_modified <-
      ifelse(
        mediator_type == "continuous",
        unadjusted_point_estimate,
        exp(unadjusted_point_estimate)
      )
    unadjusted_SE <- unadjusted_estimates[unadjusted_Tx_row, 2]
    unadjusted_lower_CI <- ifelse(
      mediator_type == "continuous",
      unadjusted_point_estimate - qnorm(0.975) * unadjusted_SE,
      exp(unadjusted_point_estimate - qnorm(0.975) * unadjusted_SE)
    )
    unadjusted_upper_CI <- ifelse(
      mediator_type == "continuous",
      unadjusted_point_estimate + qnorm(0.975) * unadjusted_SE,
      exp(unadjusted_point_estimate + qnorm(0.975) * unadjusted_SE)
    )
    
    
    print(summary(adjusted.model.fit))
    adjusted_estimates <- summary(adjusted.model.fit)$coefficients
    adjusted_coefs <- names(coef(adjusted.model.fit))
    adjusted_Tx_row <-
      as.numeric(grep(paste0("^", "TxarmDP", "$"), adjusted_coefs))
    adjusted_point_estimate <-
      adjusted_estimates[adjusted_Tx_row, 1]
    adjusted_point_estimate_modified <-
      ifelse(
        mediator_type == "continuous",
        adjusted_point_estimate,
        exp(adjusted_point_estimate)
      )
    adjusted_SE <- adjusted_estimates[adjusted_Tx_row, 2]
    adjusted_lower_CI <- ifelse(
      mediator_type == "continuous",
      adjusted_point_estimate - qnorm(0.975) * adjusted_SE,
      exp(adjusted_point_estimate - qnorm(0.975) * adjusted_SE)
    )
    adjusted_upper_CI <- ifelse(
      mediator_type == "continuous",
      adjusted_point_estimate + qnorm(0.975) * adjusted_SE,
      exp(adjusted_point_estimate + qnorm(0.975) * adjusted_SE)
    )
    
    mediator_remark <- case_when(
      mediator == "anemia_28binary" ~ "Anemia at gestational week 28",
      mediator == "anemia_36binary" ~ "Anemia at gestational week 36",
      mediator == "anyHP" ~ "Placental malaria (histopath)",
      #grepl("anymalaria", independent_var) | grepl("anymalaria", dependent_var) ~ "infant malaria",
      mediator == "LBW" ~ "Low birthweight",
      mediator == "preterm" ~ "Pre-term birth",
      mediator == "gestational_weightchange" ~ "Gestational weight change (kg)",
      mediator == "placentalLAMPdich" ~ "Placental malaria (blood)",
      mediator == "placentalmal" ~ "Placental malaria (blood or histopath)",
      mediator == "antibacterial_binary" ~ "Any antibacterial use",
      mediator == "antiparasitic_binary" ~ "Any antiparasitic",
      mediator == "antimalarial_binary" ~ "Any antimalarial",
      mediator == "antiviral_binary"  ~ "Any antiviral",
      mediator == "betalactam_binary" ~ "Beta-lactam use",
      mediator == "fluoroquinolone_binary" ~ "Fluoroquinolone",
      mediator == "sulfonamide_binary" ~ "Sulfonamide",
      mediator == "antibacterial_other_binary" ~ "Other antibiotics",
      mediator == "SCF" ~ "Stem cell factor",
      mediator == "birthweight" ~ "Birth weight (g)",
      mediator == "birthweight_kg" ~ "Birth weight (kg)",
      mediator == "birthlength" ~ "Birth length (cm)",
      TRUE ~ NA
    )
    
    
    # Build output data frame
    output <- data.frame(
      mediator = mediator,
      mediator_type = mediator_type,
      mediator_remark = mediator_remark,
      confounder = confounder,
      gravidity_strata = strata_name,
      N_from_analysis = N,
      mean.risk.DP = mean.risk.intervention,
      mean.risk.SP = mean.risk.control,
      unadjusted_point_estimate = unadjusted_point_estimate_modified,
      unadjusted_SE = unadjusted_SE,
      unadjusted_lower_95CI = unadjusted_lower_CI,
      unadjusted_upper_95CI = unadjusted_upper_CI,
      adjusted_point_estimate = adjusted_point_estimate_modified,
      adjusted_SE = adjusted_SE,
      adjusted_lower_95CI = adjusted_lower_CI,
      adjusted_upper_95CI = adjusted_upper_CI
    )
    rownames(output) <- NULL
    output_list[[i]] = output
    
  }
  
  output_df = bind_rows(output_list)
  return(output_df)
}


#--------------------------------------------------------
# function application
#--------------------------------------------------------

intervention_mediator_application <-
  function(data_set,
           crossing_set) {
    full_results_list = list()
    
    for (i in 1:nrow(crossing_set)) {
      results_list = intervention_mediator_analysis(
        data = data_set,
        mediator = as.character(crossing_set[i, "mediator"]),
        mediator_type = as.character(crossing_set[i, "mediator_type"]))
      
      full_results_list[[i]] = results_list
    }
    
    combined_results_df <- do.call(rbind, full_results_list)
    combined_results_df = data.frame(lapply(combined_results_df, function(x) if(is.numeric(x)) round(x, 2) else x))
    View(combined_results_df)
    return(combined_results_df)
  }


#--------------------------------------------------------
# create crossing sets 
#--------------------------------------------------------

crossing_IM= rbind(crossing(
  mediator = binary_mediators,
  mediator_type = "binary"
), crossing(
  mediator = continuous_mediator,
  mediator_type = "continuous"
))

crossing_IM_main= rbind(crossing(
  mediator = binary_mediators_main,
  mediator_type = "binary"
), crossing(
  mediator = continuous_mediator,
  mediator_type = "continuous"
))

crossing_IM_malaria = crossing(
  mediator = "anymalaria",
  mediator_type = "binary"
)

#--------------------------------------------------------
# save the results
#--------------------------------------------------------

#prevent using scientific notations
options(scipen = 999)

IM_results_full = intervention_mediator_application(data_set = data_zscore_quarterly, crossing_set = crossing_IM)
View(IM_results_full)
saveRDS(IM_results_full, paste0(results_path,"intervention_mediator_full_results.RDS"))

IM_results_main = intervention_mediator_application(data_set = data_zscore_quarterly, crossing_set = crossing_IM_main)
View(IM_results_main)
saveRDS(IM_results_main, paste0(results_path,"IM-MO-stratified/intervention_mediator_main_results_stratified.RDS"))

