##############################################
##############################################
# Documentation: get_prop_CI_by_age
# Usage: get_prop_CI_by_age(data, strat_varname)
# Description: calculate two-tailed Wilson confidence intervals

# Args/Options:
# data: a data frame with variables strat_var, strat_level, and agecat
# strat_varname: a string for the stratification variable

# Returns: data frame with proportion, lower bound, upper bound, and labels
# Output: none


get_prop_CI_by_age = function(data, strat_varname, agecat_col){
  
  n = data[["n_prop"]] #I will not need this for continuous 
  N = data[["nmeas"]] #it is useful to know N 
  res = BinomCI(x = n, n = N) # i will not need this 
  res_out = data.frame(cbind(agevar = as.character(agecat_col),
                             as.data.frame(res)))

  rownames(res_out) = NULL
  
  res_out = res_out %>% mutate(strat_var = strat_varname,
                               strat_level = data$strat_level) %>% 
    rename(!!sym(agecat_col) := agevar)
  
  return(res_out)
}





##############################################
##############################################
# Documentation: get_prop_CI
# Usage: get_prop_CI(data, outcome_col, strat_varname)
# Description: calculate two-tailed Wilson confidence intervals
#              for a given outcome and stratificiation 

# Args/Options:
# data: a data frame with variables strat_var, strat_level, and agecat, outcome
# outcome_col: a string for the outcome name
# strat_varname: a string for the stratification variable

# Returns: data frame with proportion, lower bound, upper bound, and labels
# Output: none

library(DescTools)

# obtain proportion and standard error 
# for binary outcomes 
get_prop_CI <- function(data, outcome_col, strat_varname, agecat_col){
  
  # filter to observations for outcome of interest
  data_nomiss <- data %>%
    mutate(outcome=!!sym(outcome_col)) ##!!sym is so that haz_s_stunt gets recognized as a column name and not a string 
  
  # calculate proportion 
  # count measurements per study by age
    if(strat_varname==agecat_col){
      data_nomiss = data_nomiss %>% 
        group_by(!!sym(strat_varname)) 
    }else{
      data_nomiss = data_nomiss %>% 
        group_by(!!sym(agecat_col), !!sym(strat_varname)) #if the grouping var is not age we would group by both 
    }
  
  #for calculating props 
  summary.data = data_nomiss %>% 
    summarise(nmeas=sum(!is.na(outcome)), # total count
              n_prop=sum(outcome==1, na.rm=TRUE),  # count when outcome=1 
              prop=mean(outcome, na.rm=TRUE)) %>% 
    filter(!is.na(prop)) 
  
  if(strat_varname==agecat_col){
    summary.data = summary.data %>% 
      mutate(strat_level = "",
             strat_var = "")
  }else{
    summary.data = summary.data    %>% 
      rename(strat_level = !!sym(strat_varname)) %>% 
      mutate(strat_var = strat_varname)
  }

  # estimate standard error within each stratum 
  if (strat_varname==agecat_col) { 
    results = get_prop_CI_by_age(data=summary.data, 
                                 strat_var = strat_varname,
                                 agecat_col = agecat_col) %>% #will be replaced by get mean CI by age func
      mutate(strat_level = "")
  } else {
    results = get_prop_CI_by_age(data=summary.data, 
                                 strat_var = strat_varname,
                                 agecat_col = agecat_col)
  }
  
  results = results %>% 
    dplyr::select(!!sym(agecat_col), strat_var, strat_level, est, lwr.ci, upr.ci) %>% 
    rename(lb = lwr.ci,
           ub = upr.ci) %>% 
    mutate(strat_level = as.character(strat_level))
  
  # format results (I will not need this unless I have props)
  results = results %>%
    mutate(est=est*100,lb=lb*100,ub=ub*100)
  
  results$ptest.f=sprintf("%0.0f",results$est)
  
  return(results)
}



##############################################
##############################################
# Documentation: get_Inc_CI
# Usage: get_Inc_CI(data, outcome_col, strat_varname, agecat_col, age_level)
# Description: calculate two-tailed Wilson confidence intervals
#              for a given incidence outcome and stratificiation 

# Args/Options:
# data: a data frame with variables strat_var, strat_level, and agecat, outcome
# outcome_col: a string for the outcome name
# strat_varname: a string for the stratification variable
# agecat_col: a string for the age category variable
# age_level: a string for the age level of the age category variable 

# Returns: data frame with proportion, lower bound, upper bound, and labels
# Output: none

get_Inc_CI <- function(data, outcome_col, strat_varname, agecat_col, age_level){
  
  # filter to the outcome and age of interest
  data_age = data %>% filter(!!sym(agecat_col)==age_level) 
  keep_cols = colnames(data_age)[grep(outcome_col, colnames(data_age))]
  cases = keep_cols[grep("incident",keep_cols)]
  atrisk = keep_cols[grep("atrisk",keep_cols)]
  
  # drop rows not at risk 
  data_sub = data_age %>% dplyr::select(all_of(c("id", strat_varname, agecat_col, keep_cols))) %>% 
    filter(!!sym(atrisk)==T)

  # calculate proportion 
  # count measurements per study by age
  if(strat_varname==agecat_col){
    data_sub = data_sub %>% 
      group_by(!!sym(strat_varname)) 
  }else{
    data_sub = data_sub %>% 
      group_by(!!sym(agecat_col), !!sym(strat_varname)) #if the grouping var is not age we would group by both 
  }
  
  #for calculating props 
  summary.data = data_sub %>% 
    summarise(nmeas=sum(!is.na(!!sym(cases))), # total count
              n_prop=sum(!!sym(cases)==1, na.rm=TRUE),  # count when outcome=1 
              prop=mean(!!sym(cases), na.rm=TRUE)) %>% 
    filter(!is.na(prop)) 
  
  if(strat_varname==agecat_col){
    summary.data = summary.data %>% 
      mutate(strat_level = "",
             strat_var = "")
  }else{
    summary.data = summary.data    %>% 
      rename(strat_level = !!sym(strat_varname)) %>% 
      mutate(strat_var = strat_varname)
  }

  # create vector of stratification variable levels 
  func = paste0("data_sub$",strat_varname)
  list_of_levels = unique(eval(parse(text=func)))
  list_of_levels = as.list(list_of_levels[!is.na(list_of_levels)]) #square brackets will only !na values 
  
  # estimate standard error within each stratum 
  if (strat_varname==agecat_col) { 
    results = get_prop_CI_by_age(data=summary.data, 
                                 strat_var = strat_varname,
                                 agecat_col = agecat_col) %>% #will be replaced by get mean CI by age func
      mutate(strat_level = "")
  } else {
    results = get_prop_CI_by_age(data=summary.data, 
                                 strat_var = strat_varname,
                                 agecat_col = agecat_col)
  }
  
  results = results %>% 
    dplyr::select(!!sym(agecat_col), strat_var, strat_level, est, lwr.ci, upr.ci) %>% 
    rename(lb = lwr.ci,
           ub = upr.ci) %>% 
    mutate(strat_level = as.character(strat_level),
           age_level = age_level,
           outcome = outcome_col)
  
  # format results (I will not need this unless I have props)
  results = results %>%
    mutate(est=est*100,lb=lb*100,ub=ub*100)
  
  results$ptest.f=sprintf("%0.0f",results$est)
  
  results = dplyr::select(results, c(outcome, agecat_col, age_level, strat_var, strat_level, everything()))
  
  return(results)
}

