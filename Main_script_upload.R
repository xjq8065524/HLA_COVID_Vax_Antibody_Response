# Load library ------------------------------------------------------------
library(tidyverse)
library(qqman)
library(manhattanly)
library(ggforce)
library(patchwork)
# Universal functions -----------------------------------------------------
temp_func <- function( raw_model){
  
  temp <- summary( raw_model )
  
  p <- temp$coefficients %>% as.data.frame()%>% janitor::clean_names() %>% pull( pr_z)
  
  Cox_results <- 
    temp$conf.int %>% as_tibble( rownames = "Vars") %>% janitor::clean_names() %>%   
    mutate( coef = temp$coefficients[,"coef"],
            se = temp$coefficients[,3],
            n = temp$n,
            nevent = temp$nevent,
            HR_CI = paste(format(round(exp_coef, 2), nsmall = 2, scientific = FALSE), " (",
                          format(round(lower_95, 2), nsmall = 2), " to ",
                          format(round(upper_95, 2), nsmall = 2), ")",sep = ""),
            p_value = p) 
  return( Cox_results)
  
}
missingness_func <- function( df){
  
  df %>% summarise( across( everything(), ~ sum( is.na(.))/ n()))
  
  
}
tidy_allele_name_fun <- function( x){
  
  tidy_output <- 
    x %>% 
    as_tibble() %>% 
    separate_wider_delim( value , delim = "_", names = c( "gene_name", "allele_name"), cols_remove = FALSE) %>% 
    mutate( allele_name_first_two_char = case_when( str_count(allele_name) == 3 ~ str_c( "0", str_sub( allele_name, 1, 1)),
                                                    str_count(allele_name) == 4 ~ str_sub( allele_name, 1, 2),
                                                    TRUE ~ "999"),
            allele_name_last_two_char = case_when( str_count(allele_name) == 3 ~ str_sub( allele_name, 2, 3),
                                                   str_count(allele_name) == 4 ~ str_sub( allele_name, 3, 4),
                                                   TRUE ~ "999")) %>% 
    mutate( tidy_allele_name = str_c( gene_name, "*", allele_name_first_two_char, ":", allele_name_last_two_char)) %>% 
    pull( tidy_allele_name)
  
  return(tidy_output)
}

# Import baseline raw datasets ---------------------------------------------------------
# UKBB England participants
baseline_df <- readRDS("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/R_datasets/COVID-19_MAIN_df/baseline_df.rds")
# Import death registry ---------------------------------------------------
death <- read_delim("D:\\DPhil\\UK_Biobank_opioid_application\\Death\\Version_20230405\\death.txt", 
                    delim = "\t",
                    col_types = cols(date_of_death = col_date(format = "%d/%m/%Y")))

death_cause <- read_delim("D:\\DPhil\\UK_Biobank_opioid_application\\Death\\Version_20230405\\death_cause.txt", delim = "\t")

# non-covid-19 infection death
complete_death <- death %>% left_join( death_cause, by = c( "eid", "ins_index"))

# Import GP vaccination (diagnosis source) -----------------------------------------
vaccine_diagnosis_source <- readRDS("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/Sub_project/Breakthrough_infection/Derived_dataset/vaccine_diagnosis_source.rds")

# Import Covid-19 testing data --------------------------------------------
#England Covid-19 test results
covid19_result_england <- 
  data.table::fread( input = "D:/DPhil/UK_Biobank_opioid_application/COVID_test_results/Version_20230405/covid19_result_england.txt") %>% 
  mutate( specdate = lubridate::dmy( specdate))

# Import HLA data ---------------------------------------------------------------------
HLA_names <- 
  read.delim("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/Sub_project/Genetics_breakthrough_infection/New Text Document.txt", header=FALSE) %>% 
  pivot_longer( cols = everything(), names_to = "Position", values_to = "Name") %>% 
  separate_wider_delim( Name, delim = "_", names = c( "gene_name", "allel_name"), cols_remove = FALSE) 


# Import antibody test data -----------------------------------------------
ukb52351 <- read.delim("D:/DPhil/UK_Biobank_opioid_application/Add_Basket_4000278/ukb52351.tab")

# Curate cohort ------------------------------------------------
# combine all cohorts
Limit_first_dose_cohort_combined <- 
  Naive_first_dose_cohort %>% 
  select( eid, event_dt) %>% 
  left_join( First_dose_covariate_cohort, by = "eid") %>% 
  rename( index_date = event_dt) %>% 
  left_join( select( death, eid, date_of_death), by = "eid") %>% 
  left_join( covid_infection_cohort_england, by = "eid") %>% 
  left_join( covid_hospitalization_cohort_england, by = "eid") %>% 
  left_join( covid_19_death_cohort_england, by = "eid") %>% 
  left_join( Tidy_binary_HLA, by = c("eid" = "f.eid")) %>% 
  mutate( follow_up_end_date = case_when( is.na( date_of_death) ~ pmin( index_date + 7*12, na.rm = TRUE), 
                                          TRUE ~ pmin( date_of_death, index_date + 7*12, na.rm = TRUE)),
          across( .cols = starts_with("covid_"), ~ case_when( .x >= index_date & .x <= follow_up_end_date ~ 1, TRUE ~ 0), .names = "incident_outcome_{.col}"),
          across( .cols = starts_with("covid_"), ~ as.numeric( pmin( follow_up_end_date, .x, na.rm = TRUE) - index_date), .names = "follow_up_days_{.col}")) %>% 
  mutate( ethnicity = replace_na( ethnicity, "Do not know"), ethnicity_binary = case_when( ethnicity %in% c( "White") ~ "0", TRUE ~ "1")) 




# second-dose covid infection naive cohort
Naive_second_dose_cohort <- 
  Second_dose_cohort %>% 
  left_join( select(covid_infection_cohort_england, eid, covid_infection), by = "eid") %>% 
  left_join( select( death, eid, date_of_death), by = "eid") %>% 
  filter( is.na( covid_infection) | covid_infection > event_dt) %>% # no previous infection
  filter( is.na( date_of_death) | date_of_death >= event_dt) # survived

# Cohort study for BT-outcomes ---------------------------------
# run cox model
cox_model_infection_func <- function( genes, df){
  
  naive_adjust_formula <- 
    as.formula( paste( 
      paste( "survival::Surv( follow_up_days_covid_infection", ",", "incident_outcome_covid_infection", ")", "~", sep = ""), 
      paste( c( genes, "age", "sex", "batch", "ethnicity_binary", str_subset( names( baseline_df), "f_22009_")), collapse="+"),
      sep = ""))
  
  HR_df <- survival::coxph( formula = naive_adjust_formula , data = df) %>% temp_func() %>% filter( row_number() == 1)
  
  return( HR_df)
  
}
cox_model_hospitalization_func <- function( genes, df){
  
  naive_adjust_formula <- 
    as.formula( paste( 
      paste( "survival::Surv( follow_up_days_covid_hospitalization", ",", "incident_outcome_covid_hospitalization", ")", "~", sep = ""), 
      paste( c( genes, "age", "sex", "batch", "ethnicity_binary", str_subset( names( baseline_df), "f_22009_")), collapse="+"),
      sep = ""))
  
  HR_df <- survival::coxph( formula = naive_adjust_formula , data = df) %>% temp_func()%>% filter( row_number() == 1)
  
  return( HR_df)
  
}

# infection
survival::coxph( formula = survival::Surv(follow_up_days_covid_infection, incident_outcome_covid_infection) ~ DQB1_6 , data = First_dose_cohort_combined) %>% temp_func() %>% filter( row_number() == 1)

# hospitalization
survival::coxph( formula = survival::Surv(follow_up_days_covid_hospitalization, incident_outcome_covid_hospitalization) ~ DQB1_6 , data = First_dose_cohort_combined) %>% temp_func() %>% filter( row_number() == 1)

# KM curve
Overall_KM_infection_plot_func <- function( main_df, upper_y_limit){
  
  KM_model <- survminer::surv_fit( survival::Surv( follow_up_days_covid_infection, incident_outcome_covid_infection) ~ DQB1_6, data = First_dose_cohort_combined)
  
  KM_plot <- survminer::ggsurvplot( KM_model,
                                    conf.int = TRUE,
                                    risk.table = TRUE,
                                    cumevents = TRUE,
                                    ylim = c(0, 0.05),
                                    xlim = c(0, 272),
                                    #ncensor.plot = TRUE,
                                    fun = "event",
                                    xlab = "Days after the first vaccination",
                                    ylab = "Cumulative incidence",
                                    break.time.by = 30,
                                    censor.shape = 124,
                                    censor = FALSE,
                                    palette = c("#999999", "#FF6666"),
                                    legend.title = "",
                                    legend.labs = c( "Non HLA-DQB1*06 carrier", "HLA-DQB1*06 carrier"),
                                    tables.y.text.col = TRUE,
                                    tables.y.text = FALSE,
                                    fontsize = unit( 2.2, "cm"),
                                    tables.height = 0.2,
                                    size = 0.5,
                                    ggtheme = theme_classic())
  
  cumulative_incidence <- 
    KM_plot$plot$data %>% 
    group_by( DQB1_6) %>% 
    summarise( max( surv))
  
  print(cumulative_incidence)
  
  
  p1 = KM_plot$plot + theme( text = element_text( family = "serif", color = "black"), 
                             axis.text = element_text( family = "serif", color = "black"),
                             panel.grid.major.y = element_line( colour = "grey", linetype = "dotted"),
                             legend.position = c(0.15,0.8))
  
  # p2 = KM_plot$table + theme( text = element_text( family = "serif"),
  #                             axis.text = element_blank(),
  #                             plot.title = element_text( size = 8, face = "bold")) + 
  #   labs( x = "")
  # 
  # p3 = KM_plot$cumevents + theme( text = element_text( family = "serif"),
  #                                 axis.text = element_blank(),
  #                                 plot.title = element_text( size = 8, face = "bold")) + 
  #   labs( x = "")
  
  # combined <- cowplot::plot_grid( p1, p2, p3, align = "v", ncol =1,rel_heights = c(2.5,1,1))
  
  return( p1)
  
}
First_window_KM_infection_plot_func <- function( main_df, upper_y_limit){
  
  KM_model <- survminer::surv_fit( survival::Surv( follow_up_days_covid_infection, incident_outcome_covid_infection) ~ DQB1_6, data = Limit_first_dose_cohort_combined)
  
  KM_plot <- survminer::ggsurvplot( KM_model,
                                    conf.int = TRUE,
                                    risk.table = TRUE,
                                    cumevents = TRUE,
                                    ylim = c(0, 0.025),
                                    xlim = c(0, 85),
                                    #ncensor.plot = TRUE,
                                    fun = "event",
                                    xlab = "Days after the first vaccination",
                                    ylab = "Cumulative incidence",
                                    break.time.by = 7,
                                    censor.shape = 124,
                                    censor = FALSE,
                                    palette = c("#999999", "#FF6666"),
                                    legend.title = "",
                                    legend.labs = c( "Non HLA-DQB1*06 carrier", "HLA-DQB1*06 carrier"),
                                    tables.y.text.col = TRUE,
                                    tables.y.text = FALSE,
                                    fontsize = unit( 2.2, "cm"),
                                    tables.height = 0.2,
                                    size = 0.5,
                                    ggtheme = theme_classic())
  
  cumulative_incidence <- 
    KM_plot$plot$data %>% 
    group_by( DQB1_6) %>% 
    summarise( max( surv))
  
  print(cumulative_incidence)
  
  
  p1 = KM_plot$plot + theme( text = element_text( family = "serif", color = "black"), 
                             axis.text = element_text( family = "serif", color = "black"),
                             panel.grid.major.y = element_line( colour = "grey", linetype = "dotted"),
                             legend.position = "none")
  
  # p2 = KM_plot$table + theme( text = element_text( family = "serif"),
  #                             axis.text = element_blank(),
  #                             plot.title = element_text( size = 8, face = "bold")) + 
  #   labs( x = "")
  # 
  # p3 = KM_plot$cumevents + theme( text = element_text( family = "serif"),
  #                                 axis.text = element_blank(),
  #                                 plot.title = element_text( size = 8, face = "bold")) + 
  #   labs( x = "")
  
  # combined <- cowplot::plot_grid( p1, p2, p3, align = "v", ncol =1,rel_heights = c(2.5,1,1))
  
  return( p1)
  
}
Second_window_KM_infection_plot_func <- function( main_df, upper_y_limit){
  
  KM_model <- survminer::surv_fit( survival::Surv( follow_up_days_covid_infection, incident_outcome_covid_infection) ~ DQB1_6, data = Second_dose_cohort_combined)
  
  KM_plot <- survminer::ggsurvplot( KM_model,
                                    conf.int = TRUE,
                                    risk.table = TRUE,
                                    cumevents = TRUE,
                                    ylim = c(0, 0.025),
                                    xlim = c(0, 85),
                                    #ncensor.plot = TRUE,
                                    fun = "event",
                                    xlab = "Days after the second vaccination",
                                    ylab = "Cumulative incidence",
                                    break.time.by = 7,
                                    censor.shape = 124,
                                    censor = FALSE,
                                    palette = c("#999999", "#FF6666"),
                                    legend.title = "",
                                    legend.labs = c( "Non HLA-DQB1*06 carrier", "HLA-DQB1*06 carrier"),
                                    tables.y.text.col = TRUE,
                                    tables.y.text = FALSE,
                                    fontsize = unit( 2.2, "cm"),
                                    tables.height = 0.2,
                                    size = 0.5,
                                    ggtheme = theme_classic())
  
  cumulative_incidence <- 
    KM_plot$plot$data %>% 
    group_by( DQB1_6) %>% 
    summarise( max( surv))
  
  print(cumulative_incidence)
  
  
  p1 = KM_plot$plot + theme( text = element_text( family = "serif", color = "black"), 
                             axis.text = element_text( family = "serif", color = "black"),
                             panel.grid.major.y = element_line( colour = "grey", linetype = "dotted"),
                             legend.position = "none")
  
  # p2 = KM_plot$table + theme( text = element_text( family = "serif"),
  #                             axis.text = element_blank(),
  #                             plot.title = element_text( size = 8, face = "bold")) + 
  #   labs( x = "")
  # 
  # p3 = KM_plot$cumevents + theme( text = element_text( family = "serif"),
  #                                 axis.text = element_blank(),
  #                                 plot.title = element_text( size = 8, face = "bold")) + 
  #   labs( x = "")
  
  # combined <- cowplot::plot_grid( p1, p2, p3, align = "v", ncol =1,rel_heights = c(2.5,1,1))
  
  return( p1)
  
}


KM_hospitalization_plot_func <- function( main_df, upper_y_limit){
  
  KM_model <- survminer::surv_fit( survival::Surv( follow_up_days_covid_hospitalization, incident_outcome_covid_hospitalization) ~ HLA_DQB1_6_index, data = main_df)
  
  KM_plot <- survminer::ggsurvplot( KM_model,
                                    conf.int = TRUE,
                                    risk.table = TRUE,
                                    cumevents = TRUE,
                                    ylim = c(0, upper_y_limit),
                                    xlim = c(0, 272),
                                    #ncensor.plot = TRUE,
                                    fun = "event",
                                    xlab = "Days after the index date",
                                    ylab = "Cumulative incidence",
                                    break.time.by = 30,
                                    censor.shape = 124,
                                    censor = FALSE,
                                    palette = c("#374E55FF", "#DF8F44FF"),
                                    legend.title = "",
                                    legend.labs = c( "Non HLA-DQB1*06 carrier", "HLA-DQB1*06 carrier"),
                                    tables.y.text.col = TRUE,
                                    tables.y.text = FALSE,
                                    fontsize = unit( 2.2, "cm"),
                                    tables.height = 0.2,
                                    size = 0.5,
                                    ggtheme = theme_classic())
  
  cumulative_incidence <- 
    KM_plot$plot$data %>% 
    group_by( HLA_DQB1_6_index) %>% 
    summarise( max( surv))
  
  print(cumulative_incidence)
  
  
  p1 = KM_plot$plot + theme( text = element_text( family = "serif", color = "black"), 
                             axis.text = element_text( family = "serif", color = "black"),
                             panel.grid.major.y = element_line( colour = "grey", linetype = "dotted"))
  
  p2 = KM_plot$table + theme( text = element_text( family = "serif"),
                              axis.text = element_blank(),
                              plot.title = element_text( size = 8, face = "bold")) + 
    labs( x = "")
  
  p3 = KM_plot$cumevents + theme( text = element_text( family = "serif"),
                                  axis.text = element_blank(),
                                  plot.title = element_text( size = 8, face = "bold")) + 
    labs( x = "")
  
  combined <- cowplot::plot_grid( p1, p2, p3, align = "v", ncol =1,rel_heights = c(2.5,1,1))
  
  return( combined)
  
}

Overall_infection_KM <- KM_infection_plot_func( main_df = First_dose_cohort_combined, upper_y_limit = 0.08)
Overall_hospitalization_KM <- KM_hospitalization_plot_func( main_df = Overall_cohort_model_cohort, upper_y_limit = 0.008)

# HLA alleles ---------------------------------------------
LR_model_func <- function( genes = "DQB1_6", df){
  
  naive_adjust_formula <- 
    as.formula( paste( 
      paste( "antibody_test_result", "~", sep = ""), 
      paste( c( genes, "age", "sex", "ethnicity_binary", str_subset( names( baseline_df), "f_22009_")), collapse="+"),
      sep = ""))
  
  LR_list <- glm( naive_adjust_formula, data = df, family = binomial( link = "logit"))
  
  #exp(confint(LR_list))

  LR_df <- summary(LR_list)$coefficients[2,]
  
  return( LR_df)
  
}

# Discovery analysis
subtype_signal_after_first_dose <- 
  map_df( set_names(full_test_list$genes), LR_model_func, df = subtype_count_only_first_batch_2, .id = "genes") %>% janitor::clean_names() %>% 
  mutate( p_after_FDR = p.adjust( pr_z, method = "fdr", n = length(pr_z)))

significance_after_first_dose_list <- subtype_signal_after_first_dose %>% filter( p_after_FDR < 0.05) %>% pull( genes)

subtype_signal_after_second_dose <- 
  map_df( set_names(full_test_list$genes), LR_model_func, df = subtype_count_first_and_second_batch_2, .id = "genes") %>% janitor::clean_names() %>% 
  mutate( p_after_FDR = p.adjust( pr_z, method = "fdr", n = length(pr_z)))

significance_after_second_dose_list <- subtype_signal_after_second_dose %>% filter( p_after_FDR < 0.05) %>% pull( genes)

# Validation analysis

# validate HLA alleles identified from the 1st dose cohort
valid_subtype_signal_after_first_dose <-
  map_df( set_names(set_names(full_test_list$genes)), LR_model_func, df = subtype_count_only_first_batch_1, .id = "genes") %>% janitor::clean_names()


only_valid_subtype_signal_after_first_dose <- 
  map_df( set_names(set_names(significance_after_first_dose_list)), LR_model_func, df = subtype_count_only_first_batch_1, .id = "genes") %>% janitor::clean_names() %>% 
  mutate( p_after_FDR = p.adjust( pr_z, method = "fdr", n = length(pr_z)))

final_validated_list_after_first_dose <- 
  select( subtype_signal_after_first_dose, genes, estimate, pr_z, p_after_FDR) %>% 
  left_join( select( only_valid_subtype_signal_after_first_dose, genes, estimate, pr_z, p_after_FDR), by = "genes") %>% 
  filter( p_after_FDR.x < 0.05, p_after_FDR.y <0.05) %>% 
  arrange(desc(estimate.y))

final_validated_independent_list_after_first_dose <- 
  final_validated_list_after_first_dose %>% 
  filter( genes %in% all_final_independent_validated_list$genes)


# validate HLA alleles identified from the 2st dose cohort
valid_subtype_signal_after_second_dose <-
  map_df( set_names(set_names(full_test_list$genes)), LR_model_func, df = subtype_count_first_and_second_batch_1, .id = "genes") %>% janitor::clean_names()

only_valid_subtype_signal_after_second_dose <- 
  map_df( set_names(significance_after_second_dose_list), LR_model_func, df = subtype_count_first_and_second_batch_1, .id = "genes") %>% janitor::clean_names()  %>% 
  mutate( p_after_FDR = p.adjust( pr_z, method = "fdr", n = length(pr_z)))

final_validated_list_after_second_dose <- 
  select( subtype_signal_after_second_dose, genes, estimate, pr_z, p_after_FDR) %>% 
  left_join( select( only_valid_subtype_signal_after_second_dose, genes, estimate, pr_z, p_after_FDR), by = "genes") %>% 
  filter( p_after_FDR.x < 0.05, p_after_FDR.y <0.05)

all_final_validated_list <- unique( c(final_validated_list_after_first_dose$genes, final_validated_list_after_second_dose$genes))

all_final_validated_list %>% str_extract(".*(?=_)") %>% unique()

#overall
overall_subtype_signal_after_first_dose <- 
  map_df( set_names(full_test_list$genes), LR_model_func, df = subtype_count_only_first, .id = "genes") %>% janitor::clean_names() %>% 
  mutate( genes_order = fct_reorder( genes, estimate, .desc = TRUE))

overall_subtype_signal_after_second_dose <- 
  map_df( set_names(full_test_list$genes), LR_model_func, df = subtype_count_first_and_second, .id = "genes") %>% janitor::clean_names() %>% 
  mutate( genes_order = fct_reorder( genes, estimate, .desc = TRUE))


# HLA alleles ---------------------------------------------
# bonferroni significant alleles
BF_func <- function( df, filter_string, coef){
  
  significance <- 
    coef %>% 
    filter( genes %in% filter_string)
  
  # a data-driven PRS
  new_one_func <- function( x, cof){
    
    y = x*cof
    
    return(y)
    
  }
  
  Score_temp <- 
    map2_df( df[ ,significance$genes], significance$estimate, new_one_func) 
  
  Score_combined <- 
    Score_temp %>% 
    mutate( score = rowSums( Score_temp[ ,significance$genes], na.rm = TRUE)) %>% 
    bind_cols( select(df, !all_of(full_test_list$genes))) %>% 
    mutate( score_cat_5 = as.factor(as.character(ntile( score, 5)))) %>% 
    mutate( score_sd = scale(score))
  
  return( Score_combined)
  
}

All_one_dose_score <- BF_func( df = subtype_count_only_first, filter_string = overall_subtype_signal_after_first_dose$genes, coef = overall_subtype_signal_after_first_dose)
P_one_dose_score <- BF_func( df = subtype_count_only_first, filter_string = "pr_z <= 0.05", doses = subtype_signal_after_first_dose)
FDR_one_dose_score <- BF_func( df = subtype_count_only_first, filter_string = "p_after_FDR < 0.05", doses = subtype_signal_after_first_dose)

All_two_dose_score <- BF_func( df = subtype_first_and_second, filter_string = "pr_z >= 0", doses = subtype_signal_after_second_dose)
P_two_dose_score <- BF_func( df = subtype_first_and_second, filter_string = "pr_z <= 0.05", doses = subtype_signal_after_second_dose)
FDR_two_dose_score <- BF_func( df = subtype_first_and_second, filter_string = "p_after_FDR < 0.05", doses = subtype_signal_after_second_dose)


# clinically count alleles 
all_cohort_score_func <- function( df = subtype_count_only_first, 
                                   filter_string = final_validated_independent_list_after_first_dose$genes, 
                                   coef = final_validated_independent_list_after_first_dose){
  
  significance <- 
    coef %>% 
    filter( genes %in% filter_string)
  
  enhancer <- significance$genes[significance$estimate.x>0]
  suppressor <- significance$genes[significance$estimate.x<0]
  
  new_one_func <- function( x, cof){
    
    y = x*cof
    
    return(y)
    
  }
  
  Score_temp <- 
    map2_df( df[ ,as.character(significance$genes)], significance$estimate.x, new_one_func) %>% 
    mutate( across( .col = all_of(enhancer), ~ case_when(. > 0 ~ 1, TRUE ~ .)),
            across( .col = all_of(suppressor), ~ case_when(. == 0 ~ 1, TRUE ~ 0)))
    
  Score_combined <- 
    Score_temp %>% 
    mutate( score = rowSums( Score_temp[ ,as.character(significance$genes)], na.rm = FALSE)) %>% 
    bind_cols( select(df, -all_of(Number_of_allels$Name))) %>% 
    mutate( score_cat = as.character(score),
            score_cat_5 = as.factor(as.character(ntile( score, 5)))) 
  

  naive_adjust_formula <- 
    as.formula( paste( 
      paste( "antibody_test_result", "~", sep = ""), 
      paste( c( "score_cat", "age", "sex", "batch", "ethnicity_binary", str_subset( names( baseline_df), "f_22009_")), collapse="+"),
      sep = ""))
    
  LR_list <- glm( naive_adjust_formula, data = Score_combined, family = binomial( link = "logit"))
  summary(LR_list)  
  
  tidy_output <- broom::tidy(LR_list, conf.int = TRUE) %>% filter( str_detect( term, "score_cat"))
    
  return(tidy_output)

}


# individual alleles effect on COVID-19 outcome ------------------------------------
# run cox model
cox_model_infection_func <- function( genes, df){
  
  naive_adjust_formula <- 
    as.formula( paste( 
      paste( "survival::Surv( follow_up_days_covid_infection", ",", "incident_outcome_covid_infection", ")", "~", sep = ""), 
      paste( c( genes, "age", "sex", "batch", "ethnicity_binary", str_subset( names( baseline_df), "f_22009_")), collapse="+"),
      sep = ""))
  
  HR_df <- survival::coxph( formula = naive_adjust_formula , data = df) %>% temp_func() %>% filter( row_number() == 1)
  
  return( HR_df)
  
}
cox_model_severity_func <- function( genes, df){
  
  naive_adjust_formula <- 
    as.formula( paste( 
      paste( "survival::Surv( follow_up_days_covid_severity", ",", "incident_outcome_covid_severity", ")", "~", sep = ""), 
      paste( c( genes, "age", "sex", "batch", "ethnicity_binary", str_subset( names( baseline_df), "f_22009_")), collapse="+"),
      sep = ""))
  
  HR_df <- survival::coxph( formula = naive_adjust_formula , data = df) %>% temp_func()%>% filter( row_number() == 1)
  
  return( HR_df)
  
}

All_allele_infection <- 
  map_df( full_test_list$genes, cox_model_infection_func, df = First_dose_cohort_combined) %>% 
  left_join( select( overall_subtype_signal_after_first_dose, genes, estimate, std_error, pr_z), by = c("vars" = "genes")) %>% 
  left_join( select( overall_subtype_signal_after_second_dose, genes, estimate, std_error, pr_z), by = c("vars" = "genes")) %>% 
  left_join( full_test_list, by = c("vars" = "genes"))

All_allele_severity <- 
  map_df( full_test_list$genes, cox_model_severity_func, df = First_dose_cohort_combined) %>% 
  left_join( select( overall_subtype_signal_after_first_dose, genes, estimate, std_error, pr_z), by = c("vars" = "genes")) %>% 
  left_join( select( overall_subtype_signal_after_second_dose, genes, estimate, std_error, pr_z), by = c("vars" = "genes")) %>% 
  left_join( full_test_list, by = c("vars" = "genes"))

# Combine alleles effect on COVID-19 outcome ---------------------------------------------------------------
# Cox model (continuous)
One_overall_cohort_model_cohort <- 
  First_dose_cohort_combined  %>% 
  left_join( select(One_dose_allele_score_cat  , eid, score, score_neg, score_sd, score_cat, contains("score_manual")), by = c("eid")) %>% 
  mutate( age = age/10, age_sd = scale(age),sex = sex)

Two_overall_cohort_model_cohort <- 
  First_dose_cohort_combined  %>% 
  left_join( select(Second_dose_allele_score_cat  , eid, score, score_neg, score_sd, score_cat, contains("score_manual")), by = c("eid")) %>% 
  mutate( age = age/10, age_sd = scale(age), sex = sex)


One_dose_severity <- 
  cox_model_severity_func( genes = "score", df = One_overall_cohort_model_cohort) %>% 
  bind_rows( cox_model_severity_func( genes = "score_manual_1", df = One_overall_cohort_model_cohort)) %>% 
  bind_rows( cox_model_severity_func( genes = "score_manual_2", df = One_overall_cohort_model_cohort)) %>% 
  bind_rows( cox_model_severity_func( genes = "score_manual_3", df = One_overall_cohort_model_cohort)) %>% 
  filter( str_detect( vars, "score")) 

Two_dose_infection <- 
  cox_model_infection_func( genes = "score", df = Two_overall_cohort_model_cohort) %>% 
  bind_rows( cox_model_infection_func( genes = "score_manual_1", df = Two_overall_cohort_model_cohort)) %>% 
  bind_rows( cox_model_infection_func( genes = "score_manual_2", df = Two_overall_cohort_model_cohort)) %>% 
  bind_rows( cox_model_infection_func( genes = "score_manual_3", df = Two_overall_cohort_model_cohort)) %>% 
  filter( str_detect( vars, "score")) 

Two_dose_severity <- 
  cox_model_severity_func( genes = "score", df = Two_overall_cohort_model_cohort) %>% 
  bind_rows( cox_model_severity_func( genes = "score_manual_1", df = Two_overall_cohort_model_cohort)) %>% 
  bind_rows( cox_model_severity_func( genes = "score_manual_2", df = Two_overall_cohort_model_cohort)) %>% 
  bind_rows( cox_model_severity_func( genes = "score_manual_3", df = Two_overall_cohort_model_cohort)) %>% 
  filter( str_detect( vars, "score")) 

