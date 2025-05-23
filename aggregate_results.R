##########################################################################
# After augsynth has been run for all cohorts, 

#### Run Analyses ####
## Install "remotes" package to download augsynth from GitHub
if (!("remotes" %in% installed.packages()))
  install.packages('remotes')
## Load augsynth
library(augsynth)
source("xaugsynth.r")
tx_states <- c("FL", "MD", "MN", "NH", "NY", "OK", "PA")
ctrl_states <- c("AL", "GA",  "IA",  "ID",  "IN",  "KS",  "KY",  "MS",  "NC",  "NE", "SC",  "SD",  "TN",  "TX",  "VA",  "WI",  "WY")
cohort.name <- "outcome"
outcomes <- outcome.vars
outcomes <- paste0(outcomes, "_proportion")

get_augsynth_avg_results <- function(){
  
 
  ## Create vectors of outcomes of interest and covariates to use in outcome model
  ## for augsynth
  ## Loop over all outcomes of interest (contained in a vector called "outcomes")
  ## and all treated states.
  for (outcome in outcomes){
    
    # Add _proportion to covariate names and add mean_age variable to the list
    covariates  <- paste0(demo.vars, "_proportion")
    covariates <- append(covariates, "mean_age")
    outcome.name <- gsub("_proportion", "", outcome)
    # Loop over all states in a cohort and for each outcome
    for (tx_state in state.list) {
      # Create a unique analysis name for each cohort, outcome and state
      analysis_name <- paste(cohort.name, outcome, tx_state, sep = "_")
      # Read augsynth results from the state-outcome working folder
      augsynth_results[[analysis_name]] <<- readRDS(paste0(base.path, "state_cohorts/", cohort.name, "/", tx_state,"/", outcome.name, "/augsynth_results_adj.xpt"))
      # Read the file with weights assigned from the augsynth run
      wt.file <- paste0(base.path, "state_cohorts/", cohort.name, "/", tx_state,"/", outcome.name, "/augsynth_weights_adj.csv")
      # Read weights from input weights file
      weights <- read.csv(wt.file)
      weights <- weights[, -1]
      colnames(weights)[colnames(weights) == "STATE_CODE"] <- "state"
      
      # Read the clean dataset used to run augsynth for each model
      data.file <- paste0(base.path, "state_cohorts/", cohort.name, "/", tx_state,"/", outcome.name, "/", cohort.name, "_", tx_state, "_clean.xpt")
      # Read the file and calculate the standard deviation for covariates
      d <- read_xpt(data.file)
      d.sd <- d %>%
        group_by(state) %>%
        summarize(across(
          all_of(covariates),
          ~ sd(.x, na.rm = TRUE),
          .names = "{.col}_sd"
        ))

      # Merge the input data for augsynth model with weights from augsynth model, and the standard deviations
      d <- merge(d, weights, by= "state", all= T)
      d <- merge(d, d.sd, by="state", all=T)
      
      # Give treated states a weight of 1 so multiplication with weights works
      d$weight[is.na(d$weight)] <- 1
      
      # Num. individuals with non-zero weights
      d$Nprime <- d$enrolled_bene * !(abs(d$weight) < 1e-16)
      # Rename enrolled bene to STATE_N
      d$STATE_N <- d$enrolled_bene
      # Create TXSTATE, POST and TIME_SINCE_TX columns
      d <- d %>%
        mutate(TXSTATE = ifelse(state==cohort_state, 1, 0)) %>%
        mutate(treat_year_month = as.yearmon(as.character(treat_year_month), format = "%Y%m")) %>%
        mutate(POST = ifelse(year_month_num>=treat_year_month, 1, 0)) %>%
        mutate(year_month_num = as.yearmon(year_month_num, format="%Y%m")) %>%
        # Calculate time difference in months
        mutate(TIME_SINCE_TX = as.integer(round(12 * (year_month_num - treat_year_month))))
      
      # Add new columns to d representing the sum of individual covariate values
      # per state-time (i.e., xbar * n) and augsynth-weighted means (note the
      # latter are means b/c weights are pre-normalized)
      # Dynamically generate SD column names from covariates
      sd_columns <- paste0(covariates, "_sd")
      d <- cbind(d,
                 # multiply mean covariates by N to get sum of indiv. values
                 d[, covariates] * d$STATE_N,
                 # multiply mean covariates by augsynth weights
                 d[, covariates] * d$weight,
                 # multiply covariate SDs(^2) by (N-1) to get SSE
                 d[, sd_columns]^2 * (d$STATE_N - 1),
                 d[, sd_columns]^2 * (d$STATE_N - 1) * d$weight
      ) |>
        structure(names = c(names(d),
                            paste0(covariates, "_n"),
                            paste0(covariates, "_wt"),
                            paste0(covariates, "_var_n"),
                            paste0(covariates, "_var_wt")))
      
      # Now aggregate by time and state: sum all of the N- and weight-multiplied
      # values over states. Gives, for each covariate
      #\sum_{time} \sum_{state} n_{state,time}\bar{X}_{state,time},
      #\sum_{time} \sum_{state} w_{state}\bar{X}_{state,time},
      #\sum_{time} \sum_{state} (n_{state,time} - 1) s^2_{state,time}
      agg <- aggregate(
        x = d[, c("STATE_N", "weight", "Nprime",
                  grep("_(n|wt)$", names(d), value = T))],
        by = list("POST" = d$POST, "TXSTATE" = d$TXSTATE),
        sum
      )
      
      # Divide overall covariate sums by either total num. observations or the
      # total weights across time
      overall_means <- cbind(
        agg[, 1:2],
        agg[, paste0(covariates, "_n")]  / agg$STATE_N,
        agg[, paste0(covariates, "_wt")] / agg$weight
      )
      
      # Calculate MEAN SSEs
      meanSSEs <- do.call(rbind, lapply(1:4, \(i) { # Loops over 4 (post-txstate) combos i.e. 00, 01, 10, 11 
        x <- subset(d, POST == overall_means$POST[i] &
                      TXSTATE == overall_means$TXSTATE[i])
        m <- subset(overall_means, POST == overall_means$POST[i] &
                      TXSTATE == overall_means$TXSTATE[i]) # subset data and aggregated means
        
        # Create variance calculations dynamically for each covariate
        var_calcs <- lapply(covariates, function(covar) {
          c(
            # Unweighted variance
            sum(x$STATE_N * (x[[covar]] - m[[paste0(covar, "_n")]])^2),
            # Weighted variance
            sum(x$weight * x$STATE_N * (x[[covar]] - m[[paste0(covar, "_wt")]])^2)
          )
        }) %>% unlist()
        
        # Create named vector for data.frame
        c(
          setNames(var_calcs, paste0(rep(covariates, each = 2),
                                     c("_var_n", "_var_wt"))),
          "STATE_N" = sum(x$STATE_N),
          "weight" = sum(x$weight * x$STATE_N),
          "Nprime" = sum(x$Nprime)
        )
      }))
      meanSSEs <- as.data.frame(meanSSEs)
      
      # Dynamic column identifiers
      id_cols <- names(agg)[1:2]
      var_cols <- grep("_var", names(agg), value = TRUE)
      unweighted_vars <- grep("_var_n$", var_cols, value = TRUE)
      weighted_vars <- grep("_var_wt$", var_cols, value = TRUE)
      
      # Overall vars list from the column identifiers
      overall_vars <- agg[, c(id_cols, var_cols)]
      
      # getting adjusted unweighted variance (variance + SSE unweighted / (state_N-1)) 
      overall_vars[, unweighted_vars] <- 
        (overall_vars[, unweighted_vars] + meanSSEs[, grep("_var_n$", names(meanSSEs))]) / 
        (meanSSEs$STATE_N - 1)
      
      # adjust weighted variances (variance+sse weighted / ((nprime-1) * weight))*nprime
      overall_vars[, weighted_vars] <- 
        (overall_vars[, weighted_vars] + meanSSEs[, grep("_var_wt$", names(meanSSEs))]) * 
        meanSSEs$Nprime / 
        (meanSSEs$weight * (meanSSEs$Nprime - 1))
      
      # Get overall standard deviations from the above variance table
      overall_sds <- cbind(
        overall_vars[, id_cols, drop = FALSE],  # Preserve dataframe structure
        sqrt(overall_vars %>% select(all_of(var_cols)))     # Apply sqrt to all variance columns dynamically
      ) %>% 
        rename_at(vars(all_of(var_cols)), ~ gsub("_var", "_sd", .))  # Rename columns appropriately
      
      # PW: Store overall means and SDs in lists
      cohort_means[[analysis_name]] <<- overall_means
      cohort_sds[[analysis_name]]   <<- overall_sds
      
      # Save the inout dataset for augsynth for each state
      d <- d[order(d$state, d$year_month_num), ]
      assign(tx_state, d, envir = .GlobalEnv)
    }
  }
  
  # All of the augsynth results will get stored in a huge list. Now we can apply
  # the averaging for the ATT stuff:
  # Store augsynth results as outcome headers and outcome_State
  # Then using grep() applies averageOutcome to each outcome-state result
  augsynth_average_results <- structure({
    lapply(outcomes, \(outcome) {
      results_for_outcome <- grep(outcome, names(augsynth_results))
      structure({
        lapply(results_for_outcome, function(i) {
          message(names(augsynth_results)[i]) # This prints which model it is evaluating
          averageOutcomes(augsynth_results[[i]]) # This function returns average results  
        })
      }, names = names(augsynth_results)[results_for_outcome] # Name of each structure
      )
    })
  }, names = outcomes)
  #
  #
  
  assign("augsynth_average_results", augsynth_average_results, envir = .GlobalEnv)
}
#############################################################################
# This function gets pre treatment and post treatment results  
# Results include Mean, difference between pre/post mean, estimated ATT, ATT std err,.
# 95% lower and upper bounds and p value
get_pre_post_augsynth_results <- function(){
  # Get, for every augsynth, the overall ATT and its standard error
  overall_average_results <- structure({
    lapply(augsynth_average_results, \(y) {
      sapply(y, function(x) {
        cbind(x["ATT", "Estimate"], x["ATT", "Std. Err"])
      })
    })
  },
  names = outcomes)
  
  # Get an inverse variance weighted average ATT
  overall_ATT <- sapply(overall_average_results, \(z) {
    # Convert vectors to 2-row matrices
    if (!is.matrix(z)) {
      z <- matrix(z, nrow = 2)
    }
    # Calculate weighted average
    weighted.mean(x = z[1, ], w = 1/z[2, ]^2)
  })
  
  # Get standard errors
  overall_SE <- sapply(overall_average_results, \(z) {
    # Convert vectors to 2-row matrices
    if (!is.matrix(z)) {
      z <- matrix(z, nrow = 2)
    }
    # Calculate sd 
    sqrt(1/sum(1/z[2, ]^2))
  })
  #-------------------------------------------------------------------
  # Aggregate across all cohort means for outcomes
  # Get weighted/unweighted averages for POST and TREATMENT combos
  cohort_means_agg <- structure({
    # Loop through each outcome in the outcomes vector
    lapply(outcomes, \(outcome) {
      # Find all indices in cohort_means that match the current outcome name
      outcomeIndices <- grep(outcome, names(cohort_means))
      # Combine columns using inverse-variance weighting:
      cbind(
        # Keep first 2 columns unchanged (POST and TREAT)
        cohort_means[[outcomeIndices[1]]][, 1:2],
        # Calculate weighted average for remaining columns:
        Reduce("+", lapply(outcomeIndices, \(i) {
          # Align with overall_average_results structure
          o_av_res_index <- i - min(outcomeIndices) + 1
          x <- cohort_means[[i]]
          # Align with overall_average_results structure
          x[, 3:ncol(x)] * 1 / overall_average_results[[outcome]][2, o_av_res_index]^2
        })
        # Normalize weights
        ) / sum(1 / overall_average_results[[outcome]][2, ]^2)
      )
    })
  }, names = outcomes)
  # Aggregate across all columns bu getting average of all
  # Rename as cohort means for "AVERAGE OVER ALL OUTCOMES"
  cohort_means_agg <- c(cohort_means_agg,
                        "AVERAGE OVER ALL OUTCOMES" = 
                          list(Reduce("+", cohort_means_agg) / length(outcomes)))
  
  # Get pre period average outcomes
  prePeriod_avg_outcomes <- lapply(augsynth_average_results, \(x) {
    # Get(select) state-specific pre-tx outcome means & SDs
    ## treated state
    m_tx <- sapply(x, \(y) y[1, 1])
    s_tx <- sapply(x, \(y) y[1, 2])
    ## synthetic control
    m_ctrl <- sapply(x, \(y) y[3, 1])
    s_ctrl <- sapply(x, \(y) y[3, 2])
    
    # ATT SE
    w <- sapply(x, \(y) y[5, 2])
    
    # weight means & SDs by inverse ATT variances
    m_tx_w <- weighted.mean(m_tx, w = 1/w^2)
    s_tx_w <- weighted.mean(s_tx, w = 1/w^2)
    m_ctrl_w <- weighted.mean(m_ctrl, w = 1/w^2)
    s_ctrl_w <- weighted.mean(s_ctrl, w = 1/w^2)
    data.frame("tx" = c(T, F),
               "mean" = c(m_tx_w, m_ctrl_w),
               "sd" = c(s_tx_w, s_ctrl_w))
  })
  # Make pre period outcomes global for latter use in table 1
  assign("prePeriod_avg_outcomes", prePeriod_avg_outcomes, envir = .GlobalEnv)
  
  # Get post period average outcomes
  postPeriod_avg_outcomes <- lapply(augsynth_average_results, \(x) {
    # Get state-specific pre-tx outcome means & SDs
    ## treated state
    m_tx <- sapply(x, \(y) y[2, 1])
    s_tx <- sapply(x, \(y) y[2, 2])
    ## synthetic control
    m_ctrl <- sapply(x, \(y) y[4, 1])
    s_ctrl <- sapply(x, \(y) y[4, 2])
    
    # ATT SE
    w <- sapply(x, \(y) y[5, 2])
    
    # weight means & SDs by inverse ATT variances
    m_tx_w <- weighted.mean(m_tx, w = 1/w^2)
    s_tx_w <- weighted.mean(s_tx, w = 1/w^2)
    m_ctrl_w <- weighted.mean(m_ctrl, w = 1/w^2)
    s_ctrl_w <- weighted.mean(s_ctrl, w = 1/w^2)
    data.frame("tx" = c(T, F),
               "mean" = c(m_tx_w, m_ctrl_w),
               "sd" = c(s_tx_w, s_ctrl_w))
  })
  
  #-------------------------------------------------------------------
  # Combine pre and post period outcomes into a dataframe
  out <- data.frame(
    # interleave outcome name and blanks
    "Outcome" = rep(outcomes, each = 2),
    "Group" = rep(c("Policy States", "Synthetic Control"), length(outcomes)),
    "prePolicyMean" = unlist(sapply(prePeriod_avg_outcomes, \(x) x$mean, simplify = F)),
    "postPolicyMean" = unlist(sapply(postPeriod_avg_outcomes, \(x) x$mean, simplify = F))
  ) |>
    transform(
      "Difference" = postPolicyMean - prePolicyMean
    )
  
  # Create a dataframe to store ATT (Average Treatment Effect) results
  outATT <- data.frame(
    "Outcome" = outcomes, # Outcome variable names 
    "EstimatedATT" = overall_ATT, # Point estimates of treatment effects
    "ATT_SE" = overall_SE) # Standard errors of ATT estimates
  
  # Calculate 95% confidence intervals using normal approximation
  outATT <- transform(outATT,
                      "ATT 95% LB" = EstimatedATT - 1.96 * ATT_SE, # Lower bound: Mean - 1.96*SE
                      "ATT 95% UB" = EstimatedATT + 1.96 * ATT_SE) # Upper bound: Mean + 1.96*SE
  # Compute two-tailed p-values using z-test
  outATT$p_value <- 2 * (1 - pnorm(abs(outATT$EstimatedATT / outATT$ATT_SE)))
  
  # Merge with existing results dataframe containing pre/post means
  out <- merge(out, outATT, all = T)
  
  # Format out to print at minimum 4 significant digits
  out[, 3:9] <- lapply(out[, 3:9], formatC, digits = 4, format = "fg")
  
  # Optionally, format columns 3 to 9 (including p_value) to 4 significant digits
  cols_to_format <- 3:9
  out[, cols_to_format] <- lapply(out[, cols_to_format], function(x) formatC(as.numeric(x), digits = 4, format = "fg"))
  
  # If p_value is not in columns 3:9, adjust accordingly
  out$p_value <- formatC(out$p_value, digits = 4, format = "fg")
  
  
  # Write the overall ATT pre/post values as a csv table
  write.csv(out, paste0(base.path, "final_exhibits/overall_ATT.csv"), row.names=FALSE)
  assign("overall_ATT", out, envir = .GlobalEnv)
}
#############################################################################
# This function generates values for Table 1 of the main exhibit 
# Table 1 includes Mean, SD and Standardized Mean Difference for all covariates and outcomes
get_table1 <- function(){
  
  # Get combined means across all states
  # Combine all cohort_means into one data frame with "analysis" identifier
  combined_means <- bind_rows(cohort_means, .id = "analysis")
  
  average_means <- combined_means %>%
    group_by(POST, TXSTATE) %>% # For each POST and TREATMENT STATE combo
    summarise(across(
      # Select columns ending with "_n" or "_wt" and calculate mean
      ends_with(c("_n", "_wt")),
      ~ mean(., na.rm = TRUE)
    ), .groups = "drop") %>%
    # Reorder columns to match original structure
    select(POST, TXSTATE, everything())
  # Add the combined means as average_across_states in the cohort means list
  cohort_means[["average_across_states"]] <- average_means
  #-------------------------------------------------------------------
  
  # Get combined SDs across all states
  # Combine all cohort_sds into one data frame
  combined_sds <- bind_rows(cohort_sds, .id = "analysis")
  
  # Calculate average SDs across states
  average_sds <- combined_sds %>%
    group_by(POST, TXSTATE) %>% # For each POST and TREATMENT STATE combo
    summarise(across(
      # Select columns ending with "_sd_n" or "_sd_wt" (adjust based on your SD column names)
      # Calculate mean of the sd columns
      ends_with(c("_sd_n", "_sd_wt")),
      ~ mean(., na.rm = TRUE)
    ), .groups = "drop") %>%
    # Reorder columns to match original structure
    select( POST, TXSTATE, everything())
  
  # Add to cohort_sds list
  cohort_sds[["average_across_states"]] <- average_sds
  #-------------------------------------------------------------------
  # Write function that reshapes the mean/sd table to wide
  reshape_cohort_table <- function(wide_data) {
    m <-  wide_data %>%
      # Keep only pre-treatment period (POST=0)
      filter(POST == 0) %>%
      # Remove POST column
      select(-POST) %>%
      # Reshape to long format
      pivot_longer(
        cols = -TXSTATE,
        names_to = c("variable", "metric"),
        names_pattern = "(.*)_(n|wt)$",  # Split at last underscore before n/wt
        values_to = "value"
      ) %>%
      # Convert metric labels and group labels
      mutate(
        metric = case_when(
          metric == "n" ~ "unweighted",
          metric == "wt" ~ "weighted",
          TRUE ~ metric
        ),
        group = ifelse(TXSTATE == 1, "treated", "synth")
      ) %>%
      # Remove original TXSTATE column
      select(-TXSTATE) %>%
      # Reshape to wide format
      pivot_wider(
        names_from = c(metric, group),
        names_glue = "{metric}_{group}",
        values_from = value
      ) %>%
      # Reorder columns
      select(variable,
             unweighted_synth, unweighted_treated,
             weighted_synth, weighted_treated)
  }
  #-------------------------------------------------------------------
  # Select the average across means and sds evaluated above
  initMeans <- cohort_means[["average_across_states"]]
  initSDs <- cohort_sds[["average_across_states"]]
  # For sds, rename column by adding _sd at the end
  colnames(initSDs) <- gsub("_sd", "", colnames(initSDs))
  
  # Reshape both tables to wide
  m <- reshape_cohort_table(initMeans)
  s <- reshape_cohort_table(initSDs)
  
  # Merge mean and sd tables 
  x <- merge(m, s, by = "variable", suffixes = c("_mean", "_sd"))
  # Add standardized mean difference as the difference of the means divided by control group's sd
  x <- transform(x,
                 "unweighted_smd" =
                   abs(unweighted_treated_mean - unweighted_synth_mean) /
                   unweighted_synth_sd,
                 "weighted_smd" =
                   abs(weighted_treated_mean - weighted_synth_mean) /
                   weighted_synth_sd
  )
  # For each covariate's smd, format to 6 decimal places
  # Note: `flag = "#"` forces printing of trailing zeros (e.g., 0.000000 instead of 0)
  x[, grep("smd", names(x))] <- apply(x[, grep("smd", names(x))], 2,
                                      formatC, format = "f", flag = "#",
                                      digits = 6)
  # Replace exact zero SMD values with "<0.01" for readability
  x$weighted_smd[x$weighted_smd == "0.00"] <- "<0.01"
  x$unweighted_smd[x$unweighted_smd == "0.00"] <- "<0.01"
  
  # Ensure demographics appear first in tables for reporting standards
  # Define key demographic covariates of interest
  covs <- c( "female", "race_white", "race_black", "race_hispanic", "race_other", "mental_illness")
  covs <- paste0(covs, "_proportion") # Append "_proportion" to match column names
  covs <- append(covs, "mean_age") # Add continuous age variable
  # Reorder rows to prioritize demographic covariates
  x <- bind_rows(
    x %>% filter(variable %in% covs) %>% slice(match(covs, variable)), # 1. Demographics first, in specified order
    x %>% filter(!variable %in% covs) # 2. Other variables follow
  )
  # Select and reorder columns for final presentation
  x <- x[,c(1, 4, 8, 5, 9, 11)]
  
  # Define Table 1 descriptions for all variables
  table1.desc <- c(
    "female_proportion"="% Female",
    "race_white_proportion"="% Race - White",
    "race_black_proportion"="% Race - Black",
    "race_hispanic_proportion"="% Race - Hispanic",
    "race_other_proportion"="% Race - Other",
    "mental_illness_proportion"="% w/mental illness diagnosis",
    "mean_age"="Mean Age",
    "in_CUD_proportion"="% of patients with any inpatient CUD treatment", 
    "ot_CUD_proportion"="% of patients with any outpatient CUD treatment",  
    "CUD_initiation_proportion"="% of patients newly initiating CUD treatment",
    "in_OUD_proportion"="% of patients with any inpatient OUD treatment", 
    "ot_OUD_proportion"="% of patients with any outpatient OUD treatment",  
    "OUD_initiation_proportion"="% of patients newly initiating OUD treatment", 
    "OUD_med_proportion"="% of patients with any medication for OUD (MOUD)",
    "ovd_overall_cannabis_proportion"="% of patients with any CUD overdose utilization", 
    "ovd_overall_opioid_proportion"="% of patients with any OUD overdose utilization"
  )
  # Add description for all variables
  x$Var <- table1.desc[x$variable]
  
  #-------------------------------------------------------------------
  # Grab mean sample size to include in Table 1 as average benes across time period
  inp.dta <- read_sas(paste0(in.path, "/twelve_m_outcome_all_state.sas7bdat"))
  inp.dta$treated_state <- ifelse(inp.dta$state %in% tx_states, 1, 0)
  mean.ss <- inp.dta %>%
    group_by(treated_state) %>%  # replace 'group' with your actual column for treatment/control
    summarise(mean_enr_bene = mean(enrolled_bene, na.rm = TRUE))
  #-------------------------------------------------------------------
  # Table 1 formatting
  x_fin <- x %>%
    mutate(
      # FOrmat to 6 decimal places and show as Mean(SD), and rename the variables 
      synth_mean = sprintf("%.6f (%.6f)", weighted_synth_mean, weighted_synth_sd),
      treated_mean = sprintf("%.6f (%.6f)", weighted_treated_mean, weighted_treated_sd),
      smd = weighted_smd
    ) %>%
    # Now drop the original mean/sd columns
    select(-weighted_synth_mean, -weighted_synth_sd, -weighted_treated_mean, -weighted_treated_sd, -weighted_smd, -variable)
  #-------------------------------------------------------------------
  # Now that we got Table 1 for covariates,
  # Attach table 1 values for outcomes from Pre period avg outcome values
  prePeriod_avg_outcomes <- Map(function(df, nm) {
    df$outcome <- nm
    df
  }, prePeriod_avg_outcomes, names(prePeriod_avg_outcomes))
  
  # Combine all dataframes for each state and pivot to wide format
  prePeriod_avg_wide <- prePeriod_avg_outcomes %>%
    bind_rows() %>%
    pivot_wider(
      id_cols = outcome,  # Group by outcome
      names_from = tx,
      values_from = c(mean, sd),
      names_glue = "{ifelse(tx, 'treated', 'control')}_{.value}"
    ) %>%
    # Ensure no list columns remain
    unnest(everything())
  # Calculate SMD and format to 6 decimal places
  prePeriod_avg_wide$smd <- abs(prePeriod_avg_wide$treated_mean - prePeriod_avg_wide$control_mean) /prePeriod_avg_wide$control_sd
  prePeriod_avg_wide$smd <- sprintf("%.6f", as.numeric(prePeriod_avg_wide$smd))
  
  prePeriod_avg_wide <- prePeriod_avg_wide %>%
    # Format numeric columns to 6 decimal places as numeric (for further calculations)
    mutate(across(where(is.numeric), ~ round(.x, 6))) %>%
    # Create columns with mean(sd) format
    mutate(
      synth_mean = sprintf("%.6f (%.6f)", control_mean, control_sd),
      treated_mean = sprintf("%.6f (%.6f)", treated_mean, treated_sd)
    ) %>%
    rename(Var = outcome) %>%
    # Now drop the original mean/sd columns
    select(Var, synth_mean, treated_mean, smd)
  #-------------------------------------------------------------------
  # Attach description to outcomes too
  prePeriod_avg_wide$Var <- table1.desc[prePeriod_avg_wide$Var]
  # Combine with the covariate mean and sds from above and store as table 1
  table1 <- bind_rows(x_fin, prePeriod_avg_wide)
  return (table1)
}
#############################################################################
# Get Appendix plots for all outcomes
# That includes figure with error bar by state and att values and p values for each state
get_appx_exhibits <- function(){
  # Table for Appendix B
  for (outcome in outcomes) {
    
    fin.out.path <- paste0(base.path, "final_exhibits/", outcome)
    dir.create(fin.out.path)
    # Delete all files in the folder
    unlink(file.path(fin.out.path, "*"), recursive = FALSE, force = TRUE)
    
    # Extract state-specific ATT results and format
    stateSpecificATTs <-
      do.call(rbind,
              lapply(1:length(augsynth_average_results[[outcome]]), \(i) {
                # Extract state code (last 2 chars)
                state <- substr(names(augsynth_average_results[[outcome]])[i],
                                nchar(names(augsynth_average_results[[outcome]])[i])-1, nchar(names(augsynth_average_results[[outcome]])[i]))
                # Combine state ID with ATT results
                cbind(data.frame("state" = state),
                      augsynth_average_results[[outcome]][[i]]["ATT", ])
              }))
    # Sort alphabetically by state
    stateSpecificATTs <- stateSpecificATTs[order(stateSpecificATTs$state), ]
    
    # Table Formatting & Export
    df_formatted <- stateSpecificATTs
    # Format all numeric columns to 4 decimal places
    df_formatted[] <- lapply(df_formatted, function(x) {
      if (is.numeric(x)) sprintf("%.4f", x) else x
    })
    names(df_formatted)[1] <- "State"
    
    # Create table with gridExtra
    table <- tableGrob(df_formatted,
                       rows = NULL,
                       theme = ttheme_default(
                         core = list(bg_params = list(fill = "white"), fg_params = list(cex = 0.7)),   # body font size
                         colhead = list(bg_params = list(fill = "white"), fg_params = list(cex = 0.8)), # header font size
                         rowhead = list(padding = unit(c(0, 0), "mm"))
                       ))
    # Add a horizontal line below the column headers
    # The line is added after the first row (headers)
    line <- segmentsGrob(
      x0 = unit(0, "npc"), x1 = unit(1, "npc"),
      y0 = unit(0, "npc"), y1 = unit(0, "npc"),
      gp = gpar(lwd = 2)
    )
    # Add custom border line
    table <- gtable_add_grob(
      table, line, t = 1, l = 1, r = ncol(table)
    )
    # Save as high-res PNG
    gg_tbl <- as.ggplot(table) +
      theme(plot.margin = margin(0, 0, 0, 0))
    ggsave(paste0(fin.out.path, "/appx_e_table.png"), gg_tbl, width = 8, height = 4, units = "in", dpi = 300, bg = "white")
    
    # Calculate max and min values from confidence intervals
    max_val <- max(stateSpecificATTs$`95% Upper Bound`, na.rm = TRUE)
    min_val <- min(stateSpecificATTs$`95% Lower Bound`, na.rm = TRUE)
    
    # Find the symmetric limit for symmetrical display
    sym_limit <- max(max_val, abs(min_val))
    
    # Add 15% padding
    ylim_upper <- sym_limit + 0.15 * sym_limit
    ylim_lower <- -ylim_upper
    
    # Effect Estimate Visualization
    # Create the plot and save 
    appx.plot <- ggplot(stateSpecificATTs, aes(x = state, y = Estimate)) + # Point estimates
      geom_point(size = 0.7) +
      geom_errorbar(aes(ymin = `95% Lower Bound`, ymax = `95% Upper Bound`), width = 0.1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth=0.1) + # Null effect line
      theme_minimal(base_size = 10) +  # Reduce base font size
      theme(
        panel.border = element_blank(),
        axis.line.y.left = element_line(color = "black", linewidth = 0.1),
        axis.line.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 6), # Reduce title size
        axis.text.x = element_text(angle = 0,hjust = 1, size = 4, face="bold"),       # Optionally reduce axis text
        axis.text.y = element_text(size = 4, face="bold"),
        axis.title.x = element_text(size = 5),   # x-axis label font size
        axis.title.y = element_text(size = 5),   # y-axis label font size
        plot.margin = margin(10, 10, 10, 10)
      ) +
      labs(
        title = "Estimates by State with 95% Confidence Intervals",
        x = "",
        y = "Effect Estimate"
      ) +
      scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
      coord_cartesian(ylim = c(ylim_lower, ylim_upper)) # Symmetric axis
    ggsave(paste0(fin.out.path, "/appx_e_plot.png"), appx.plot, width = 8, height = 6, units = "cm", dpi = 300, bg = "white")
    
    print(paste0("Appendix generated for: ", outcome))
  }
}
#############################################################################
# This function generates time series graphs by outcome, for control vs treatment states combined 
# Total number of graphs=Total number of outcomes
get_time_series_graphs <- function(){
  # --------------------------------------------------------------
  # --- Part 1: Extract Time-Specific ATT Estimates ---
  # For each outcome, extract ATT estimates at each time point from augsynth results
  timeSpecificATTs <- structure({
    lapply(outcomes, \(outcome) {
      # Find all results matching each outcome 
      outcomeIndices <- grep(outcome, names(augsynth_results))
      # For each matching result, extract ATT estimates and calculate SE from CI
      lapply(augsynth_results[outcomeIndices], \(y) summary(y)$att |>
               transform("sd" = (lower_bound - Estimate) / -1.96)) # SE = (CI_width)/1.96
    })
  }, names = outcomes) # Preserve outcome names
  
  # --- Part 2: Enrich ATT Data with Metadata ---
  timeSpecificATTs <-
    structure({
      lapply(outcomes, \(outcome) {
        structure({
          lapply(1:length(timeSpecificATTs[[outcome]]), \(i) {
            # Get state identifier from result name (e.g., "in_CUD_proportion_CA" → "CA")
            x <- names(timeSpecificATTs[[outcome]])[i]
            tx_state <- substr(x, nchar(x) - 1, nchar(x)) # Last 2 chars = state code
            
            # Get ATT estimates for this state-outcome pair
            att <- timeSpecificATTs[[outcome]][[i]]
            # Aggregate control group sample sizes by month
            control_n <- aggregate(STATE_N ~ year_month_num,
                                   data = subset(get(tx_state), state != tx_state),
                                   sum) # Total control units per month
            # Get time intervals (TIME_SINCE_TX) for merging
            time_int <- subset(get(tx_state), state != tx_state) %>% distinct(year_month_num, TIME_SINCE_TX)
            
            # Merge three components:
            # 1. ATT estimates (att)
            # 2. Control group sample sizes (control_n)
            # 3. Time mapping (time_int)
            merge(
              merge(att, control_n,
                    by.x = "Time", by.y = "year_month_num", all = T),
              time_int, by.x="Time", by.y="year_month_num"
            )
          })
        }, names = names(timeSpecificATTs[[outcome]])) # Preserve state identifiers
      })
    }, names = outcomes) # Preserve outcome names
  
  # Stack state-specific results into outcome-specific dataframes
  timeSpecificATTs_stacked <- lapply(timeSpecificATTs, \(x) {
    do.call(rbind, x) # Combine state results vertically
  })
  
  # Calculate Overall ATTs (Average Treatment Effects)
  overallATTs <- lapply(timeSpecificATTs_stacked, \(z) {
    sapply(-36:35, \(tp) { # For each month relative to treatment (-36 to +35)
      tpIndices <- which(abs(z$TIME_SINCE_TX - tp) < .001) # Find matching timepoints
      x <- z$Estimate[tpIndices] # Extract ATT estimates
      # Weighting strategy:
      if (tp < 0) { # Pre-treatment period
        w <- z$STATE_N[tpIndices] # Weight by control group sample size
        stdev <- NA # No SE calculated pre-treatment
      } else { # Post-treatment period
        w <- z$sd[tpIndices]^(-2) # Inverse variance weighting (1/SE²)
        stdev <- sqrt(1 / sum(w, na.rm = T)) # SE of weighted mean
      }
      c(tp, weighted.mean(x, w), stdev) # Return: [Time, ATT, SE]
    })
  })
  
  # Calculate Confidence Limits
  overall_limits <- lapply(overallATTs, \(z) {
    rbind(z[2, ] - 1.96 * z[3, ], # Lower bound (ATT - 1.96*SE)
          z[2, ] + 1.96 * z[3, ]) # Upper bound (ATT + 1.96*SE)
  })
  
  # Compute Treatment Group Means
  overallTxMeans <- structure({
    lapply(outcomes, \(outcome) {
      # Reshape treated states' data to wide format (time × state)
      out <- do.call(rbind, lapply(tx_states, \(tx_state) {
        subset(get(tx_state), state == tx_state,
               select = c("TIME_SINCE_TX", "state", outcome))
      })) |>
        reshape(direction = "wide", idvar = "TIME_SINCE_TX",
                timevar = "state")
      
      outN <- do.call(rbind, lapply(tx_states, \(tx_state) {
        subset(get(tx_state), state == tx_state,
               select = c("TIME_SINCE_TX", "state", "STATE_N"))
      })) |>
        reshape(direction = "wide", idvar = "TIME_SINCE_TX",
                timevar = "state")
      
      out <- out[order(out$TIME_SINCE_TX), ]
      outN <- outN[order(outN$TIME_SINCE_TX), ]
      
      # Calculate sample-size weighted average across treated states
      avg_tx_state <- as.data.frame(
        cbind(out[, 1],
              do.call(rbind, lapply(1:nrow(out), \(i) {
                weighted.mean(x = out[i, -1], w = outN[i, -1], na.rm = T)
              }))))
      names(avg_tx_state) <- c("TIME_SINCE_TX", outcome)
      avg_tx_state
    })
  },
  names = outcomes)
  
  #-------------------------------------------------------------------
  
  # Function that plots time series graph for each outcome aggregated across states
  plotFun <- function(outcome, dis, overallTxMeans,
                      overall_limits, overallATTs) {
    
    # Adjust time index to start at 0 
    overallTxMeans[[outcome]][["TIME_SINCE_TX"]] <- overallTxMeans[[outcome]][["TIME_SINCE_TX"]] + 36
    # Calculate y-axis limits conditionally based on outcome type
    max_y <- max(overallTxMeans[[outcome]][, outcome], na.rm = TRUE)
    ylim_upper <- 10 * max_y # Set upper limit to 10 times the max y value
    ylim_lower <- 0 # Always start y-axis at 0
    
    # Configure plot layout
    # Set small margins: bottom, left, top, right
    par(mar = c(2.5, 2.5, 1.5, 1), oma = c(0, 0, 0, 0))
    plot(
      x = 0:71, y = overallTxMeans[[outcome]][, outcome],
      type = "n", lwd = 2, # Empty plot frame
      ylim = c(ylim_lower, ylim_upper),
      xlab = "Months",
      ylab = "Proportion of patients in a given month",
      main = "",
      xaxt = "n",          # suppress x-axis, will add manually
      yaxt = "n",
      las = 1,             # Horizontal y-axis labels
      cex.lab = 0.7,       # axis label size (x/y labels)
      cex.axis = 0.7,      # y-axis tick label size (x handled below)
      cex.main = 0.6,      # main title size
      cex.sub = 0.6,       # subtitle size, if used
      mgp = c(1.3, 1, 0),  # Axis label positioning
      bty = "n"            # No box around plot
    )
    
    # Add policy implementation marker
    # Vertical line at law implementation
    abline(v = 36, lty = "dotted")
    
    # Draw confidence interval polygon
    conflim <-
      data.frame(x = c(0:71, rev(0:71)),
                 y = c(overallTxMeans[[outcome]][, outcome] -
                         overall_limits[[outcome]][1, ],
                       rev(overallTxMeans[[outcome]][, outcome] -
                             overall_limits[[outcome]][2, ]))) |>
      na.omit()
    polygon(conflim$x, conflim$y,
            col = "#0000FF33", border = NA) # Semi-transparent blue
    
    # Plot observed and counterfactual trends
    # Observed trend
    lines(x = 0:71, y = overallTxMeans[[outcome]][, outcome], lwd=1)
    # Counterfactual trend
    lines(x = 0:71, y = overallTxMeans[[outcome]][, outcome] -
            overallATTs[[outcome]][2, ],
          col = "blue", lty = 2, lwd=1) # Semi-transparent blue
    par(mgp = c(1.3, 0.3, 0))  # Adjust as needed
    
    # Customize axes
    axis(1,
         at = c(seq(0, 24, 12), seq(36, 72, 12)),
         labels = c("0", "12", "24", "36 (Law)", "48", "60", "72"),
         cex.axis = 0.7,   # x-axis tick label size
         lwd = 0.5,        # thin axis line
         lwd.ticks = 0.5)   # thin tick marks
    par(mgp = c(1.3, 0.5, 0))  # Adjust as needed
    axis(2, cex.axis = 0.5, lwd = 0.5, lwd.ticks = 0.5)
    # Add title closer to the plot (e.g., line = 1)
    title(main = dis, line = 1, cex.main = 0.7)  # smaller line value = closer
    # Draw only left and bottom borders, thin
    box(bty = "l", lwd = 0.5)
    legend("top",
           legend = c("Treatment States",
                      "Comparison States"),
           lty = c(1, 2),
           col = c("black", "blue"), bg = NA,
           horiz = TRUE, cex = 0.5)
  }
  
  # Description for each outcome to include in the figure
  outcome.fig.desc <- c("in_CUD"="Any inpatient CUD treatment", 
                        "ot_CUD"="Any outpatient CUD treatment",  
                        "CUD_initiation"="CUD treatment Initiation",
                        "in_OUD"="Any inpatient OUD treatment", 
                        "ot_OUD"="Any outpatient OUD treatment",  
                        "OUD_initiation"="OUD treatment Initiation", 
                        "OUD_med"="Any medication for OUD (MOUD)",
                        "ovd_overall_cannabis"="CUD overdose", 
                        "ovd_overall_opioid"="OUD overdose")
  
  # Loop across outcomes and generate and save the figure 
  for (outcome in outcomes){
    # Define path for each outcome
    fin.out.path <- paste0(base.path, "final_exhibits/", outcome)
    # Create a plot 
    png(paste0(fin.out.path,"/time_Series_plot.png"), width = 3.4, height = 3.4, units = "in", res = 300)
    # Get the figure and add to the plot
    plotFun(outcome, dis =outcome.fig.desc[[gsub("_proportion", "", outcome)]],
            overallTxMeans = overallTxMeans, overall_limits = overall_limits,
            overallATTs = overallATTs)
    # Clear the plot space
    dev.off()
    print(paste0("Time series plot generated for", outcome))
  }
}
#############################################################################


