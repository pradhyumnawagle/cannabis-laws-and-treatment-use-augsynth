## List all packages to install
packs <- c("haven", "dplyr", "tidyverse", "ggplot2", "zoo", "readxl", "re",
           "xlsx", "skimr", "qcc", "magick", "openxlsx", "here", "ggpubr", "webshot2", "ggplotify",
           "gridExtra", "grid", "gtable", "gt") 
# Install all packages
for (p in setdiff(packs, installed.packages())) {
  install.packages(p)
}
## Load all libraries
library(zoo)
library(haven)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(re)
library(readxl)
library(skimr)
library(qcc)
library(magick)
library(openxlsx)
library(gridExtra)
library(grid)
library(gtable)
library(gt)
library(webshot2)
library(ggplotify)

## Declare and create(if/as needed) input and output paths for datasets
base.path <- "C:/Users/prw4002/Weill Cornell Medicine/Kayla Tormohlen - State Medical Cannabis Laws- VRDC/State Medical Cannabis Laws-Working Files/"
in.path   <- paste0(base.path, "/12m CE_SOURCE_COPY") # Input files path
oinp.path <- paste0(base.path, "/other_inputs") # Other inputs path
outc.path <- paste0(base.path, "/outcome_workbooks") # Outcome workbooks path
log.path  <- paste0(base.path, "/logs") # Log path

##################### Parameters ###############################################
## List states to run
state.list <- as.list(c( "FL", "MD", "MN", "NH", "NY", "OK", "PA"))
## List cohorts to run
cohort.list <- c("outcome")
## List outcomes to evaluate
outcome.vars <- c("in_CUD", "ot_CUD", "CUD_initiation",
                  "in_OUD", "ot_OUD", "OUD_initiation", "OUD_med",
                  "ovd_overall_cannabis", "ovd_overall_opioid" )
## List all additional parameters needed to control the augsynth run
small.cell.rep = 5 # Replace small cell suppression by 5
################################################################################
## Evaluate all variables needed for the model using parameters above
## List cohort names in to show in figures
cohort.names <- c("ARTHRITIS"="Arthritis","HEADACHE"="Headache","outcome"="Overall", 
                  "FIBROMYALGIA"="Fibromyalgia", "NEUROPATHIC"="Neuropathic", "LBP"="LBP", 
                  "black" = "Race-black","hispanic"="Race-hispanic","other"="Race-other","white"="Race-white")
## Create a list of cohorts and states combined to run loops 
cohort.state.list <- as.list(outer(cohort.list, state.list, paste))
## List of race specific cohorts
race.cohort.list <- c("black", "hispanic", "other", "white")
## List of non-race specific cohorts
non.race.props.list <- c("female_proportion", "mental_illness_proportion")
## List of demograpic variables
demo.vars     <- c( "female", "mental_illness", "race_white", "race_black", "race_other", "race_hispanic")
## Proportion names of demographic variables
demo.props  <- paste0(demo.vars, "_proportion")

## Define outcome variable proportions and description for each outcome
outcome.props <- paste0(outcome.vars, "_proportion")
outcome.desc <- c("in_CUD"="% of patients with any inpatient CUD treatment", 
                   "ot_CUD"="% of patients with any outpatient CUD treatment",  
                   "CUD_initiation"="% of patients newly initiating CUD treatment",
                   "in_OUD"="% of patients with any inpatient OUD treatment", 
                   "ot_OUD"="% of patients with any outpatient OUD treatment",  
                   "OUD_initiation"="% of patients newly initiating OUD treatment", 
                   "OUD_med"="% of patients with any medication for OUD (MOUD)",
                   "ovd_overall_cannabis"="% of patients with any CUD overdose utilization", 
                   "ovd_overall_opioid"="% of patients with any OUD overdose utilization")

###################################################################
## Run - stratify all cohorts 
## This script will stratify all input files by state and cohort

source("stratify_state_cohorts.R")
#------------Open a log file--------------------------
log_file <- paste0(log.path, "/statify_state_cohorts_log.txt")
log_conn <- file(log_file, open = "wt")
sink(log_conn, append = TRUE, type = "output")  # Capture normal output
sink(log_conn, append = TRUE, type = "message")  # Capture warnings and messages
# Redirect errors to the same file
options(warn = 1)  # This makes sure warnings are written out as they occur
#------------Run the program-----------------------------
stratify_state_cohorts()
#------------Close the log file--------------------------
sink()  # Reset the output capture
sink(type = "message")  # Reset message capture
close(log_conn)
#####################################################################
## Run augsynth for all cohorts
## This script will run augsynth, save results, average outcomes, and trends figure in the respective state-cohort-outcome subfolder

source("cohort_augsynth.r")
#------------Open a log file--------------------------
start.time <- Sys.time()
log_file <- paste0(log.path, "/augsynth_log.txt")
log_conn <- file(log_file, open = "wt")
sink(log_conn, append = TRUE, type = "output")  # Capture normal outputs
sink(log_conn, append = TRUE, type = "message")  # Capture warnings and messages
# Redirect errors to the same file
options(warn = 1)  # This makes sure warnings are written out as they occur
#------------Run the program-----------------------------
run_augsynth_for_cohorts() 
#------------Close the log file--------------------------
end.time <- Sys.time()
runtime <- end.time - start.time
cat("Task total runtime:", runtime, "minutes\n")
sink()  # Reset the output capture
sink(type = "message")  # Reset message capture
close(log_conn)
#####################################################################
# Generate Results
# 1) Get average augsunth results for all outcmoes and states within a cohort
# 2) Get pre and post treatement average augsynth results 
# 3) Get table 1 for the main exhibit
# 4) Get the time seried graphs for the main exhibit
# 5) Get Appendices for the main exhibit

source("aggregate_results.r")
#------------Open a log file--------------------------
start.time <- Sys.time()
log_file <- paste0(log.path, "/aggregate_results.txt")
log_conn <- file(log_file, open = "wt")
sink(log_conn, append = TRUE, type = "output")  # Capture normal outputs
sink(log_conn, append = TRUE, type = "message")  # Capture warnings and messages
# Redirect errors to the same file
options(warn = 1)  # This makes sure warnings are written out as they occur
#------------Run the program-----------------------------
## Instantiate lists to store results
augsynth_results <- list()
cohort_means <- list()
cohort_sds   <- list()
#------------------------------------------------------
get_augsynth_avg_results()
get_pre_post_augsynth_results()
table1 <- get_table1()
get_appx_exhibits()
get_time_series_graphs()
#------------Close the log file--------------------------
end.time <- Sys.time()
runtime <- end.time - start.time
cat("Task total runtime:", runtime, "minutes\n")
sink()  # Reset the output capture
sink(type = "message")  # Reset message capture
close(log_conn)
##################--------------
# Optional
# Save all results in workbooks by each cohort
# That includes adjusted and unadjusted models in the workbook
# Each workbook is for a cohort, and each tab is an outcome
source("save_results_to_excel.r")
save_results_to_workbook()
#####################################################

