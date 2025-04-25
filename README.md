# cannabis-laws-and-treatment-use-augsynth
This repository stores code for the study: The impact of medical cannabis laws on cannabis and opioid use disorder treatment and overdose-related healthcare utilization among adults with chronic non-cancer pain. 
Included are R scripts for data cleaning, running the Augmented Synthetic Control model, generating summary statistics, aggregating results across states and time, and plotting time series graphs.


- [main_caller.R](main_caller.R): 
Primary caller file where all variables, parameters and libraries are defined. 
This program calls all other program in order(top-down).
Run all parts from top-bottom to run the entire program.


- [stratify_state_cohorts.R](stratify_state_cohorts.R): 
This is the first program called. This program organizes all cohort data files based on the study_type variable. The study_type variable combines the specific type of study with the name of each treatment state. Each cohort file contains a total of 9,072 records, with 1,296 entries for each of the 7 study types. Within each study type, there are 72 records for each of the 17 control states, as well as 72 records for the corresponding treatment state. These 72 records represent each month in the study period. For example, for the study type overall_FL, there will be 72 records for each of the 17 comparison states and 72 for Florida itself, totaling 1,296 records. This same structure applies to all treatment states.


- [cohort_augsynth.R](cohort_augsynth.R):
This program runs augsynth for all stratified cohort state files. This program calls run_and_plot to run augsynth and plot the augsynth results.  


- [run_and_plot.R](run_and_plot.R):
This utility program runs augsynth, saves results and plots ATT figures over time.


- [xaugsynth.r](xaugsynth.r):
This program is a utility program which generates average ATT values from augsynth results, and creates differences plot. This code is derived from [this](https://github.com/nickseewald/opioid-prescribing-augsynth/blob/main/xaugsynth.R) original code, which extends the [augsynth](https://github.com/ebenmichael/augsynth) package. 


- [aggregate_results.R](aggregate_results.R):
This program aggregates all augsynth results and generates Table 1, time series graphs and appendix for the the main exhibits.


- [save_results_to_excel.R](save_results_to_excel.R):
This optional program saves all augsynth results and differences plot in excel file for each cohort.
