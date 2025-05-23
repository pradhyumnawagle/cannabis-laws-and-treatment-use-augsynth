
## This function will stratify sate cohorts by study type
stratify_state_cohorts <- function(){
  
  ## Define and create output directories for each cohort-state combo
  out.path  <- paste0(base.path, "/state_cohorts")
  dir.create(out.path)
  
  # Initialize an empty list to store summaries
  summary_list <- list()
  # Function to Create directory if not present
  create_dir <- function(path, folder){
    # Create a file path if not already there
    dir <- file.path(path, folder)
    if (!dir.exists(dir)) dir.create(dir)
  }
  
  # Function to extract the dynamic string
  extract_string <- function(filename) {
    # Use regex to capture the string between "twelve_m_" and "_all_state"
    sub(".*twelve_m_(.*?)_all_state.*", "\\1", filename)
  }
  
  # Read the implementation data dataset to get cannabis law implementation for each state
  imp.dates.exc <- paste0(oinp.path, "/state_cannabis_laws_implementation_dates.xlsx")
  imp.dates <- read_excel(imp.dates.exc) %>% drop_na()
  
  # List all files to break into cohorts
  files <- list.files(path=in.path, pattern="*.sas7bdat", full.names=TRUE, recursive=FALSE)
  # Work through one file at a time and break it into state-cohort level files
  for (x in files) {
    # Read the table
    data <- read_sas(x)
    cohort.name <- extract_string(x)
    print(cohort.name)
    # Keep the columns in the table or columns to generate proportion
    vars.to.eval <- c(demo.vars, outcome.vars)
    curr.vars <- intersect(vars.to.eval, names(data))
    
    ## Output a summary table with counts of total benes, outcomes and small cell suppression
    ## Done here so that we capture small cell suppressed values before replacing
    outc.summ <- data %>%
      filter(state %in% state.list) %>%
      mutate(
        ot_cr_CUD = ot_CUD + cr_CUD,
        ot_cr_OUD = ot_OUD + cr_OUD,
        across(all_of(outcome.vars), ~ ifelse(. == 99999, 1, .))) %>%
      group_by(state) %>%
      summarise(
        mean_bene_count = mean(enrolled_bene),
        min_bene_count = min(enrolled_bene),
        max_bene_count = max(enrolled_bene),

        across(all_of(outcome.vars),
               list(mean= ~mean(., na.rm=TRUE),
                    scs = ~ sum(. == 1, na.rm = TRUE),
                    scs_pct = ~ sum(. == 1, na.rm = TRUE) / 72),
               .names = "{col}_{fn}"),
      # Add the count of rows where sud_overall == 99999 and sud_init == 0
      sud_diff = sum(SUD_overall == 1 & (SUD_initiation != 0 & SUD_initiation != 1 )),
      cud_diff = sum(CUD_overall == 1 & (CUD_initiation != 0 & CUD_initiation != 1)),
      oud_diff = sum(OUD_overall == 1 & (OUD_initiation != 0 & OUD_initiation != 1)),

      cohort = cohort.name,
      .groups = 'drop')
    # Append to summary list
    summary_list[[x]] <- outc.summ
    
    # This is where actual data manipulation happens before stratifying
    # Replace small cell suppressed values in each of those columns as NA and generate proportion out of total benes
    data <- data %>%
      mutate(across(all_of(curr.vars), ~ ifelse(. == 99999, small.cell.rep, .))) %>%
      mutate(across(all_of(curr.vars), ~ (. / enrolled_bene)*100, .names = "{.col}_proportion")) %>%
      mutate(mean_age = (age_at_end_1+age_at_end_2+age_at_end_3+age_at_end_4
                         +age_at_end_5+age_at_end_6)/enrolled_bene) %>%
      # Create a column that has the sum of outpatient plus carrier claims
      # Calculate the proportion for that variable and perform small cell reps
      mutate(
        across(c("cr_CUD", "cr_OUD", "ot_OUD", "ot_CUD"), ~ ifelse(. == 99999, small.cell.rep, .)),
        ot_cr_CUD = ot_CUD + cr_CUD,
        ot_cr_OUD = ot_OUD + cr_OUD,
        
        ot_cr_CUD_proportion = (ot_cr_CUD/enrolled_bene)*100,
        ot_cr_OUD_proportion = (ot_cr_OUD/enrolled_bene)*100,
      )
    
    
    # Split the table into tiny tibbles by study type using the split function
    strat_data <- split(data, data$study_type)
    # For each state-cohort, loop through and work as needed
    for (name in names(strat_data)) {

      # Extract the cohort name from cohort using substring
      cohort.name<- substr(name, 1, nchar(name)-3)
      state.name <- substr(name, nchar(name)-1, nchar(name))

      # Create the cohort-state folder

      dir <- create_dir(out.path, cohort.name)
      dir <- create_dir(paste(out.path, cohort.name, sep="/"), state.name)
      ## Delete all files in the folder
      unlink(file.path(paste0(out.path, "/", cohort.name, "/", state.name), "*"), recursive = FALSE, force = TRUE)

      # Write the file into the path
      write_xpt(strat_data[[name]], path = paste0(out.path, "/", cohort.name, "/", state.name, "/", name, ".xpt"))
      #save(strat_data[[name]], file = paste0(out.path, "/", cohort.name, "/", state.name, "/", name, ".xpt"))
      cat("File saved: ",name , " with ", nrow(strat_data[[name]]),  "observations. \n")
    #   
    #   ##########################################################################
    ## Add cohort specific state in the dataset and add if each row is for treatment state or comparison pool
      strat_data[[name]]$cohort_state <- state.name
      strat_data[[name]]$treatment_state <- ifelse(strat_data[[name]]$state==strat_data[[name]]$cohort_state,strat_data[[name]]$state,"Comparison Pool")

      ## Get the variable distribution for each file and output the distribution
      count_zeros <- function(x) {
        sum(x == 0, na.rm = TRUE)
      }
      skim_func <- skim_with(numeric=sfl(hist=NULL, n_total=length, n_zero = ~ count_zeros(.)) )
      var.dist.table <- strat_data[[name]]  %>%
        group_by(treatment_state) %>%
        select(any_of(vars.to.eval), treatment_state)
      var.dist.table <- skim_func(var.dist.table)
      colnames(var.dist.table) <- gsub("^(numeric\\.|skim\\_)", "", colnames(var.dist.table))
      write_csv(var.dist.table, paste0(out.path, "/", cohort.name, "/", state.name, "/", name, "_variable_distribution.csv"))
      ##########################################################################
      
      ## The code below produces pre-period trends graph for overall cohort
      # if (cohort.name=="outcome"){ 
      #   # Merge the implemented dates to keep only pre period data
      #   imp_date_merge <- strat_data[[name]] %>%
      #                     left_join(imp.dates, by = c("cohort_state" = "Treatment_State")) %>%
      #                     mutate(treat_year_month = as.integer(format(Implementation_date_rounded, "%Y%m")))
      #   # Group by the treatment vs comparison state and summarize the demographics
      #   demo.summary <- imp_date_merge %>%
      #                   filter(year_month_num<treat_year_month) %>%
      #                   select(all_of(demo.props), mean_age, treatment_state) %>%  
      #                   group_by(treatment_state) %>%
      #                   summarize(across(everything(), \(x) mean(x, na.rm = TRUE))) %>%
      #                   mutate(across(where(is.numeric), ~ round(.x, 2)))
      #   # Write the demgraphic summary file into the path
      #   #write_xpt(demo.summary, path = paste0(out.path, "/", cohort.name, "/", state.name, "/", name, "_demographic_summary.xpt"))
      #   write_csv(demo.summary, paste0(out.path, "/", cohort.name, "/", state.name, "/", name, "_demographic_summary.csv"))
      #   
      #   outcomes.summary <- imp_date_merge %>%
      #                       select(all_of(outcome.props), treatment_state, year_month_num, treat_year_month) %>%  
      #                       group_by(treatment_state, year_month_num) %>%
      #                       summarize(across(everything(), mean, na.rm = TRUE))
      #   
      #   state.treat.year.month <- imp.dates %>% filter(Treatment_State==state.name) %>% pull(Implementation_date_rounded)
      #   state.treat.year.month <- format(state.treat.year.month, "%Y%m")
      #   tick_labels <- as.character(outcomes.summary$year_month_num[seq(1, length(outcomes.summary$year_month_num), length.out = 25)])
      #   
      #   for (var in outcome.props){
      #     out.summ.plot <- outcomes.summary%>%
      #       arrange(year_month_num, treatment_state) %>%
      #       mutate(year_month_num = as.character(year_month_num)) %>%
      #       ggplot(aes_string(x = "year_month_num", y=var, color = "treatment_state", group = "treatment_state")) +
      #       geom_line() +
      #       scale_x_discrete(breaks = tick_labels) +
      #       geom_vline(xintercept = state.treat.year.month, color = "black") +
      #       theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
      #       labs(title = paste0(var, " by Treatment State vs Comparison Pool"),
      #            x = "Year-Month",
      #            y = var)
      #     out.summ.plot <- out.summ.plot+
      #                     annotate(x=state.treat.year.month,y=+Inf,label=paste0("Treatment (", state.treat.year.month, ")"),vjust=2,geom="label")
      #     ggsave(paste0(out.path, "/", cohort.name, "/", state.name, "/", var, ".jpeg"), plot = out.summ.plot, width = 8, height = 6, dpi = 300)
      #   }
      #   # Write the demographic summary file into the path
      #   #write_xpt(outcomes.summary, path = paste0(out.path, "/", cohort.name, "/", state.name, "/", name, "_outcomes_summary.xpt"))
      #   outcomes.summary <- outcomes.summary %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
      #   write_csv(outcomes.summary, paste0(out.path, "/", cohort.name, "/", state.name, "/", name, "_outcomes_summary.csv"))
      #   
       }
  }
  
  
  # # Combine all summaries into one data frame
  final_summary <- bind_rows(summary_list)
  cat("rows:", nrow(final_summary), "\n")
  output_file <- paste0(outc.path, "/Outcomes_summary.xlsx")
  write.xlsx(final_summary, output_file)
  # Inform the user
  cat("Excel file saved as:", output_file, "\n")
}

  


 




