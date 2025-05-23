
## Run augsynth model for all cohorts
run_augsynth_for_cohorts <- function(){
  # Augsynth runs from this R program
  source("run_and_plot.r")
  
  # Read the implementation date dataset
  imp.dates.exc <- paste0(oinp.path, "/state_cannabis_laws_implementation_dates.xlsx")
  imp.dates <- read_excel(imp.dates.exc) %>% drop_na()
  
  # Run through all state cohort list
  lapply(cohort.state.list, function(x){
    # Get the current state name and cohort name
    cohort.name <- strsplit(x, " ")[[1]][1]
    state.name <- strsplit(x, " ")[[1]][-1]
    # Declare and create(if needed) input and output paths for datasets
    in.path   <- paste0(base.path, "state_cohorts/", cohort.name, "/", state.name)
    
    # Run the loop for each outcome variable
    for (oc in outcome.vars){
      # Create a folder to store results and figures for each outcome
      out.path  <- paste0(base.path, "/state_cohorts/", cohort.name, "/", state.name,"/", oc)
      dir.create(out.path)
      # Delete all files in the folder
      unlink(file.path(out.path, "*"), recursive = FALSE, force = TRUE)
      
      outcome <- paste0(oc, "_proportion")
      
      # Read the stratified file
      inp.file <- paste0(in.path, "/",cohort.name, "_", state.name, ".xpt")
      data <- read_xpt(inp.file)
      # Assign cohort state as the state for which the file is read
      data$cohort_state <- state.name
      # Create a final dataset before running augsynth
      # Merge implementation dates
      # Assign each row with treatment indicator(1/0)
      # Change year month num to date format
      # Replace all small cell suppression across demographic variables by the number assigned
      # Replace all small cell suppression for the outcome variable
      data.fin <- data %>%
        left_join(imp.dates, by = c("cohort_state" = "Treatment_State")) %>%
        mutate(treat_year_month = as.integer(format(Implementation_date_rounded, "%Y%m"))) %>%
        mutate(treatment = ifelse(year_month_num>=treat_year_month & state==cohort_state, 1, 0)) %>%
        mutate(year_month_num = as.yearmon(as.character(year_month_num), format = "%Y%m"))
        
      
      # # Write the final file into the output folder for aggregation purpose
      write_xpt(data.fin, path = paste0(in.path, "/", oc , "/", cohort.name, "_", state.name, "_clean.xpt"))

      # # Create a list of valid columns based on columns present in the dataset to generate a formula
      if (cohort.name %in% race.cohort.list){
        non.race.props.list <- append(non.race.props.list, "mean_age")
        valid.columns <- paste(intersect(non.race.props.list, colnames(data.fin)), collapse=" + ")
      } else {
        demo.props <- append(demo.props, "mean_age")
        valid.columns <- paste(intersect(demo.props, colnames(data.fin)), collapse=" + ")
      }

      ## Run augmented synthetic control method using ridge regression to augment, as
      ## recommended by Ben-Michael, et al. (2021).
      ## Use the run_augsynth function from run_and_plot script to run augsynth, save results and plots
      ## Run twice: one for adjusted and one for unadjusted model
      model <- as.formula(paste(outcome, "~ treatment |", valid.columns))
      print(paste0("adj model - ", cohort.name, " " , state.name, " ", outcome))
      print(model)
      run_augsynth(dataset=data.fin,
                   eqn=model,
                   adjustment="adj",
                   state.name=state.name,
                   cohort=cohort.names[[cohort.name]],
                   outcome=oc,
                   outcome_desc=outcome.desc[[oc]],
                   output_path = out.path
      )
      print("##########################################")
      
      # Unadjusted model for augsynth: commented out now
      # model2 <- as.formula(paste0(outcome, " ~ treatment"))
      # print(paste0("unadj model - ", cohort.name, " " , state.name, " ", outcome))
      # print(model2)
      # run_augsynth(dataset=data.fin,
      #              eqn=model2,
      #              adjustment="unadj",
      #              state.name=state.name,
      #              cohort=cohort.names[[cohort.name]],
      #              outcome=oc,
      #              outcome_desc=outcome.desc[[oc]],
      #              output_path = out.path
      # )
      #print("####################################################################################")
      print("####################################################################################")
    }
  })
}

