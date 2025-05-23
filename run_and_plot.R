## Run augsynth function below runs Augmented synthetic Control model, saves outputs and plots
run_augsynth <- function(dataset, eqn, adjustment, state.name, cohort, outcome, outcome_desc, output_path) {
  

  #################Run Augsynth #######################################  
  source("xaugsynth.r")
  
  # Run augmented synthetic model
  covsyn <-   augsynth(eqn,
                       unit = state,
                       time = year_month_num,
                       data = dataset,
                       progfunc = "ridge",  #use ridge regression to impute control outcomes
                       scm = TRUE,          #use scm weighting function
                       fixedeff = TRUE      #include unit fixed effects
  )
  
  # Write the final file into the output folder for aggregation purpose
  saveRDS(covsyn, paste0(output_path, "/augsynth_results_", adjustment,".xpt"))
  # Store results for aggregation
  analysis_name <- paste(state.name, cohort, sep = "_")
  # Store augsynth weights
  weights <- covsyn$weights
  augsynth_weights <- data.frame("STATE_CODE" = rownames(weights),
                                 "weight" = weights)
  write.csv(augsynth_weights, paste0(output_path, "/augsynth_weights_", adjustment,".csv"), row.names=TRUE)
  
  ##################Plot difference graph ###########################
  # Plot the percent points difference plot and save
  # Set global plotting parameters
  par(ylim = c(-1, 1))
  plot_path <- paste0(output_path, "/synth_differences_", adjustment,".jpeg")
  jpeg(plot_path, width=600, height=430)
  plot_obj <- plot(covsyn) +
    labs(x = "Month-Year",
         y = "Percentage point change",
         title=paste0("Percent Point difference: ", state.name, " and \nSynthetic Control for cohort: ", cohort)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 6),
          plot.margin = margin(t = 40, r = 10, b = 20, l = 20),
          axis.title.x = element_text(vjust = -8, size=12),  # Move the x-axis title below the ticks
          axis.text.x = element_text(vjust = -3, size=12),
          axis.title.y = element_text(vjust = 6, size=12))
  # Modify y-axis limits using ggplot2
  plot_obj <- plot_obj + ggplot2::scale_y_continuous(limits = c(-1, 1))
  print(plot_obj)
  # Center the title
  dev.flush()  # Force the plot to render before closing
  dev.off()
  # Reset global parameters after plotting
  par(ylim = NULL)
  
  #############Save Augsynth results #################################
  # Write and save the average outcomes dataset from the model
  avg.outcomes <- averageOutcomes(covsyn)
  write.csv(avg.outcomes, paste0(output_path, "/augsynth_summ_", adjustment,".csv"), row.names=TRUE)
  
  ###########Plot trends graph ########################################
  # Calculate the max y for plotting purpose
  avg.outcomes[, 3][is.infinite(avg.outcomes[, 3])] <- NA
  att_summ <- summary(covsyn)$att
  # Calculate the max value, and add 0.5
  y.max <- max(max(att_summ$Estimate, na.rm = TRUE)) + 1
  y.min <- min(att_summ$Estimate, na.rm = TRUE) - 1
  # Plot the difference in values with time between treatment state and synthetic control
  res_path <- paste0(output_path, "/results_", adjustment,".jpeg")
  jpeg(res_path, width = 600, height = 450)
  potentialplot(covsyn,
                xlab = "Month-Year",
                ylab = outcome_desc,
                main = paste0(state.name, ": ", outcome_desc, " \nper month for Cohort: ", cohort),
                legendParams = list(legend = c(state.name, "Synthetic Control"),
                                    horiz = T, cex = 1))
  dev.off()
  
  ########Combine 2 figures for workbook purpose(Optional)#############
  # Combine 2 plots in a single one to store it on results workbook
  q1 <- image_read(plot_path)
  q2 <- image_read(res_path)
  img <- c(q1, q2)
  comb <- image_append(img)
  comb_path <- paste0(output_path, "/combined_", adjustment,".png")
  image_write(comb, path=comb_path)
  
}
