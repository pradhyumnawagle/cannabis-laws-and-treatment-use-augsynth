# Write results to excel workbook
save_results_to_workbook <- function(){
  
  
  # This function sets border for a table where the outer borders are bold
  set_table_border <- function(start.row, n.rows, start.col, n.col, fill.inner=TRUE, style.headers=TRUE ){
    
    # Define top border style
    top_border_style <- createStyle(
      border = "top", 
      borderColour = "black", 
      borderStyle = "medium"
    )
    
    # Define bottom border style
    bottom_border_style <- createStyle(
      border = "bottom", 
      borderColour = "black", 
      borderStyle = "medium"
    )
    
    # Define left border style
    left_border_style <- createStyle(
      border = "left", 
      borderColour = "black", 
      borderStyle = "medium"
    )
    
    # Define dotted inner border style (for the inside of the table)
    inner_border_style <- createStyle(
      border = "TopBottomLeftRight", 
      borderColour = "black", 
      borderStyle = "dotted",
      numFmt = "0.000",
      wrapText = TRUE,
      halign = "center",
      valign = "center"
    )
    
    # Define dotted header border style (for the header of the table)
    header_style <- createStyle(
      wrapText = TRUE,
      halign = "center",
      valign = "center",
      textDecoration = "italic",
      fgFill = "#E0E0E0"
    )
    
    
    # Apply the solid outer border to the entire table range (first and last row and column)
    addStyle(wb, sheet = outcome, style = left_border_style, rows = (start.row+1):(start.row+n.rows+1), cols = n.col+1, gridExpand = TRUE)
    addStyle(wb, sheet = outcome, style = bottom_border_style, rows = start.row, cols = start.col:n.col, gridExpand = TRUE)
    if (fill.inner) {
      addStyle(wb, sheet = outcome, style = inner_border_style, rows = (start.row+1):(start.row+n.rows+1), cols = start.col:n.col, gridExpand = TRUE)
    }
    addStyle(wb, sheet = outcome, style = top_border_style, rows = start.row+n.rows+2, cols = start.col:n.col, gridExpand = TRUE)
    if (style.headers) {
      addStyle(wb, sheet = outcome, style = header_style, rows = start.row+1, cols = start.col:n.col, gridExpand = TRUE)
      addStyle(wb, sheet = outcome, style = header_style, rows = (start.row+1):(start.row+n.rows+1), cols = start.col, gridExpand = TRUE)
    }
  }
  
  # Function to add header to the table 
  add_header <- function(row, col, text){
    
    # Define border style
    header_style <- createStyle(
      fontSize=12,
      textDecoration = "bold",
      borderColour = "black",
      border= "bottom",
      borderStyle = "medium"
    )
    writeData(wb, sheet=outcome,x=text, startCol = col, startRow = row)
    addStyle(wb, sheet = outcome, style = header_style, rows = row, cols = col)
  }
  
  # Output all files to excel 
  for(cohort in cohort.list){
    
    # Create a new Excel workbook
    wb <- createWorkbook()
    
    for(outcome in outcome.vars){
      # Add a worksheet to the workbook
      addWorksheet(wb, outcome, gridLines = FALSE )
      # Assign start row
      wks.row <- 3
      # Write tables for each state
      for (state in state.list){
        in.path  <- paste0(base.path, "state_cohorts/", cohort, "/", state,"/", outcome)
        
        # Write the tables to the worksheet
        tab.path.adj <- paste0(in.path, "/augsynth_summ_adj.csv")
        tab.path.unadj <- paste0(in.path, "/augsynth_summ_unadj.csv")
        df1 <- read_csv(tab.path.adj, show_col_types = FALSE)
        # Assign first colname as empty
        colnames(df1)[1] <- ""
        df2 <- read_csv(tab.path.unadj, show_col_types = FALSE)
        colnames(df2)[1] <- ""
        
        # Write 2 datasets: adjusted and unadjusted
        writeData(wb, sheet = outcome, x = df1, startCol = 1, startRow = wks.row+1)
        # Apply the dotted inner border to the inside cells (excluding the outermost cells)
        set_table_border(wks.row, nrow(df1), 1, ncol(df1))
        # Write the state and adjusted/unadjusted name
        add_header(wks.row, 1, paste0(state, "- adjusted"))
        
        writeData(wb, sheet = outcome, x = df2, startCol = 1, startRow = wks.row + 12)
        # Apply the dotted inner border to the inside cells (excluding the outermost cells)
        set_table_border(wks.row+11, nrow(df2), 1, ncol(df2))
        add_header(wks.row+11, 1, paste0(state, "- unadjusted"))
        
        # Write the plots to the worksheet
        plot_path_adj <- paste0(in.path, "/combined_adj.png")
        plot_path_unadj <- paste0(in.path, "/combined_unadj.png")
        insertImage(wb, sheet=outcome, plot_path_adj, startRow=wks.row, startCol = 6, width = 8, height = 3)
        insertImage(wb, sheet=outcome, plot_path_unadj, startRow=wks.row + 11 , startCol = 6, width = 8, height = 3)
        
        set_table_border(wks.row-2, 21 , 1, 15, FALSE, FALSE)
        wks.row <- wks.row+24
      }
    }
    # # Save the Excel file
    cohort.name <- cohort.names[[cohort]]
    output_file <- paste0(outc.path, "/", cohort.name, ".xlsx")
    saveWorkbook(wb, file = output_file, overwrite = TRUE)
    # Inform the user
    cat("Excel file saved as:", output_file, "\n")
  }
}


