#### Set up command line: ####
args = commandArgs(trailingOnly = TRUE)

accession_ID = args[1]
virus_name = args[2]
filter_value = as.integer(args[3])
summary_file = args[4]

# Check if arguments were parsed:
if (is.na(accession_ID)){
  cat("\nERROR! No accession ID is provided", "\n")
  quit(save = 'default', status = 1, runLast = TRUE)
}
if (is.na(virus_name)){
  cat("\nERROR! No viral name is provided", "\n")
  quit(save  = 'default', status = 1, runLast = TRUE)
}
if (is.na(filter_value)){
  cat("\nWARNING! No filter value is selected, default value will be set to 100.", "\n")
  filter_value = as.integer(100)
}
if (is.na(summary_file)){
  cat("\nWARNING! No summary file is selected, default value will be set to Summary_matrix.csv\n")
  summary_file = "Summary_matrix.csv"
}



# #### Libraries: ####
library(ggplot2)
library(dplyr)
library(gridExtra)
library(devEMF)
library(devtools)
library(ComplexHeatmap)
library(circlize)



#### Pathways: ####
serratus = "../Serratus"
acc_dir = paste(serratus, accession_ID, sep="/")
plot_dir = paste(acc_dir, "Plots", sep="/")
summary_path = paste(acc_dir, summary_file, sep="/")

# Terminal output:
cat("\n")
cat("Accession ID:\t", accession_ID, "\n")
cat("Virus name:\t", virus_name, "\n")
cat("Filter value:\t", filter_value, "\n")
cat("Path to summary matrix:\t", summary_path, "\n\n")

# Check if directory and files exist:
if (!dir.exists(acc_dir)){
  cat("ERROR! Incorrect accession id has been entered: ", accession_ID, "\n")
  quit(save  = 'default', status = 1, runLast = TRUE)
}
if (!dir.exists(plot_dir)){
  dir.create(plot_dir)
}
if (!file.exists(summary_path)){
  cat("ERROR! The summary matrix does not exist in this path: ", summary_path, "\n")
  quit(save  = 'default', status = 1, runLast = TRUE)
}



#### Functions: ####
# Opens a df, with specific families (default removes selected families, if not write except False):
open_df = function(file_name, fam=NULL, except=TRUE){
  df = read.table(file_name, header=TRUE, row.names=1, sep=",")
  # Specify which families to remove/select:
  if (!is.null(fam)){
    if(isTRUE(except)){
      df = df[, !names(df) %in% fam]
    }
    else{
      df = subset(df, select = c(fam))
    }
  }
  return(df)
}

# Creates two plots and combines them:
two_plots = function(df, title1 = NULL, title2 = NULL, loging1=FALSE, loging2=TRUE, filter=100){
  # Create first table, that will sum number of samples per each family, according to the filter criteria:
  non_zero_samples = apply(df > filter, 2, sum)
  df1 = data.frame(Family = colnames(df), NonZeroSamples = non_zero_samples)
  df1 = df1[rev(order(df1$NonZeroSamples)),]
  
  # Create second table, that will show a sum of all reads that passes the filter for each family:
  total = apply(df, 2, function(x) sum(x[x > filter]))
  df2 = data.frame(Family = colnames(df), Total = total)
  df2 = df2[rev(order(df2$Total)),]
  
  # If loging1 TRUE, log by the base 10:
  if(isTRUE(loging1)){
    df1$NonZeroSamples = log10(df1$NonZeroSamples + 1)
  }
  
  # If loging2 TRUE, log by the base 10:
  if(isTRUE(loging2)){
    df2$Total = log10(df2$Total + 1)
  }
  
  # Order column names for graphs:
  column_names = df1$Family
  df1$Family = factor(df1$Family, levels = column_names)
  df2$Family = factor(df2$Family, levels = column_names)
  
  # Plot graphs:
  # Create custom theme
  my_theme = theme(
    legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic'),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.title = element_text(size = 13), axis.text = element_text(size = 11),
    plot.title = element_text(size = 16, face = 'bold'),
    axis.ticks.length = unit(0.2, 'cm')
  )
  
  plot1 = ggplot(df1, aes(x=Family, y=NonZeroSamples, fill = Family)) +
    geom_bar(stat='identity', color = "black") +
    labs(x = NULL, y = "Number of samples", title = title1) +
    theme_minimal() +
    my_theme
  
  plot2 = ggplot(df2, aes(x = Family, y = Total, fill = Family)) +
    geom_bar(stat = "identity", color = "black") +
    labs(title = title2, x = NULL, y = "Total reads (log10)") +
    theme_minimal() +
    my_theme
  
  # Combine both plots:
  plot = grid.arrange(plot1, plot2, nrow=2)
  return(plot)
}

# Parse through selected dataframe, by filtering read values, and then log-transforming them into matrix:
df_to_matrix = function(df, filter=100, loging=TRUE, delete_0s=FALSE){
  # First filter the values, if they do not pass, replace to 0:
  df[df<filter] = 0
  
  # Delete non-zero rows, if option is enabled:
  if(isTRUE(delete_0s)){
    rows = apply(df != 0, 1, any)
    columns = apply(df != 0, 2, any)
    df = df[rows, columns]
  }
  
  # Then log10 the values:
  if(isTRUE(loging)){
    df = log(df + 1, 10)
  }
  
  # Finally return the matrix:
  df = as.matrix(df)
  return(df)
}

# Saves plot (can optimise height and width):
save_plot = function(plot, name, width=10, height=8){
  name = paste(name, ".emf", sep="")
  ggsave(name, plot, width = width, height = height)
}

# Save heatmap:
save_heatmap = function(plot, name){
  name = paste(name, ".emf", sep="")
  emf(name)
  draw(plot)
  dev.off()
}



#### Script: ####
## For generating barplots:
# Open up df, remove AMR and Unknown:
acc_df = open_df(summary_path, fam=c("AMR", "Unknown"))

# Generate plot:
barplot_title1 = paste("Number of samples per family in ", virus_name, " containing samples (>", filter_value, " reads)", sep="")
barplot_title2 = paste("Total number of reads per family in ", virus_name, " containing samples (>", filter_value, " reads)", sep="")
acc_plots = two_plots(acc_df, title1 = barplot_title1, title2 = barplot_title2, filter = filter_value)

# Save plot:
plot_name = paste("Plots", virus_name, filter_value, sep="_")
plot_name = paste(plot_dir, plot_name, sep="/")
save_plot(acc_plots, plot_name)

cat("\nBarplots saved at: ", plot_name, ".emf", sep = "")

## For generating Heatmaps:
# Convert df to a matrix:
acc_matrix = df_to_matrix(acc_df, filter = filter_value)

# Create heatmap:
col_fun = colorRamp2(c(0,8), c("whitesmoke", "darkblue"))
col_fun(seq(0, 8))

heatmap_title = paste("Number of reads per family in ", virus_name, " containing samples (>", filter_value, " reads)", sep="")


acc_heatmap = Heatmap(acc_matrix, name = "log(10)",
                      column_title = heatmap_title,
                      col = col_fun, show_row_names = FALSE,
                      border_gp = gpar(col = "black", lty=1, lwd=2),
                      column_names_gp = gpar(fontface = 'italic'),
                      column_title_gp = gpar(fontface = 'bold', fontsize = 14))

# Save plot:
heatmap_name = paste("Heatmap", virus_name, filter_value, sep="_")
heatmap_name = paste(plot_dir, heatmap_name, sep="/")
save_heatmap(acc_heatmap, heatmap_name)

cat("\nHeatmap saved at: ", heatmap_name, ".emf", sep="")


