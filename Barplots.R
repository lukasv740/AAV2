#### Libraries: ####
library(ggplot2)
library(dplyr)
library(gridExtra)
library(devEMF)



#### Pathways: ####
aav2_file = "./data/Summary_matrix_aav2.csv"
aav3_file = "./data/Summary_matrix_aav3.csv"
havf_file = "./data/Summary_matrix_havf.csv"
asf_file = "./data/Summary_matrix_asf.csv"



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


# Saves plot (can optimise height and width):
save_plot = function(plot, name, width=10, height=8){
  name = paste(name, ".emf", sep="")
  ggsave(name, plot, width = width, height = height)
}



#### Script: ####
# Opening a summary matrix:
aav2_df = open_df(aav2_file, fam=c("AMR", "Unknown"))
aav3_df = open_df(aav3_file, fam=c("AMR", "Unknown"))
havf_df = open_df(havf_file, fam=c("AMR", "Unknown"))
asf_df = open_df(asf_file, fam=c("AMR", "Unknown"))


# Creating plots:
# Filter 100:
aav2_plots = two_plots(aav2_df, title1 = "Number of samples per family in AAV2 (>100 reads)",
                       title2 = "Total number of reads per family in AAV2 (>100 reads)", filter = 100)
aav3_plots = two_plots(aav3_df, title1 = "Number of samples per family in AAV3 (>100 reads)",
                       title2 = "Total number of reads per family in AAV3 (>100 reads)", filter = 100)
havf_plots = two_plots(havf_df, title1 = "Number of samples per family in HAVF (>100 reads)",
                       title2 = "Total number of reads per family in HAVF (>100 reads)", filter = 100)
asf_plots = two_plots(asf_df, title1 = "Number of samples per family in ASF (>100 reads)",
                       title2 = "Total number of reads per family in ASF (>100 reads)", filter = 100)

# Filter 1000:
aav2_plots_1000 = two_plots(aav2_df, title1 = "Number of samples per family in AAV2 (>1000 reads)",
                       title2 = "Total number of reads per family in AAV2 (>1000 reads)", filter = 1000)
aav3_plots_1000 = two_plots(aav3_df, title1 = "Number of samples per family in AAV3 (>1000 reads)",
                       title2 = "Total number of reads per family in AAV3 (>1000 reads)", filter = 1000)
havf_plots_1000 = two_plots(havf_df, title1 = "Number of samples per family in HAVF (>1000 reads)",
                       title2 = "Total number of reads per family in HAVF (>1000 reads)", filter = 1000)
asf_plots_1000 = two_plots(asf_df, title1 = "Number of samples per family in ASF (>1000 reads)",
                      title2 = "Total number of reads per family in ASF (>1000 reads)", filter = 1000)


# Save plots:
save_plot(aav2_plots, "./Plots/aav2_plots_100")
save_plot(aav3_plots, "./Plots/aav3_plots_100")
save_plot(havf_plots, "./Plots/havf_plots_100")
save_plot(asf_plots, "./Plots/asf_plots_100")

save_plot(aav2_plots_1000, "./Plots/aav2_plots_1000")
save_plot(aav3_plots_1000, "./Plots/aav3_plots_1000")
save_plot(havf_plots_1000, "./Plots/havf_plots_1000")
save_plot(asf_plots_1000, "./Plots/asf_plots_1000")
