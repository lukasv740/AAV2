#### Libraries: ####
library(devtools)
library(ComplexHeatmap)
library(circlize)
library(dplyr)



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
save_plot = function(plot, name){
  name = paste(name, ".emf", sep="")
  emf(name)
  draw(plot)
  dev.off()
}



#### Script: ####
# Opening a summary table, and excluding families:
aav2_df = open_df(aav2_file, fam=c("AMR", "Unknown"))
aav3_df = open_df(aav3_file, fam=c("AMR", "Unknown"))
havf_df = open_df(havf_file, fam=c("AMR", "Unknown"))
asf_df = open_df(asf_file, fam=c("AMR", "Unknown"))


# Converting them into loged matrices with specified filters:
aav2_matrix_100 = df_to_matrix(aav2_df, filter=100)
aav3_matrix_100 = df_to_matrix(aav3_df, filter=100)
havf_matrix_100 = df_to_matrix(havf_df, filter=100)
asf_matrix_100 = df_to_matrix(asf_df, filter=100)

aav2_matrix_1000 = df_to_matrix(aav2_df, filter=1000)
aav3_matrix_1000 = df_to_matrix(aav3_df, filter=1000)
havf_matrix_1000 = df_to_matrix(havf_df, filter=1000)
asf_matrix_1000 = df_to_matrix(asf_df, filter=1000)


# Creating heatmaps:
col_fun = colorRamp2(c(0,8), c("whitesmoke", "darkblue"))
col_fun(seq(0, 8))

# For AAV2:
aav2_heatmap_100 = Heatmap(aav2_matrix_100, name = "log(10)",
                           column_title = "Number of reads per family in AAV2 samples (>100 reads)",
                           col = col_fun, show_row_names = FALSE,
                           border_gp = gpar(col = "black", lty=1, lwd=2),
                           column_names_gp = gpar(fontface = 'italic'),
                           column_title_gp = gpar(fontface = 'bold', fontsize = 14))
aav2_heatmap_100

aav2_heatmap_1000 = Heatmap(aav2_matrix_1000, name = "log(10)",
                            column_title = "Number of reads per family in AAV2 samples (>1000 reads)",
                            col = col_fun, show_row_names = FALSE,
                            border_gp = gpar(col = "black", lty=1, lwd=2),
                            column_names_gp = gpar(fontface = 'italic'),
                            column_title_gp = gpar(fontface = 'bold', fontsize = 14))
aav2_heatmap_1000

# For AAV3:
aav3_heatmap_100 = Heatmap(aav3_matrix_100, name = "log(10)",
                           column_title = "Number of reads per family in AAV3 samples (>100 reads)",
                           col = col_fun, show_row_names = FALSE,
                           border_gp = gpar(col = "black", lty=1, lwd=2),
                           column_names_gp = gpar(fontface = 'italic'),
                           column_title_gp = gpar(fontface = 'bold', fontsize = 14))
aav3_heatmap_100

aav3_heatmap_1000 = Heatmap(aav3_matrix_1000, name = "log(10)",
                            column_title = "Number of reads per family in AAV3 samples (>1000 reads)",
                            col = col_fun, show_row_names = FALSE,
                            border_gp = gpar(col = "black", lty=1, lwd=2),
                            column_names_gp = gpar(fontface = 'italic'),
                            column_title_gp = gpar(fontface = 'bold', fontsize = 14))
aav3_heatmap_1000

# For HAVF:
havf_heatmap_100 = Heatmap(havf_matrix_100, name = "log(10)",
                           column_title = "Number of reads per family in HAVF samples (>100 reads)",
                           col = col_fun, show_row_names = FALSE,
                           border_gp = gpar(col = "black", lty=1, lwd=2),
                           column_names_gp = gpar(fontface = 'italic'),
                           column_title_gp = gpar(fontface = 'bold', fontsize = 14))
havf_heatmap_100

havf_heatmap_1000 = Heatmap(havf_matrix_1000, name = "log(10)",
                            column_title = "Number of reads per family in HAVF samples (>1000 reads)",
                            col = col_fun, show_row_names = FALSE,
                            border_gp = gpar(col = "black", lty=1, lwd=2),
                            column_names_gp = gpar(fontface = 'italic'),
                            column_title_gp = gpar(fontface = 'bold', fontsize = 14))
havf_heatmap_1000

# For ASF:
asf_heatmap_100 = Heatmap(asf_matrix_100, name = "log(10)",
                          column_title = "Number of reads per family in ASF samples (>100 reads)",
                          col = col_fun, show_row_names = FALSE,
                          border_gp = gpar(col = "black", lty=1, lwd=2),
                          column_names_gp = gpar(fontface = 'italic'),
                          column_title_gp = gpar(fontface = 'bold', fontsize = 14))
asf_heatmap_100

asf_heatmap_1000 = Heatmap(asf_matrix_1000, name = "log(10)",
                           column_title = "Number of reads per family in ASF samples (>1000 reads)",
                           col = col_fun, show_row_names = FALSE,
                           border_gp = gpar(col = "black", lty=1, lwd=2),
                           column_names_gp = gpar(fontface = 'italic'),
                           column_title_gp = gpar(fontface = 'bold', fontsize = 14))
asf_heatmap_1000


# Save plots:
save_plot(aav2_heatmap_100, "./Plots/aav2_heatmap_100")
save_plot(aav3_heatmap_100, "./Plots/aav3_heatmap_100")
save_plot(havf_heatmap_100, "./Plots/havf_heatmap_100")
save_plot(asf_heatmap_100, "./Plots/asf_heatmap_100")

save_plot(aav2_heatmap_1000, "./Plots/aav2_heatmap_1000")
save_plot(aav3_heatmap_1000, "./Plots/aav3_heatmap_1000")
save_plot(havf_heatmap_1000, "./Plots/havf_heatmap_1000")
save_plot(asf_heatmap_1000, "./Plots/asf_heatmap_1000")