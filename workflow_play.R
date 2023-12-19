
# get the metadata subset that corresponds to the small HMP sample
# did this just once - not meant for repeating
metadata_subset <- matrix( nrow=30, ncol=76 )
colnames(metadata_subset) <- colnames(my_metadata)
rownames <- list()
for ( i in 1:ncol(my_data)){ # my data is the subset of HMP
  #print(colnames(my_data)[i]) 
  metadata_subset[i,] <- my_metadata[ colnames(my_data)[i], ]
  rownames[i] <- colnames(my_data)[i]
  }
rownames(metadata_subset) <- rownames
# export data as a text file
source("~/Documents/GitHub/Kevin_R_scripts/export_data.r")
export_data(metadata_subset, "filtered_counts.metadata.txt")
##################################################################


# Start here --------------------------------------------------------------


# set the working directory
setwd("~/Documents/GitHub/PCA_tools_for_R/")

## Load some sample data to play with
#source("~/Documents/GitHub/Kevin_R_scripts/import_data.r")
#my_data <- import_data("filtered_counts.txt")
#my_metadata <- import_data("filtered_counts.metadata.txt")


# Preprocess the data
source("~/Documents/GitHub/PCA_tools_for_R/preprocessing_tool.r")
preprocessing_tool("filtered_counts.txt")

# calculate PCoA on the preprocessed data
source("~/Documents/GitHub/PCA_tools_for_R/calculate_pco.r")
calculate_pco(file_in="filtered_counts.txt.standardize.PREPROCESSED.txt")

# render interactive 3d PCoA plot from the PCoA and the corresponding metadata
source("~/Documents/GitHub/PCA_tools_for_R/render_calculated_pcoa.r")
plot_interactive_colored_3d_pcoa(
  pcoa_data_file = "HMP.Jumpstart.DESeq_normed.euclidean.PCoA", #my_test_pcoa,
  selected_eigen_vectors = c(1,2,3),
  pcoa_metadata_file = "HMP_jumpstart_metadata.txt",
  metadata_column = "env_package.data.body_site"#metadata_colors
)

# iterate through the metadata creating a static 3d plot for each metadata column
plot_static_colored_3d_pcoas(
  pcoa_filename = "HMP.Jumpstart.DESeq_normed.euclidean.PCoA",
  selected_eigen_vectors = c(1,2,3),
  metadata_filename = "HMP_jumpstart_metadata.txt",
  debug = TRUE
)


