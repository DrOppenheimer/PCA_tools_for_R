# set the working directory
setwd("~/Documents/GitHub/PCA_tools_for_R/HMP_example")

# render interactive 3d PCoA plot from the PCoA and the corresponding metadata
source("render_calculated_pcoa.r")
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
system("open *.png")




