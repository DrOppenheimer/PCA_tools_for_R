install.packages("umap")
install.packages("ggplot2")
install.packages("gg3D")
install.packages("plotly")


# Load necessary libraries
library(umap)
library(ggplot2)
library(plotly)


# Example data
set.seed(42)
data <- matrix(runif(100 * 50), nrow = 100, ncol = 50)  # 100 samples with 50 features

# Perform PCA
pca_result <- prcomp(data, scale. = TRUE)
pca_components <- pca_result$x[, 1:10]  # Keep the top 10 principal components

# Apply UMAP
umap_result <- umap(pca_components, n_components = 3)

# Convert to data frame for plotting
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2", "UMAP3")

# Plot UMAP result in 3D using plotly
plot_ly(umap_df, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3, type = "scatter3d", mode = "markers") %>%
  layout(title = "3D UMAP Projection of PCA Components",
         scene = list(xaxis = list(title = "UMAP 1"),
                      yaxis = list(title = "UMAP 2"),
                      zaxis = list(title = "UMAP 3")))

# Extract eigenvalues
eigenvalues <- pca_result$sdev^2

# Calculate the sum of the original eigenvalues
total_eigenvalues <- sum(eigenvalues)

# Scale the eigenvalues so that their sum equals 100
scaled_eigenvalues <- (eigenvalues / total_eigenvalues) * 100

# Check sum of eigenvalues
sum(scaled_eigenvalues)

# Print scaled eigenvalues
print(scaled_eigenvalues)

######################################################################################
######################################################################################
### WITH HMP DATA
library(umap)
library(plotly)

setwd("C:/Users/kevin/OneDrive/Documents/GitHub/PCA_tools_for_R/HMP_example/")

## ######################
## # SUB(1): FUNCTION TO LOAD A PRECALCULATED *.PCoA
## ######################
load_pcoa_data <- function(PCoA_in){
  # create two connections to the file
  con_1 <- file(PCoA_in)
  con_2 <- file(PCoA_in)
  # read through the first time to get the number of samples
  open(con_1);
  num_values <- 0
  data_type = "NA"
  while ( length(my_line <- readLines(con_1,n = 1, warn = FALSE)) > 0) {
    if ( length( grep("PCO", my_line) ) == 1  ){
      num_values <- num_values + 1
    }
  }
  close(con_1)
  # create object for values
  eigen_values <- matrix("", num_values, 1)
  dimnames(eigen_values)[[1]] <- 1:num_values
  # create object for vectors
  eigen_vectors <- matrix("", num_values, num_values)
  dimnames(eigen_vectors)[[1]] <- 1:num_values
  # read through the input file a second time to populate the R objects
  value_index <- 1
  vector_index <- 1
  open(con_2)
  current.line <- 1
  data_type = "NA"
  
  while ( length(my_line <- readLines(con_2,n = 1, warn = FALSE)) > 0) {
    if ( length( grep("#", my_line) ) == 1  ){
      if ( length( grep("EIGEN VALUES", my_line) ) == 1  ){
        data_type="eigen_values"
      } else if ( length( grep("EIGEN VECTORS", my_line) ) == 1 ){
        data_type="eigen_vectors"
      }
    }else{
      split_line <- noquote(strsplit(my_line, split="\t"))
      if ( identical(data_type, "eigen_values")==TRUE ){
        dimnames(eigen_values)[[1]][value_index] <- noquote(split_line[[1]][1])
        eigen_values[value_index,1] <- noquote(split_line[[1]][2])       
        value_index <- value_index + 1
      }
      if ( identical(data_type, "eigen_vectors")==TRUE ){
        dimnames(eigen_vectors)[[1]][vector_index] <- noquote(split_line[[1]][1])
        for (i in 2:(num_values+1)){
          eigen_vectors[vector_index, (i-1)] <- as.numeric(noquote(split_line[[1]][i]))
        }
        vector_index <- vector_index + 1
      }
    }
  }
  close(con_2)
  
  # finish labeling of data objects and return them in a single list object
  dimnames(eigen_values)[[2]] <- "EigenValues"
  dimnames(eigen_vectors)[[2]] <- dimnames(eigen_values)[[1]]
  class(eigen_values) <- "numeric"
  class(eigen_vectors) <- "numeric"
  return(list(eigen_values=eigen_values, eigen_vectors=eigen_vectors))
  
}
## ######################
## ######################

## ######################
## # IMPORT PRECALCULATED *.PCoA
## ######################
my_pcoa <- load_pcoa_data("HMP.Jumpstart.DESeq_normed.euclidean.PCoA")

## Example data
#set.seed(42)
#data <- matrix(runif(100 * 50), nrow = 100, ncol = 50)  # 100 samples with 50 features

# Perform PCA
#pca_result <- prcomp(data, scale. = TRUE)
pcoa_components <- my_pcoa$eigen_vectors[, 1:10]  # Keep the top 10 principal components

# Extract eigenvalues and compute variance explained
#eigenvalues <- pca_result$sdev^2

total_eigenvalues <- sum(my_pcoa$eigen_values)
variance_explained <- (my_pcoa$eigen_values / total_eigenvalues) * 100
variance_first_ten <- sum(variance_explained[1:10]) 

# Apply UMAP
umap_result <- umap(pcoa_components, n_components = 3)

# Convert to data frame for plotting
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2", "UMAP3")



## ######################
## # SUB(2): FUNCTION TO LOAD METADATA FILE
## ######################
load_metadata <- function(metadata_file){
  metadata_matrix <- as.matrix( # Load the metadata table (same if you use one or all columns)
    read.table(
      file=metadata_file,row.names=1,header=TRUE,sep="\t",
      colClasses = "character", check.names=FALSE,
      comment.char = "",quote="",fill=TRUE,blank.lines.skip=FALSE
    )
  ) 
  return(metadata_matrix)
}
## ######################
## ######################

# Import metadata
my_metadata <- load_metadata("HMP_jumpstart_metadata.txt")

# select a single metadata column
metadata_column <- as.matrix(my_metadata[,"env_package.data.body_site"]) 


## ######################
## # SUB(3): FUNCTION TO GENERATE COLORS FOR THE PLOT 
## ######################
create_colors <- function(metadata_column, color_mode = "auto"){ # function to     
  my_data.color <- data.frame(metadata_column)
  column_factors <- as.factor(metadata_column[,1])
  column_levels <- levels(as.factor(metadata_column[,1]))
  num_levels <- length(column_levels)
  color_levels <- col.wheel(num_levels)
  levels(column_factors) <- color_levels
  my_data.color[,1]<-as.character(column_factors)
  
  my_colors_legText_legColor <- list(colors = my_data.color, legText = column_levels, legColors = color_levels)
  # combined_object <- list(data_frame = df, list1 = char_list1, list2 = char_list2)
  return(my_colors_legText_legColor)
}

######################
# SUB(8): Create optimal contrast color selection using a color wheel
# adapted from https://stat.ethz.ch/pipermail/r-help/2002-May/022037.html 
######################
col.wheel <- function(num_col, my_cex=0.75) {
  cols <- rainbow(num_col)
  col_names <- vector(mode="list", length=num_col)
  for (i in 1:num_col){
    col_names[i] <- getColorTable(cols[i])
  }
  cols
}
######################
######################

######################
# SUB(9): The inverse function to col2rgb()
# adapted from https://stat.ethz.ch/pipermail/r-help/2002-May/022037.html
######################
rgb2col <- function(rgb) {
  rgb <- as.integer(rgb)
  class(rgb) <- "hexmode"
  rgb <- as.character(rgb)
  rgb <- matrix(rgb, nrow=3)
  paste("#", apply(rgb, MARGIN=2, FUN=paste, collapse=""), sep="")
}
######################
######################

######################
# SUB(10): Convert all colors into format "#rrggbb"
# adapted from https://stat.ethz.ch/pipermail/r-help/2002-May/022037.html
######################
getColorTable <- function(col) {
  rgb <- col2rgb(col);
  col <- rgb2col(rgb);
  sort(unique(col))
}
######################
######################


# Get colors for selected column
my_colors <- create_colors(metadata_column)
#return(my_data.color, column_levels, color_levels)



########################################
########################################
### Plot without legend
# Plot UMAP result in 3D using plotly
plot <- plot_ly(umap_df, 
                x = ~UMAP1, 
                y = ~UMAP2, 
                z = ~UMAP3, 
                type = "scatter3d", 
                mode = "markers",
                marker = list(color = my_colors$colors[,1])) %>%
  layout(title = paste0("3D UMAP Projection of PCoA Components (Variance Explained: ",
                        round(variance_first_ten, 2), "%)"),
         scene = list(xaxis = list(title = "UMAP 1",
                      yaxis = list(title = "UMAP 2",
                      zaxis = list(title = "UMAP 3")))))
         #annotations = list(text=my_colors$legText),
         #showlegend = FALSE
         
plot


########################################
########################################
### Plot with legend play


# create a data frame with a category axis

my_df <- cbind(metadata_column, umap_df)

# Define custom legend text and colors
#legend_labels <- c("Alpha", "Beta", "Gamma")  # Custom legend text
#legend_colors <- c("red", "blue", "green")    # Custom legend colors
categories <- unique(my_df$metadata_column)  # Extract unique categories

# Create an empty plot
fig <- plot_ly()
# Add each category as a separate trace with custom legend name and color
for (i in seq_along(categories)) {
  cat_data <- my_df[my_df$metadata_column == categories[i], ]  # Filter category data
  fig <- fig %>%
    add_trace(
      x = cat_data$UMAP1, 
      y = cat_data$UMAP2, 
      z = cat_data$UMAP3, 
      type = "scatter3d", 
      mode = "markers",
      marker = list(color = my_colors$legColors[i], size = 6),  # Custom color
      name = legend_labels[i]  # Custom legend text
    )
}

# Customize layout
fig <- fig %>%
  layout(
    title = "Custom Legend in 3D Plot",
    scene = list(xaxis = list(title = "UMAP1"), 
                 yaxis = list(title = "UMAP2"), 
                 zaxis = list(title = "UMAP3")),
    legend = list(title = list(text = "Test Legend"))
  )

fig







# Define custom legend text and colors
legend_labels <- c("Alpha", "Beta", "Gamma")  # Custom legend text
legend_colors <- c("red", "blue", "green")    # Custom legend colors
categories <- unique(df$category)  # Extract unique categories

# Create an empty plot
fig <- plot_ly()

# Add each category as a separate trace with custom legend name and color
for (i in seq_along(categories)) {
  cat_data <- df[df$category == categories[i], ]  # Filter category data
  fig <- fig %>%
    add_trace(
      x = cat_data$x, 
      y = cat_data$y, 
      z = cat_data$z, 
      type = "scatter3d", 
      mode = "markers",
      marker = list(color = legend_colors[i], size = 6),  # Custom color
      name = legend_labels[i]  # Custom legend text
    )
}

# Customize layout
fig <- fig %>%
  layout(
    title = "Custom Legend in 3D Plot",
    scene = list(xaxis = list(title = "X"), 
                 yaxis = list(title = "Y"), 
                 zaxis = list(title = "Z")),
    legend = list(title = list(text = "Legend Title"))
  )

fig







########################################
########################################

library(plotly)

# Sample data
df <- data.frame(
  category = c("A", "B", "C", "A", "B", "C"),
  x = c(1, 2, 3, 4, 5, 6),
  y = c(3, 2, 4, 5, 6, 7),
  z = c(7, 8, 5, 4, 3, 6)
)

# Define custom legend text and colors
legend_labels <- c("Alpha", "Beta", "Gamma")  # Custom legend text
legend_colors <- c("red", "blue", "green")    # Custom legend colors
categories <- unique(df$category)  # Extract unique categories

# Create an empty plot
fig <- plot_ly()

# Add each category as a separate trace with custom legend name and color
for (i in seq_along(categories)) {
  cat_data <- df[df$category == categories[i], ]  # Filter category data
  fig <- fig %>%
    add_trace(
      x = cat_data$x, 
      y = cat_data$y, 
      z = cat_data$z, 
      type = "scatter3d", 
      mode = "markers",
      marker = list(color = legend_colors[i], size = 6),  # Custom color
      name = legend_labels[i]  # Custom legend text
    )
}

# Customize layout
fig <- fig %>%
  layout(
    title = "Custom Legend in 3D Plot",
    scene = list(xaxis = list(title = "X"), 
                 yaxis = list(title = "Y"), 
                 zaxis = list(title = "Z")),
    legend = list(title = list(text = "Legend Title"))
  )

fig
















########################################
########################################


kaleido(plot, file="test_plot.png")



## Create legend from legend data

# Load the grid package
library(grid)

# Define the legend names and colors
legend_names <- c(my_colors$legText)
legend_colors <- c(my_colors$legColors)

# Create a function to draw the custom legend
create_legend <- function(names, colors) {
  # Set up the viewport
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(length(names), 2, widths = unit(c(1, 4), "lines"))))
  
  # Draw the colored points and text
  for (i in seq_along(names)) {
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
    grid.rect(gp = gpar(fill = colors[i]), width = unit(1, "lines"), height = unit(1, "lines"))
    upViewport()
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 2))
    grid.text(names[i], x = unit(0.1, "npc"), just = "left")
    upViewport()
  }
}

# Create and display the custom legend
create_legend(legend_names, legend_colors)












library(gridExtra)

# Combine the custom legend with other plots (if needed)
grid.arrange(legend_plot)
















######################################################################################
######################################################################################
######################################################################################
### JUNK BELOW HERE


# Plot with amount of variance in each UMAP coordinate explained 
# NOPE - just shows the first 3 sorted eigenvalues as if they were
# the variation in the UMAP projection.
# Load necessary libraries
library(umap)
library(plotly)

# Example data
set.seed(42)
data <- matrix(runif(100 * 50), nrow = 100, ncol = 50)  # 100 samples with 50 features

# Perform PCA
pca_result <- prcomp(data, scale. = TRUE)
pca_components <- pca_result$x[, 1:10]  # Keep the top 10 principal components

# Extract eigenvalues and compute variance explained
eigenvalues <- pca_result$sdev^2
total_eigenvalues <- sum(eigenvalues)
variance_explained <- (eigenvalues / total_eigenvalues) * 100

# Apply UMAP
umap_result <- umap(pca_components, n_components = 3)

# Convert to data frame for plotting
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2", "UMAP3")

# Plot UMAP result in 3D using plotly
plot <- plot_ly(umap_df, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3, type = "scatter3d", mode = "markers") %>%
  layout(title = paste0("3D UMAP Projection of PCA Components (Variance Explained: ",
                        round(variance_explained[1], 2), "%, ", 
                        round(variance_explained[2], 2), "%, ", 
                        round(variance_explained[3], 2), "%)"),
         scene = list(xaxis = list(title = paste0("UMAP 1 (", round(variance_explained[1], 2), "%)")),
                      yaxis = list(title = paste0("UMAP 2 (", round(variance_explained[2], 2), "%)")),
                      zaxis = list(title = paste0("UMAP 3 (", round(variance_explained[3], 2), "%)"))))

plot
# This was from Copilot -- not correct, just shows the eigenvalues as from PCA as if
# They are variation from UMAP





