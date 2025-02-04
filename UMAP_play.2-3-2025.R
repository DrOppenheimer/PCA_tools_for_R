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





