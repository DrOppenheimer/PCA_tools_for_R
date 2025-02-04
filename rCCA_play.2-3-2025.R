
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("mixOmics")

# Here’s a complete example code that performs PCA on two datasets 
# and then applies rCCA to the principal components:

# Load necessary library
library(mixOmics)

# Example data (two datasets)
set.seed(42)
data1 <- matrix(runif(100 * 50), nrow = 100, ncol = 50)  # Dataset 1: 100 samples with 50 features
data2 <- matrix(runif(100 * 60), nrow = 100, ncol = 60)  # Dataset 2: 100 samples with 60 features

# Perform PCA on both datasets
pca_result1 <- prcomp(data1, scale. = TRUE)
pca_result2 <- prcomp(data2, scale. = TRUE)

# Extract the top 10 principal components from both datasets
pca_components1 <- pca_result1$x[, 1:10]
pca_components2 <- pca_result2$x[, 1:10]

# Perform rCCA on the principal components
rcca_result <- rcc(X = pca_components1, Y = pca_components2, lambda1 = 0.1, lambda2 = 0.1)

# Display the canonical correlations
print(rcca_result$cor)

# Define the grouping of samples
group <- rep(1:2, each = nrow(pca_components1) / 2)  # Adjust according to the number of samples

# Plot the canonical variables
plotIndiv(rcca_result, comp = 1:2, group = group, legend = TRUE)

# In this code:

# We generate two example datasets data1 and data2.
# Perform PCA on both datasets and extract the top 10 principal components.
# Apply rCCA using the rcc function from the mixOmics package.
# Display the canonical correlations.
# Plot the canonical variables.
# You can adjust the lambda1 and lambda2 parameters to control the 
# regularization for each dataset.