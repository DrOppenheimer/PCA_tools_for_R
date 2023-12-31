# Assuming your count data is stored in 'counts'
# 'counts' should be a matrix or data frame with genes as rows and samples as columns

# Calculate library size factors
libSizes <- colSums(counts)

# Calculate counts per million (CPM)
cpm <- cpm(counts)

# Calculate effective transcript length (e.g., gene length) for each gene
# Replace 'gene_lengths' with a vector or data frame containing lengths for each gene
effective_lengths <- gene_lengths

# Calculate TPM
tpm <- cpm / (effective_lengths / 1000)   # Dividing by 1000 to convert lengths to kilobases

# Normalize TPM to sum up to a constant value, e.g., 1e6
tpm_scaled <- tpm * (1e6 / colSums(tpm))

# 'tpm_scaled' now contains the TPM-normalized data