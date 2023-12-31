
# Data scientists use various models to analyze and extract insights from data. 
# Some of the most common models include: 
# 1. Linear Regression: A basic statistical model used to understand the relationship 
# between independent variables and a dependent variable, assuming a linear relationship. 
# 2. Logistic Regression: It’s used for binary classification problems, predicting 
# the probability of a certain event occurring. 
# 3. Decision Trees: These models use a tree-like structure to make decisions based 
# on feature values. Random Forests and Gradient Boosting are extensions of this model, 
# improving accuracy by combining multiple decision trees. 
# 4. Support Vector Machines (SVM): A supervised learning model used for classification 
# and regression analysis. SVMs aim to find the best possible boundary between classes. 
# 5. Neural Networks: Including models like Multi-layer Perceptrons (MLP), 
# Convolutional Neural Networks (CNN), and Recurrent Neural Networks (RNN), 
# these models are part of deep learning methods and are used for various tasks such as 
# image recognition, natural language processing, and time series analysis. 
# 6. Clustering Algorithms: K-means, hierarchical clustering, and DBSCAN are used to group 
# similar data points together based on certain characteristics. 
# 7. Naive Bayes: Based on Bayes’ theorem, this model is particularly useful for text 
# classification and sentiment analysis. 
# 8. Ensemble Methods: Techniques like Bagging (Bootstrap Aggregating), Boosting, and 
# Stacking combine multiple models to improve overall performance and reduce overfitting. 
# 9. Time Series Models: ARIMA (AutoRegressive Integrated Moving Average) and 
# LSTM (Long Short-Term Memory) are commonly used for analyzing and forecasting 
# time-dependent data. The choice of model depends on various factors like the nature of 
# the data, the problem being addressed (classification, regression, clustering, etc.), 
# the size of the dataset, and the desired outcome. Data scientists often experiment with 
# multiple models to find the one that best fits the particular task at hand.



# Import PCoA
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
## # SUB(2): FUNCTION TO LOAD METADATA
## ######################
import_metadata <- function(group_table){ #, group_column, sample_names){
  metadata_matrix <- as.matrix( # Load the metadata table (same if you use one or all columns)
    read.table(
      file=group_table,row.names=1,header=TRUE,sep="\t",
      colClasses = "character", check.names=FALSE,
      comment.char = "",quote="",fill=TRUE,blank.lines.skip=FALSE
    )
  )
}


library(matlab)
library(tidyverse)
library(e1071)
library(caret)
library(brm)

# SETWD
setwd("~/Documents/GitHub/PCA_tools_for_R/")

# LOAD PCoA
my_pcoa <- load_pcoa_data("HMP.Jumpstart.DESeq_normed.euclidean.PCoA")
my_pcoa_vectors <- my_pcoa$eigen_vectors
my_pcoa_vectors <- my_pcoa_vectors[sort(rownames(my_pcoa_vectors)),] # sort by mgid
my_pcoa_vectors <- as.data.frame(my_pcoa_vectors)
my_pcoa_values <- my_pcoa$eigen_values
rm(my_pcoa) # cleanup

# SEE HOW MUCH VARIATION IN THE FIRST FEW COORDINATES
round(my_pcoa_values[1:10,]*100,digits=2)
sum(round(my_pcoa_values[1:10,]*100,digits=2))
# JUST USE THE FIRST 3 COORDINATES, accounts for 58% of the variation
my_pcoa_vectors <- my_pcoa_vectors[,1:10]
# OR FIRST 10, 75% variation

# LOAD METADATA
my_metadata <- import_metadata("HMP_jumpstart_metadata.txt")
my_metadata <- my_metadata[sort(rownames(my_metadata)),] # sort by mgid
my_metadata <- as.data.frame(my_metadata)

# MAKE PCOA DATA AND METADATA TIDY -- just the env_package.data.body_site
my_tidy_data <- cbind(my_pcoa_vectors, my_metadata[,1] )
my_tidy_data <- as.data.frame(my_tidy_data)
#colnames(my_tidy_data[ncol(my_tidy_data)]) <- ""
#colnames(my_tidy_data[ncol(my_tidy_data)]) <- "env_package.data.body_site" # hmmm this doesn't work
dim(my_tidy_data)
colnames(my_tidy_data)[11] <- "env_package.data.body_site" # This works
rm(my_pcoa_vectors,my_metadata)

# SPLIT INTO TRAINING AND TEST DATA
# Create indices for splitting (70% training, 30% test) 
set.seed(123) # For reproducibility 
indices <- createDataPartition(my_tidy_data$env_package.data.body_site, p = 0.7, list = FALSE) 
# Split the dataset using the indices 
training_data <- my_tidy_data[indices, ] 
test_data <- my_tidy_data[-indices, ]
# encode metadata as factor
training_data$env_package.data.body_site <- factor(training_data$env_package.data.body_site)
test_data$env_package.data.body_site <- factor(test_data$env_package.data.body_site)
rm(my_tidy_data) #cleanup



# CREATE SVM MODEL 
# FAILS - need numeric dep variable, hot encode data as follows
#training_data$`my_metadata[, "env_package.data.body_site"]` <- as.numeric(factor(training_data$`my_metadata[, "env_package.data.body_site"]`))
svm_model <- svm( env_package.data.body_site ~ ., data = training_data)
# Predict using the SVM model
predictions <- predict(svm_model, test_data)
# check the accuracy
##accuracy <- mean(predictions == test_data$env_package.data.body_site) # This doesn't work because the levels are different between the two
##summary(predictions == test_data$env_package.data.body_site)
##print(paste("Accuracy of SVM model:", round(accuracy * 100, 2), "%"))

count_differences(predictions, test_data$env_package.data.body_site)
# 83.19% - first 3 coordinates
# 91.6 % - first 10 coordinates

# SO use this functions instead
# returns percent similiarity between two same length lists
count_differences <- function(list1, list2) {
  # Initialize a counter for differences
  num_differences <- 0
  
  # Loop through the lists and count differences
  for (i in seq_along(list1)) {
    if (list1[i] != list2[i]) {
      num_differences <- num_differences + 1
    }
  }
  
  return(100 - round(num_differences/i*100, digits=2))
}



# Try with a Baysian model 
# CREATE BAYSIAN MODEL
bayesian_model <- naiveBayes(env_package.data.body_site ~ ., data = training_data)
# Predict using the SVM model
predictions <- predict(bayesian_model, test_data)
# check the accuracy
##accuracy <- mean(predictions == test_data$env_package.data.body_site) # This doesn't work because the levels are different between the two
##summary(predictions == test_data$env_package.data.body_site)
##print(paste("Accuracy of SVM model:", round(accuracy * 100, 2), "%"))
count_differences(predictions, test_data$env_package.data.body_site)
# 75.21 % - first three coordinates
# 85.08 % - first 10 coordinates
