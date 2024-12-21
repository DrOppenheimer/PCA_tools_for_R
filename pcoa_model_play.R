
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


# Function to compare two lists
# returns percent similarity between two same length lists
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


# START HERE --------------------------------------------------------------

library(matlab)
library(tidyverse)
library(e1071)
library(caret)
library(brm)
library(neuralnet)

# SETWD
#setwd("~/Documents/GitHub/PCA_tools_for_R/")
setwd("~/GitHub/PCA_tools_for_R/")

# LOAD PCoA
my_pcoa <- load_pcoa_data("HMP.Jumpstart.DESeq_normed.euclidean.PCoA")
my_pcoa_vectors <- my_pcoa$eigen_vectors
my_pcoa_vectors <- my_pcoa_vectors[sort(rownames(my_pcoa_vectors)),] # sort by mgid
my_pcoa_vectors <- as.data.frame(my_pcoa_vectors)
my_pcoa_values <- my_pcoa$eigen_values
rm(my_pcoa) # cleanup

# Choose number of coordinates to use below
num_coord <- 15

# SEE HOW MUCH VARIATION IN THE FIRST FEW COORDINATES
round(my_pcoa_values[1:num_coord,]*100,digits=2)
sum(round(my_pcoa_values[1:num_coord,]*100,digits=2))
# JUST USE THE FIRST SELECTED COORDINATES
my_pcoa_vectors <- my_pcoa_vectors[,1:num_coord]
# FIRST  3  58% variation
# FIRST  4, 62% variation
# FIRST 10, 75% variation
# FIRST 20, 80% variation
# LOAD METADATA
my_metadata <- import_metadata("HMP_jumpstart_metadata.txt")
my_metadata <- my_metadata[sort(rownames(my_metadata)),] # sort by mgid
my_metadata <- as.data.frame(my_metadata)

# MAKE PCOA DATA AND METADATA TIDY -- just the env_package.data.body_site
# my_tidy_data <- cbind(my_pcoa_vectors, my_metadata[,"env_package.data.body_site"] ) # won't work, creates colname "my_metadata[, \"env_package.data.body_site\"]"
# Do it this way
my_tidy_data <- cbind(my_pcoa_vectors, env_package.data.body_site=my_metadata[,"env_package.data.body_site"] )
my_tidy_data <- as.data.frame(my_tidy_data)
colnames(my_tidy_data) # This works
##colnames(my_tidy_data[ncol(my_tidy_data)]) <- "" # misplaced ")"
##colnames(my_tidy_data[ncol(my_tidy_data)]) <- "env_package.data.body_site" # hmmm this doesn't work # misplaced ")"
dim(my_tidy_data)
###colnames(my_tidy_data)[11] <- "env_package.data.body_site" # This works
##colnames(my_tidy_data)[ncol(my_tidy_data)] <- "env_package.data.body_site" # This works

# clean house
# rm(my_pcoa_vectors,my_metadata)

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
# rm(my_tidy_data) #cleanup

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
# 88.24% - first 4 coordinates
# 91.6% - first 10 coordinates
# 92.86% - first 20 coordinates

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
# 80.67 % - first four coordinates
# 85.08 % - first 10 coordinates
# 86.13 % - first 20 coordinates 

# https://www.datacamp.com/tutorial/neural-network-models-r
# from # https://medium.com/geekculture/introduction-to-neural-network-2f8b8221fbd3#:~:text=Number%20of%20Neurons%20and%20Number%20of%20Layers%20in%20Hidden%20Layer&text=The%20number%20of%20hidden%20neurons,size%20of%20the%20output%20layer.
# The number of hidden neurons should be between the size of the input layer and the size of the output layer.
# The number of hidden neurons should be 2/3 the size of the input layer, plus the size of the output layer.
# The number of hidden neurons should be less than twice the size of the input layer.
num_neurons <- ceiling( (2/3*num_coord) + length(unique(training_data[,"env_package.data.body_site"])) )
num_neurons
# Try Neural Network model -- with 3 and 10 coordinates predicted the same for everything
# then did 20 coordinates with 30 neurons got 92% accuracy
# then 20 coordinates with two layers, 10 in each (10,10) got 91% accuracy
# "                   with one layer, 10 neurons got 92% accuracy
# then 10 coordinates with one layer, 10 or 23 neurons got 25%
# then 15 coordinates with one layer, 10 neurons got 25%
# then 15 coordinates with one layer, 26 neurons got 91% (this and all above with sigmoid - default)
# then 15 coordinates with one layer, 26 neurons got __ (with tanh) does not converge
# then 15 coordinates with one layer, 26 neurons got __ (with softmax)

NN_model = neuralnet(
  env_package.data.body_site~.,
  data=training_data,
  #act.fct = sigmoid, # 'logistic' by default (same as sigmoid)
  hidden=c(26), # 4,2 # for two hidden layers # 
  linear.output = FALSE,  # For classification tasks, set to FALSE
  threshold = 0.01  # Adjust the threshold for classification decisions
  )
# view model
plot(NN_model,rep = "best")
# Make predictions on the test set
raw_predictions <- predict(NN_model, newdata = test_data, type = "response")
# Convert raw predictions to categorical labels
predicted_labels <- levels(training_data$env_package.data.body_site)[max.col(raw_predictions, "first")]
# Convert predicted_labels to a factor with the same levels as the training set
predicted_labels <- factor(predicted_labels, levels = levels(training_data$env_package.data.body_site))
# Stat of the results
count_differences(predicted_labels, test_data$env_package.data.body_site)
# 92.39% with first 20 coordinates and 30 neurons
# Combine the original test_data with the predicted labels
results <- cbind(test_data, PredictedSpecies = predicted_labels)
# Display the results
print(results)

# other activation functions (from https://rpubs.com/shailesh/activation-functions )
# https://www.v7labs.com/blog/neural-networks-activation-functions#:~:text=drive%20V7's%20tools.-,What%20is%20a%20Neural%20Network%20Activation%20Function%3F,prediction%20using%20simpler%20mathematical%20operations.
softplus <- function(x) log(1+exp(x))
softmax <- function(x) exp(x)/sum(exp(x))
tan(x) # built in
atan(x) # built in
sin(x) # built in
#sigmoid <- function(x) 1/(1 + exp(-x)). # 1/(1 + e^-x)) # already built into R, same as "Logistic"
softsign <- function(x) x / (abs(x) + 1)
swish <- function(x) x * sigmoid(x)
binary_step <- function(x) ifelse(x >= 0, 1, 0)
rectified_linear_unit <- function(x) ifelse(x < 0 , 0, x )
leaky_rectified_linear_unit <- function(x) ifelse(x < 0 , 0.01 *x , x )
bent_identity <- function(x) (sqrt(x^2 + 1) - 1)/2 + x
sinc <- function(x) ifelse(x == 0, 1, sin(x) / x)
guassian <- function(x) exp(-x^2)

# Fun
plot(x=-10:10, y=sigmoid(-10:10))
plot(x=-10:10, y=softmax(-10:10)) # <- for multiple classification
plot(x=-10:10, y=softplus(-10:10)) # not scaled to 1 
plot(x=-10:10, y=tanh(-10:10))
plot(x=-10:10, y=softsign(-10:10))
plot(x=-10:10, y=swish(-10:10))
plot(x=-10:10, y=binary_step(-10:10))
plot(x=-10:10, y=rectified_linear_unit(-10:10))
plot(x=-10:10, y=leaky_rectified_linear_unit(-10:10))
plot(x=-10:10, y=bent_identity(-10:10))
plot(x=-10:10, y=sinc(-10:10))
plot(x=-10:10, y=guassian(-10:10))

#check 
net$act.fct(x)
# from https://stackoverflow.com/questions/47818075/how-to-use-custom-activation-function-in-neuralnet-in-r




# Convolutional Neural Network in R with Keras











