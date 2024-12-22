
# INSTALLATION ------------------------------------------------------------

# On the console
# Install Python 

# Install pip3

# Install keras version 2.12.0
pip3 install keras==2.12.0
# Install tensorflow version 2.12.0
pip3 install tensorflow==2.12.0

#devtools::install_github('rstudio/tensorflow')
install.packages("reticulate")
devtools::install_github("rstudio/keras")

# which -a python python3
library(reticulate)
use_python("/usr/bin/python3")

library(keras)
is_keras_available()


# TEST --------------------------------------------------------------------

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


# Function to scale values in a matrix between 0 and 1
scale_matrix <- function(mat) {
  scaled <- (mat - min(mat)) / (max(mat) - min(mat))
  return(scaled)
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


# Function to replace caret::createDataPartition
createDataPartition <- function(y, times = 1, p = 0.5, list = TRUE, groups = min(5, length(y))) {
  if (class(y)[1] == "factor") {
    y <- as.integer(y)
  }
  
  # Calculate the number of samples in each group
  cuts <- unique(quantile(y, probs = seq(0, 1, length = groups + 1), na.rm = TRUE, names = FALSE))
  y <- cut(y, cuts, include.lowest = TRUE)
  
  # Partition the data
  out <- vector(mode = "list", length = times)
  for (j in 1:times) {
    for (i in levels(y)) {
      inTrain <- sample(which(y == i), size = ceiling(p * length(which(y == i))))
      out[[j]] <- c(out[[j]], inTrain)
    }
  }
  
  if (!list) {
    out <- do.call("cbind", out)
  }
  
  return(out)
}












# which -a python python3
library(reticulate)
use_python("/usr/bin/python3")

library(keras)
is_keras_available()

library(tidyverse)
library(caret)

# SETWD
setwd("~/Documents/GitHub/PCA_tools_for_R/")


# LOAD PCoA
my_pcoa <- load_pcoa_data("HMP.Jumpstart.DESeq_normed.euclidean.PCoA")
my_pcoa_vectors <- my_pcoa$eigen_vectors
my_pcoa_vectors <- my_pcoa_vectors[sort(rownames(my_pcoa_vectors)),] # sort by mgid
my_pcoa_vectors <- as.data.frame(my_pcoa_vectors)
my_pcoa_values <- my_pcoa$eigen_values
#rm(my_pcoa) # cleanup

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
# FIRST 15, 78% variation
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

# SPLIT INTO TRAINING AND TEST DATA
# Create indices for splitting (70% training, 30% test)
set.seed(123) # For reproducibility
# mix the data by rows
my_tidy_data <- my_tidy_data[sample(nrow(my_tidy_data)), ]
# create index
indices <- createDataPartition(my_tidy_data$env_package.data.body_site, p = 0.7, list = FALSE)
# Split the dataset using the indices
training_data <- my_tidy_data[indices, ]
testing_data     <- my_tidy_data[-indices, ]
# encode metadata as factor
training_data$env_package.data.body_site <- factor(training_data$env_package.data.body_site)
 testing_data$env_package.data.body_site <- factor( testing_data$env_package.data.body_site)

# Prepare the data
train_labels <- to_categorical(as.numeric(training_data$env_package.data.body_site) - 1)
train_data <- as.matrix(training_data %>% select(-env_package.data.body_site))
rownames(train_data) <- 1:nrow(train_data) # Keras doesn't like rownames

test_labels <-  to_categorical(as.numeric(testing_data$env_package.data.body_site) - 1)
test_data <- as.matrix(testing_data %>% select(-env_package.data.body_site))
rownames(test_data) <- 1:nrow(test_data)

# Don't need to split manually, Keras can do it automatically with the "validation_split" arg under fit
my_tidy_data$env_package.data.body_site <- factor(my_tidy_data$env_package.data.body_site)
all_labels <- to_categorical(as.numeric(my_tidy_data$env_package.data.body_site) - 1)
all_data <- as.matrix(my_tidy_data %>% select(-env_package.data.body_site))
rownames(all_data) <- 1:nrow(all_data)

# try scaling train_data values between 0 and 1
all_data <- scale_matrix(all_data)

# Figure out the number of neurons in each layer
# how big to make the input layer # https://stats.stackexchange.com/questions/181/how-to-choose-the-number-of-hidden-layers-and-nodes-in-a-feedforward-neural-netw
num_input_neurons <- ncol(all_data) - 1
# How big to make the hidden layer# # from # https://medium.com/geekculture/introduction-to-neural-network-2f8b8221fbd3#:~:text=Number%20of%20Neurons%20and%20Number%20of%20Layers%20in%20Hidden%20Layer&text=The%20number%20of%20hidden%20neurons,size%20of%20the%20output%20layer.
num_hidden_neurons <- ceiling( (2/3*num_coord) + length(unique(my_metadata[,"env_package.data.body_site"])) )
# How big to make the output layer
num_output_neurons <- length(unique(my_metadata[,"env_package.data.body_site"]))

# Define your model 
model <- keras_model_sequential() %>% 
  layer_dense(units = num_input_neurons, activation = 'relu', input_shape = ncol(all_data)) %>% 
  layer_dropout(rate = 0.3) %>% 
  layer_dense(units = num_hidden_neurons, activation = 'relu') %>% #num_hidden_neurons
  layer_dropout(rate = 0.3) %>% 
  layer_dense(units = num_output_neurons, activation = 'softmax') #num_output_neurons

# Compile the model
model %>% compile(
  loss = "categorical_crossentropy",
  optimizer = optimizer_rmsprop(),
  metrics = c("accuracy")
)

# Train the model 
#history <- model %>% fit( train_data, train_labels, epochs = 1000 )
history <- model %>% fit( all_data, all_labels, epochs = 5000, batch_size = 128, validation_split = 0.2 ) 
#history <- model %>% fit( X, y, epochs = 20, validation_split = 0.2 ) 
# plot(history)

# Evaluate the model 
evaluation <- model %>% evaluate(all_data, all_labels) 
print(evaluation)



test_data <- scale_matrix(test_data)


# predict classes
prediction <- model %>% predict(test_data) %>% k_argmax()





test <- to_categorical(as.numeric(prediction))


# Compare the matrices element-wise
comparison <- test == test_labels

# View the comparison result
print(comparison)

# Count the number of "FALSE" values in the matrix
false_count <- sum(comparison == FALSE)

# Print the count
print( 1 - (false_count/2)/nrow(comparison) )








# Entire model
# The entire model can be saved to a file that contains the weight values, 
# the model’s configuration, and even the optimizer’s configuration. This allows 
# you to checkpoint a model and resume training later —from the exact same state 
# —without access to the original code.

# Save entire model to the SavedModel format
model %>% save_model_tf('my_model/')

# Recreate the exact same model, including weights and optimizer.
model <- load_model_tf('my_model/')
###################################################################
###################################################################
X <- as.matrix(training_data %>% select(-env_package.data.body_site))
y <- to_categorical(as.numeric(training_data$env_package.data.body_site) - 1)

rownames(X) <- 1:nrow(X)

# Define your model 
model <- keras_model_sequential() %>% 
  layer_dense(units = 16, activation = 'relu', input_shape = ncol(X)) %>% 
  layer_dropout(rate = 0.1) %>% 
  layer_dense(units = 26, activation = 'relu') %>% 
  layer_dropout(rate = 0.1) %>% 
  layer_dense(units = num_output_neurons, activation = 'softmax') 

# Compile the model 
#model %>% compile( loss = 'binary_crossentropy', optimizer = optimizer_rmsprop(), metrics = c('accuracy') ) 

# Compile the model
model %>% compile(
  loss = "categorical_crossentropy",
  optimizer = optimizer_rmsprop(),
  metrics = c("accuracy")
)



# Train the model 
history <- model %>% fit( X, y, epochs = 10, batch_size = 128, validation_split = 0.2 ) 
history <- model %>% fit( X, y, epochs = 20, validation_split = 0.2 ) 


# Evaluate the model 
evaluation <- model %>% evaluate(X_test, y_test) 
print(evaluation)

###################################################################
###################################################################




# Load required libraries
library(keras)
library(dplyr)

# Load the Iris dataset
data(iris)

# Shuffle the dataset
set.seed(123)
iris <- iris[sample(nrow(iris)), ]

# Split data into training and test sets
train_index <- sample(nrow(iris), 0.7 * nrow(iris))
train_data <- iris[train_index, ]
test_data <- iris[-train_index, ]

# Preprocess the data
x_train <- as.matrix(select(train_data, -Species))  # Features
y_train <- to_categorical(as.numeric(train_data$Species) - 1)  # Labels

x_test <- as.matrix(select(test_data, -Species))    # Features
y_test <- to_categorical(as.numeric(test_data$Species) - 1)    # Labels

# Define the model architecture
model <- keras_model_sequential() %>%
  layer_dense(units = 16, activation = "relu", input_shape = ncol(x_train)) %>%
  layer_dense(units = 3, activation = "softmax")

# Compile the model
model %>% compile(
  loss = "categorical_crossentropy",
  optimizer = optimizer_rmsprop(),
  metrics = c("accuracy")
)

# Train the model
history <- model %>% fit(
  x_train, y_train,
  epochs = 100,
  batch_size = 32,
  validation_split = 0.2
)

# Evaluate the model on test data
score <- model %>% evaluate(x_test, y_test)
cat("Test loss:", score["loss"], "\n")
cat("Test accuracy:", score["accuracy"], "\n")












# 
# model %>%
#   
#   # Adds a densely-connected layer with 64 units to the model:
#   layer_dense(units = num_input_neurons, activation = 'relu') %>%
#   
#   # Add another:
#   layer_dense(units = num_hidden_neurons, activation = 'relu') %>%
#   
#   # Add a softmax layer with 10 output units:
#   layer_dense(units = num_output_neurons, activation = 'softmax')


# Configure a model for categorical classification.
model %>% compile(
  optimizer = optimizer_rmsprop(learning_rate = 0.01),
  loss = "categorical_crossentropy",
  metrics = list("categorical_accuracy")
)


data <- keras_training_data <- as.matrix(training_data[,-ncol(training_data)])
labels <- keras_training_labels <- as.numeric(training_data[,ncol(training_data)])/max(as.numeric(training_data[,ncol(training_data)]))

data <- data[,1]

model <- keras_model_sequential()

model %>%
  
  # Adds a densely-connected layer with 64 units to the model:
  layer_dense(units = 15, activation = 'relu') %>%
  
  # Add another:
  layer_dense(units = 15, activation = 'relu') %>%
  
  # Add a softmax layer with 10 output units:
  layer_dense(units = 16, activation = 'softmax')













inputs <- layer_input(shape = (num_input_neurons))  # Returns a placeholder tensor
# layer_dense(units = 256, activation = 'relu', input_shape = c(784))
predictions <- inputs %>%
  layer_dense(units = num_hidden_neurons, activation = 'relu') %>%
  layer_dense(units = num_hidden_neurons, activation = 'relu') %>%
  layer_dense(units = num_output_neurons, activation = 'softmax')

# Instantiate the model given inputs and outputs.
model <- keras_model(inputs = inputs, outputs = predictions)


# Configure a model for categorical classification.
model %>% compile(
  optimizer = optimizer_rmsprop(learning_rate = 0.01),
  loss = "categorical_crossentropy",
  metrics = list("categorical_accuracy")
)







# Trains for 5 epochs
model %>% fit(
  data,
  labels,
  #batch_size = 32,
  epochs = 10
)
