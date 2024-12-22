# DOCKER
# Install docker
# Install wsl first
# In powershell as admin: 
# wsl --install
# Then run the docker installer from here:
# https://docs.docker.com/desktop/setup/install/windows-install/
# Run docker image
# start docker desktop
# then in powershell: 
# docker run -d -e USER=rstudio -v C:/Users/kosso/Downloads/:/home/rstudio -e PASSWORD=test -p 8787:8787 --name rstudio jamestripp/rkerasrstudio
docker_key <- system("docker run -d -e USER=rstudio -v C:/Users/kosso/Downloads/:/home/rstudio -e PASSWORD=test -p 8787:8787 --name rstudio jamestripp/rkerasrstudio", intern = TRUE) 
system(paste("docker stop", docker_key)) # to stop the docker
system(paste("docker remove", docker_key)) # to remove the docker
# then in a web browser:
# http://localhost:8787
# system("start duckduckgo http://localhost:8787") # This failed with a 127
shell("start duckduckgo http://localhost:8787") # this worked
##################################
##################################
# Save docker locally
docker save -o [filename].tar [image_name]
# Run docker from local image
docker load -i [filename].tar
docker images
docker run -it [image_name] /bin/bash
##################################
##################################

# LOCAL
# From here: https://towardsdatascience.com/installing-keras-tensorflow-using-anaconda-for-machine-learning-44ab28ff39cb
# Install Anoconda: https://www.anaconda.com/download-success
# Do the "Just me" install
# Start the Anoconda prompt and do the following
conda create --name PythonCPU
activate PythonCPU
conda install python=3.6
conda install -c anaconda keras
conda install tensorflow
conda install spyder
conda install -c anaconda pandas
conda install -c anaconda xlrd
conda install -c anaconda xlwt
conda install -c anaconda seaborn
conda install -c anaconda scikit-learn
conda install pillow
spyder
# Then in spyder run the following
import numpy as np # For numerical fast numerical calculations
import matplotlib.pyplot as plt # For making plots
import pandas as pd # Deals with data
import seaborn as sns # Makes beautiful plots
from sklearn.preprocessing import StandardScaler # Testing sklearn
import tensorflow # Imports tensorflow
import keras # Imports keras
# in R run
library(keras)
library(reticulate)
# in case you run into error run this : reticulate::py_discover_config("keras") 
# use_python("<yourpath>/Anaconda3/envs/r-tensorflow/Scripts/python.exe")
use_python("C:/Users/kosso/AppData/Local/r-miniconda/envs/PythonCPU")
# change <yourpath> approriately( "conda env list" to see all envs )
# write all the codes for building model in keras (or tensorflow) e.g. mnist<-dataset_mnist()


#############################
#############################

# MODEL TEST

# Prepare the data
mnist <- dataset_mnist()
x_train <- mnist$train$x
y_train <- mnist$train$y
x_test <- mnist$test$x
y_test <- mnist$test$y

x_train <- array_reshape(x_train, c(nrow(x_train), 784)) / 255
x_test <- array_reshape(x_test, c(nrow(x_test), 784)) / 255

y_train <- to_categorical(y_train, 10)
y_test <- to_categorical(y_test, 10)

# Build the model
model <- keras_model_sequential() %>%
  layer_dense(units = 256, activation = 'relu', input_shape = c(784)) %>%
  layer_dropout(rate = 0.4) %>%
  layer_dense(units = 128, activation = 'relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 10, activation = 'softmax')

# Compile the model
model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

# Train the model
history <- model %>% fit(
  x_train, y_train,
  epochs = 30,
  batch_size = 128,
  validation_split = 0.2
)

# Evaluate the model
model %>% evaluate(x_test, y_test)

#############################
#############################

# MODEL TEST TWO


# which -a python python3
#library(reticulate)
#use_python("/usr/bin/python3")

#library(keras)
#is_keras_available()

library(tidyverse)
library(caret)

# Source all additional functions from Kevin's github repository ----------------------------------------------

source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/import_data.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/import_metadata.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/preprocessing_tool.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/export_data.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/calculate_pco.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/render_calculated_pcoa.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/sigtest.R")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/remove_last_n_columns.R")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/heatmap_dendrogram.r")

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
scale_matrix <- function(mat) {
  # Find the minimum and maximum values of the matrix
  min_val <- min(mat)
  max_val <- max(mat)
  # Scale the matrix values between 0 and 1
  scaled_mat <- (mat - min_val) / (max_val - min_val)
  return(scaled_mat)
}
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






##################################################
##################################################
# Bench
install.packages("microbenchmark") 
install.packages("ggplot2") 
install.packages("tictoc")
library(microbenchmark) 
library(ggplot2)
library(tictoc)


tic()
history <- model %>% fit(
  x_train, y_train,
  epochs = 30,
  batch_size = 128,
  validation_split = 0.2
)
toc()


results <- microbenchmark(
 
  history <- model %>% fit(
    x_train, y_train,
    epochs = 30,
    batch_size = 128,
    validation_split = 0.2
  ),
  times=2
  
)

# Print the results 
print(results) 

# Visualize the results 
autoplot(results)

##################################################
##################################################


##################################################
##################################################
# OTHER

# More shell stuff (specifically, powershell):
# A string with spaces
cmd <- "echo Hello World"

# Without shQuote
system(cmd)  # This might cause an error or unexpected behavior

# With shQuote
system(paste("powershell -Command", shQuote(cmd)))  # This ensures the string is properly quoted






# replace caret::createDataPartition	
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