# Source all additional functions from Kevin's github repository ----------------------------------------------

source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/import_data.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/import_metadata.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/preprocessing_tool.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/export_data.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/calculate_pco.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/render_calculated_pcoa.r")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/sigtest.R")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/remove_last_n_columns.R")
source("https://raw.githubusercontent.com/DrOppenheimer/workflow_play/master/heatmap_dendrogram.r")# Install required packages

#########################################################################################################

# Install packages and set the environment
install.packages("keras")
remotes::install_github("rstudio/keras")
install.packages("kerastuneR")
install.packages("tfruns")


# Load the packages
library(keras)
#library(kerastuneR)
library(reticulate)
library(tfruns)
# Load the correct python environment (see keras_model_play.PC.R)
use_python("C:/Users/kosso/AppData/Local/r-miniconda/envs/PythonCPU")


# Prepare the Data (Using MNIST Dataset)
mnist <- dataset_mnist()
x_train <- mnist$train$x / 255
y_train <- mnist$train$y
x_val <- mnist$test$x / 255
y_val <- mnist$test$y

x_train <- array_reshape(x_train, c(nrow(x_train), 28, 28, 1))
x_val <- array_reshape(x_val, c(nrow(x_val), 28, 28, 1))

#####################################################################
# Define the Model with Hyperparameters Create a script file build_model.R with the following content:
FLAGS <- flags(
  flag_integer("units", 128),
  flag_numeric("learning_rate", 0.001)
)

build_model <- function() {
  model <- keras_model_sequential() %>%
    layer_flatten(input_shape = c(28, 28, 1)) %>%
    layer_dense(units = FLAGS$units, activation = 'relu') %>%
    layer_dense(units = 10, activation = 'softmax')
  
  model %>% compile(
    optimizer = optimizer_adam(lr = FLAGS$learning_rate),
    loss = 'sparse_categorical_crossentropy',
    metrics = c('accuracy')
  )
  
  return(model)
}

model <- build_model()
history <- model %>% fit(
  x_train, y_train,
  epochs = 10,
  validation_data = list(x_val, y_val)
)
#######################################################################

# Set Up Hyperparameter Tuning Define the hyperparameter search space and perform the search:
runs <- tuning_run(
  "build_model.R",
  flags = list(
    units = c(32, 64, 128, 256, 512),
    learning_rate = c(1e-4, 1e-3, 1e-2)
  )
)


# Running the Hyperparameter Tuning
# To execute the tuning, run the following command:
runs

# Convert to data frame and save as CSV 
#runs_df <- as.data.frame(runs) 
#write.csv(runs_df, file = "tuning_results.csv", row.names = FALSE)
# Runs already is a dataframe
write.csv(runs, file = "tuning_results_2.csv", row.names = FALSE)



######################################################################
# Example 2 with PCoA data
######################################################################

### First part, wrangle the data
library(tidyverse)
library(caret)
library(reticulate)
library(tfruns)
# Load the correct python environment (see keras_model_play.PC.R)
use_python("C:/Users/kosso/AppData/Local/r-miniconda/envs/PythonCPU")

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

## Prepare the data
#train_labels <- to_categorical(as.numeric(training_data$env_package.data.body_site) - 1)
#train_data <- as.matrix(training_data %>% select(-env_package.data.body_site))
#rownames(train_data) <- 1:nrow(train_data) # Keras doesn't like rownames

#test_labels <-  to_categorical(as.numeric(testing_data$env_package.data.body_site) - 1)
#test_data <- as.matrix(testing_data %>% select(-env_package.data.body_site))
#rownames(test_data) <- 1:nrow(test_data)

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
num_input_neurons
# How big to make the hidden layer# # from # https://medium.com/geekculture/introduction-to-neural-network-2f8b8221fbd3#:~:text=Number%20of%20Neurons%20and%20Number%20of%20Layers%20in%20Hidden%20Layer&text=The%20number%20of%20hidden%20neurons,size%20of%20the%20output%20layer.
num_hidden_neurons <- ceiling( (2/3*num_coord) + length(unique(my_metadata[,"env_package.data.body_site"])) )
num_hidden_neurons
# How big to make the output layer
num_output_neurons <- length(unique(my_metadata[,"env_package.data.body_site"]))
num_output_neurons

### Part two, build model that will go into my_model.R
###########################################################################
# Define your model 
# --> See my_model.R

###########################################################################

# Part Three: set up function to run the tuning and run it
runs <- tuning_run(
  "my_model.R",
  flags = list(
    units = c(26, 50),
    learning_rate = c(1e-4, 1e-3, 1e-2),
    activation = c("relu", "tanh", "sigmoid"),
    epochs = c(100, 2000)
  )
)

###########################################################################

# save the final output
write.csv(runs, file = "tuning_results_3.csv", row.names = FALSE)
