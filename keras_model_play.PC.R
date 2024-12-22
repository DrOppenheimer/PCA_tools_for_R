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
