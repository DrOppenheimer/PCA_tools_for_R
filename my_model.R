# Figure out the number of neurons in each layer
# how big to make the input layer # https://stats.stackexchange.com/questions/181/how-to-choose-the-number-of-hidden-layers-and-nodes-in-a-feedforward-neural-netw
num_input_neurons <- ncol(all_data) - 1
# How big to make the hidden layer# # from # https://medium.com/geekculture/introduction-to-neural-network-2f8b8221fbd3#:~:text=Number%20of%20Neurons%20and%20Number%20of%20Layers%20in%20Hidden%20Layer&text=The%20number%20of%20hidden%20neurons,size%20of%20the%20output%20layer.
num_hidden_neurons <- ceiling( (2/3*num_coord) + length(unique(my_metadata[,"env_package.data.body_site"])) )
# How big to make the output layer
num_output_neurons <- length(unique(my_metadata[,"env_package.data.body_site"]))

FLAGS <- flags(
  flag_integer("units", 128),
  flag_numeric("learning_rate", 0.001),
  flag_string("activation", "relu"),
  flag_integer("epochs", 10)
)

build_model <- function() {
  model <- keras_model_sequential() %>%
    layer_dense(units = num_input_neurons, activation = 'relu', input_shape = ncol(all_data)) %>%
    layer_dropout(rate = 0.3) %>%
    layer_dense(units = FLAGS$units, activation = FLAGS$activation) %>% # num_hidden_neurons
    layer_dropout(rate = 0.3) %>%
    layer_dense(units = num_output_neurons, activation = 'softmax')
    
  model %>% compile(
    loss = "categorical_crossentropy",
    optimizer = optimizer_rmsprop(),
    metrics = c("accuracy")
  )
  
  return(model)
}

#print(paste("Dim all_labels", dim(all_labels)))

model <- build_model()
history <- model %>% fit(
  all_data, all_labels,
  epochs = FLAGS$epochs,
  validation_split = 0.2
)

## Calculate some custom metrics validation accuracy 
#avg_val_accuracy <- mean(history$metrics$val_accuracy) # avg value accuracy
#write_run_metadata("avg_val_accuracy", avg_val_accuracy)
#max_val_accuracy <- max(history$metrics$val_accuracy) # max value accuracy
#write_run_metadata("max_val_accuracy", max_val_accuracy)
#min_val_accuracy <- min(history$metrics$val_accuracy) # min value accuracy
#write_run_metadata("min_val_accuracy", min_val_accuracy)