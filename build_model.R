library(keras)

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