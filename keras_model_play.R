

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

# From https://tensorflow.rstudio.com/guides/keras/basics

model <- keras_model_sequential()

model %>%
  
  # Adds a densely-connected layer with 64 units to the model:
  layer_dense(units = 64, activation = 'relu') %>%
  
  # Add another:
  layer_dense(units = 64, activation = 'relu') %>%
  
  # Add a softmax layer with 10 output units:
  layer_dense(units = 10, activation = 'softmax')

  # Create a sigmoid layer:
  layer_dense(units = 64, activation ='sigmoid')

  # A linear layer with L1 regularization of factor 0.01 applied to the kernel matrix:
  layer_dense(units = 64, kernel_regularizer = regularizer_l1(0.01))

  # A linear layer with L2 regularization of factor 0.01 applied to the bias vector:
  layer_dense(units = 64, bias_regularizer = regularizer_l2(0.01))

  # A linear layer with a kernel initialized to a random orthogonal matrix:
  layer_dense(units = 64, kernel_initializer = 'orthogonal')

  # A linear layer with a bias vector initialized to 2.0:
  layer_dense(units = 64, bias_initializer = initializer_constant(2.0))

  # After the model is constructed, configure its learning process by calling the compile method:
    model %>% compile(
      optimizer = 'adam',
      loss = 'categorical_crossentropy',
      metrics = list('accuracy')
    )
  # compile takes three important arguments:
    # optimizer: This object specifies the training procedure. Commonly used optimizers are e.g.
    # adam, rmsprop, or sgd.
    # loss: The function to minimize during optimization. Common choices include mean square error (mse), categorical_crossentropy, and binary_crossentropy.
    # metrics: Used to monitor training. In classification, this usually is accuracy.
     
    # Configure a model for mean-squared error regression.
    model %>% compile(
      optimizer = 'adam',
      loss = 'mse',           # mean squared error
      metrics = list('mae')   # mean absolute error
    )
    
    # Configure a model for categorical classification.
    model %>% compile(
      optimizer = optimizer_rmsprop(learning_rate = 0.01),
      loss = "categorical_crossentropy",
      metrics = list("categorical_accuracy")
    )
  
    # Input data
    # 
    # You can train keras models directly on R matrices and arrays (possibly created from R data.frames). 
    # A model is fit to the training data using the fit method:
    #   
    data <- matrix(rnorm(1000 * 32), nrow = 1000, ncol = 32)
    labels <- matrix(rnorm(1000 * 10), nrow = 1000, ncol = 10)
 
    model %>% fit(
      data,
      labels,
      epochs = 10,
      batch_size = 32
    )
    # fit takes three important arguments:
    #   
    # epochs: Training is structured into epochs. An epoch is one iteration over the entire input 
    #  data (this is done in smaller batches).
    # batch_size: When passed matrix or array data, the model slices the data into smaller batches 
    #  and iterates over these batches during training. This integer specifies the size of each batch. 
    #  Be aware that the last batch may be smaller if the total number of samples is not divisible by the batch size.
    # validation_data: When prototyping a model, you want to easily monitor its performance on some validation data. 
    #  Passing this argument — a list of inputs and labels — allows the model to display the loss and metrics in inference mode for the passed data, at the end of each epoch.
    # 

    #Here’s an example using validation_data:
      
      data <- matrix(rnorm(1000 * 32), nrow = 1000, ncol = 32)
    labels <- matrix(rnorm(1000 * 10), nrow = 1000, ncol = 10)
    
    val_data <- matrix(rnorm(1000 * 32), nrow = 100, ncol = 32)
    
    val_labels <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
    
    model %>% fit(
      data,
      labels,
      epochs = 10,
      batch_size = 32,
      validation_data = list(val_data, val_labels)
    )
    
    # Same as fit, the evaluate and predict methods can use raw R data as well as a dataset.
    # 
    # To evaluate the inference-mode loss and metrics for the data provided:
      
      model %>% evaluate(test_data, test_labels, batch_size = 32)
    
    model %>% evaluate(test_dataset, steps = 30)
    
    # And to predict the output of the last layer in inference for the data provided, again as R data as well as a dataset:
    #   
      model %>% predict(test_data, batch_size = 32)
    
    model %>% predict(test_dataset, steps = 30)
    
    
    
    # Functional API
    # 
    # The sequential model is a simple stack of layers that cannot represent arbitrary models. Use the Keras functional API to build complex model topologies such as:
    #   
    #   multi-input models,
    # multi-output models,
    # models with shared layers (the same layer called several times),
    # models with non-sequential data flows (e.g., residual connections).
    # 
    # Building a model with the functional API works like this:
    #   
    #   A layer instance is callable and returns a tensor.
    # Input tensors and output tensors are used to define a keras_model instance.
    # This model is trained just like the sequential model.
    # 
    # The following example uses the functional API to build a simple, fully-connected network:
    #   
      inputs <- layer_input(shape = (32))  # Returns a placeholder tensor
    
    predictions <- inputs %>%
      layer_dense(units = 64, activation = 'relu') %>%
      layer_dense(units = 64, activation = 'relu') %>%
      layer_dense(units = 10, activation = 'softmax')
    
    # Instantiate the model given inputs and outputs.
    model <- keras_model(inputs = inputs, outputs = predictions)
    
    # The compile step specifies the training configuration.
    model %>% compile(
      optimizer = optimizer_rmsprop(lr = 0.001),
      loss = 'categorical_crossentropy',
      metrics = list('accuracy')
    )
    
    # Trains for 5 epochs
    model %>% fit(
      data,
      labels,
      batch_size = 32,
      epochs = 5
    )
    

# SECOND EXAMPLE ----------------------------------------------------------
# from https://cran.r-project.org/web/packages/keras/vignettes/
    # We can learn the basics of Keras by walking through a simple example: recognizing 
    # handwritten digits from the MNIST dataset. MNIST consists of 28 x 28 grayscale images 
    # of handwritten digits like these:
    #   
    #   
    #   
    #   The dataset also includes labels for each image, telling us which digit it is. For example, 
    # the labels for the above images are 5, 0, 4, and 1.
    # Preparing the Data
    # The MNIST dataset is included with Keras and can be accessed using the dataset_mnist() 
    # function. Here we load the dataset then create variables for our test and training data:
      
      library(keras)
    mnist <- dataset_mnist()
    x_train <- mnist$train$x
    y_train <- mnist$train$y
    x_test <- mnist$test$x
    y_test <- mnist$test$y   
    
    # The x data is a 3-d array (images,width,height) of grayscale values . 
    # To prepare the data for training we convert the 3-d arrays into matrices by reshaping 
    # width and height into a single dimension (28x28 images are flattened into length 784 vectors). 
    # Then, we convert the grayscale values from integers ranging between 0 to 255 into floating point 
    # values ranging between 0 and 1:
    
    # reshape
    x_train <- array_reshape(x_train, c(nrow(x_train), 784))
    x_test <- array_reshape(x_test, c(nrow(x_test), 784))
    # rescale
    x_train <- x_train / 255
    x_test <- x_test / 255
    
    # Note that we use the array_reshape() function rather than the dim<-() function to reshape the array. 
    # This is so that the data is re-interpreted using row-major semantics (as opposed to R’s default 
    # column-major semantics), which is in turn compatible with the way that the numerical libraries 
    # called by Keras interpret array dimensions.
    # 
    # The y data is an integer vector with values ranging from 0 to 9. To prepare this data for training 
    # we one-hot encode the vectors into binary class matrices using the Keras to_categorical() function:
    
    y_train <- to_categorical(y_train, 10)
    y_test <- to_categorical(y_test, 10)
    
    # Defining the Model
    # The core data structure of Keras is a model, a way to organize layers. The simplest type of model 
    # is the Sequential model, a linear stack of layers.
    # 
    # We begin by creating a sequential model and then adding layers using the pipe (%>%) operator:
    #   
      model <- keras_model_sequential() 
    model %>% 
      layer_dense(units = 256, activation = 'relu', input_shape = c(784)) %>% 
      layer_dropout(rate = 0.4) %>% 
      layer_dense(units = 128, activation = 'relu') %>%
      layer_dropout(rate = 0.3) %>%
      layer_dense(units = 10, activation = 'softmax')
    # The input_shape argument to the first layer specifies the shape of the input data 
    # (a length 784 numeric vector representing a grayscale image). The final layer outputs a 
    # length 10 numeric vector (probabilities for each digit) using a softmax activation function.
    # 
    # Use the summary() function to print the details of the model:
      
      summary(model)
    
      # Next, compile the model with appropriate loss function, optimizer, and metrics:
        
        model %>% compile(
          loss = 'categorical_crossentropy',
          optimizer = optimizer_rmsprop(),
          metrics = c('accuracy')
        )
        # Training and Evaluation
        # Use the fit() function to train the model for 30 epochs using batches of 128 images:
          
          history <- model %>% fit(
            x_train, y_train, 
            epochs = 30, batch_size = 128, 
            validation_split = 0.2
          )
          # The history object returned by fit() includes loss and accuracy metrics which we can plot:
            
            plot(history)
    
            # Evaluate the model’s performance on the test data:
              
              model %>% evaluate(x_test, y_test)
            # $loss
            # [1] 0.1149
            # 
            # $acc
            # [1] 0.9807
            # Generate predictions on new data:
              
              model %>% predict_classes(x_test)
            # [1] 7 2 1 0 4 1 4 9 5 9 0 6 9 0 1 5 9 7 3 4 9 6 6 5 4 0 7 4 0 1 3 1 3 4 7 2 7 1 2
            # [40] 1 1 7 4 2 3 5 1 2 4 4 6 3 5 5 6 0 4 1 9 5 7 8 9 3 7 4 6 4 3 0 7 0 2 9 1 7 3 2
            # [79] 9 7 7 6 2 7 8 4 7 3 6 1 3 6 9 3 1 4 1 7 6 9
            # [ reached getOption("max.print") -- omitted 9900 entries ]
            # Keras provides a vocabulary for building deep learning models that is simple, elegant, and intuitive. Building a question answering system, an image classification model, a neural Turing machine, or any other model is just as straightforward.

# restart R session programmatically --------------------------------------
.rs.restartR() 



mnist <- dataset_mnist()
x_train <- mnist$train$x
y_train <- mnist$train$y
x_test <- mnist$test$x
y_test <- mnist$test$y




# CRAP BELOW HERE ---------------------------------------------------------



reticulate::install_miniconda()



library(reticulate)
library(keras3)
library(tensorflow)