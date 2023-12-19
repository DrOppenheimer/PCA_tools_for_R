# Load the plotly library
library(plotly)

# Generate sample data
x <- rnorm(100)
y <- rnorm(100)
z <- rnorm(100)


setwd()
# Import PCoA
## ######################
## # SUB(2): FUNCTION TO LOAD A PRECALCULATED *.PCoA
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

# import the calculated PCoA
my_test <- load_pcoa_data("HMP.Jumpstart.DESeq_normed.euclidean.PCoA")


# Import Metadata 
metadata_matrix <- as.matrix( # Load the metadata table (same if you use one or all columns)
  read.table(
    file=metadata_table,row.names=1,header=TRUE,sep="\t",
    colClasses = "character", check.names=FALSE,
    comment.char = "",quote="",fill=TRUE,blank.lines.skip=FALSE
  )
)  

metadata_table <- "HMP_jumpstart_metadata.txt"

# covert single column to colors
metadata_colors <- create_colors(
  metadata_column = metadata_matrix[,"env_package.data.body_site", 
                                    drop=FALSE]
  )

# Create an interactive 3D scatter plot
plot_ly(
  x = my_test[["eigen_vectors"]][,1], 
  y = my_test[["eigen_vectors"]][,2], 
  z = my_test[["eigen_vectors"]][,3], 
  type = "scatter3d", 
  mode = "markers",
  marker = list(color = metadata_colors[,1])) %>%
  #marker = list(color = "red")) %>%
  layout(
    title = "Interactive 3D Scatter Plot", 
    scene = list(xaxis = list(title = "X Axis"),
                 yaxis = list(title = "Y Axis"),
                 zaxis = list(title = "Z Axis"))
)




# Load the plotly library
library(plotly)

# Generate sample data
x <- rnorm(100)
y <- rnorm(100)
z <- rnorm(100)

# Create an interactive 3D scatter plot
plot_ly(
  x = x, 
  y = y, 
  z = z, 
  type = "scatter3d", 
  mode = "markers", 
  marker = list(color = "blue")) %>%
  layout(title = "Interactive 3D Scatter Plot", scene = list(xaxis = list(title = "X Axis"),
                                                             yaxis = list(title = "Y Axis"),
                                                             zaxis = list(title = "Z Axis")))
