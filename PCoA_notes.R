# Load the plotly library
library(plotly)

# Generate sample data
x <- rnorm(100)
y <- rnorm(100)
z <- rnorm(100)


setwd("~/Documents/GitHub/PCA_tools_for_R/")



## Workflow

# WORKFLOW ----------------------------------------------------------------


# import the calculated PCoA
my_test_pcoa <- load_pcoa_data("HMP.Jumpstart.DESeq_normed.euclidean.PCoA")

# import metadata
my_metadata <- load_metadata("HMP_jumpstart_metadata.txt")

# select a single column of metadata for color generation
column_name = "env_package.data.body_site"
metadata_column = metadata_matrix[,column_name, drop=FALSE]
# generate colors for selected column
metadata_colors <- create_colors(metadata_column)

# generate the interactive 3d plot
plot_colored_3d_pcoa(
    pcoa_data = my_test_pcoa,
    selected_eigen_vectors = c(1,2,3),
    my_metadata_colors = metadata_colors
)





####################################

# FUNCTIONS ---------------------------------------------------------------



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



### Functions

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



# Import metadata
## ######################
## # SUB(2): FUNCTION TO LOAD METADATA FILE
## ######################
load_metadata <- function(metadata_file){
  metadata_matrix <- as.matrix( # Load the metadata table (same if you use one or all columns)
    read.table(
      file=metadata_file,row.names=1,header=TRUE,sep="\t",
      colClasses = "character", check.names=FALSE,
      comment.char = "",quote="",fill=TRUE,blank.lines.skip=FALSE
    )
  ) 
  return(metadata_matrix)
}




## ######################
## # SUB(3): FUNCTION TO GENERATE COLORS FOR THE PLOT 
## ######################
create_colors <- function(metadata_column, color_mode = "auto"){ # function to     
  my_data.color <- data.frame(metadata_column)
  column_factors <- as.factor(metadata_column[,1])
  column_levels <- levels(as.factor(metadata_column[,1]))
  num_levels <- length(column_levels)
  color_levels <- col.wheel(num_levels)
  levels(column_factors) <- color_levels
  my_data.color[,1]<-as.character(column_factors)
  return(my_data.color)
}


## ######################
## # SUB(4): FUNCTION TO GENERATE THE INTERACTIVE 3D PLOT 
## ######################
# Create an interactive 3D scatter plot using the selected_eigen_vectors and metadata_colors
plot_colored_3d_pcoa <- function(
    pcoa_data = my_test_pcoa,
    selected_eigen_vectors = c(1,2,3),
    my_metadata_colors = metadata_colors
){
  plot_ly(
    x = pcoa_data[["eigen_vectors"]][,selected_eigen_vectors[1]], 
    y = pcoa_data[["eigen_vectors"]][,selected_eigen_vectors[2]], 
    z = pcoa_data[["eigen_vectors"]][,selected_eigen_vectors[3]], 
    type = "scatter3d", 
    mode = "markers",
    marker = list(color = my_metadata_colors[,1])) %>%
    layout(
      title = "Interactive 3D Scatter Plot", 
      scene = list(xaxis = list(title = "X Axis"),
                   yaxis = list(title = "Y Axis"),
                   zaxis = list(title = "Z Axis"))
    )
}













############ JUST CRAP BELOW HERE
############ JUST CRAP BELOW HERE
############ JUST CRAP BELOW HERE





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
