# Load the plotly library
library(plotly)
library(scatterplot3d)

setwd("~/Documents/GitHub/PCA_tools_for_R/")




# WORKFLOW ----------------------------------------------------------------


# import the calculated PCoA
my_test_pcoa <- load_pcoa_data("HMP.Jumpstart.DESeq_normed.euclidean.PCoA")

# import metadata
my_metadata <- load_metadata("HMP_jumpstart_metadata.txt")

# select a single column of metadata for color generation
column_name = "env_package.data.body_site"
#column_name = "mixs.seq_method"
#column_name = "library.data.seq_center"
metadata_column = metadata_matrix[,column_name, drop=FALSE]
# generate colors for selected column
metadata_colors <- create_colors(metadata_column)

# generate the interactive 3d plot
plot_interactive_colored_3d_pcoa(
    pcoa_data = my_test_pcoa,
    selected_eigen_vectors = c(1,2,3),
    my_metadata_colors = metadata_colors
)

# iterate through the metadata creating a static 3d plot for each metadata column
plot_static_colored_3d_pcoas(
    pcoa_filename = "HMP.Jumpstart.DESeq_normed.euclidean.PCoA",
    selected_eigen_vectors = c(1,2,3),
    metadata_filename = "HMP_jumpstart_metadata.txt",
    debug = TRUE
)




####################################

# FUNCTIONS ---------------------------------------------------------------

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
plot_interactive_colored_3d_pcoa <- function(
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



## ######################
## # SUB(5): FUNCTION TO GENERATE STATIC 3D PLOTS FROM ALL OF THE METADATA 
## ######################
######################
# SUB(1): Workhorse function that creates the 3d plot
######################
plot_static_colored_3d_pcoas <- function(
  pcoa_filename = my_test_pcoa,
  selected_eigen_vectors = c(1,2,3),
  metadata_filename = my_test_metadata,
  debug = TRUE
){
  
  # import the calculated PCoA
  pcoa_data <- load_pcoa_data(pcoa_filename)
  
  # import metadata
  my_metadata <- load_metadata(metadata_filename)
  
  # iterate through the metadata
  # create plot and color it automatically
  for (i in 1:ncol(my_metadata)){
    
    metadata_column_name <- colnames(my_metadata)[i]
    metadata_column = my_metadata[,metadata_column_name, drop=FALSE]
    my_output_png_filename <- paste(pcoa_filename, ".", metadata_column_name, ".png", sep = "")
    print(metadata_column_name)
    column_color = create_colors(metadata_column)
    if(debug==TRUE){ print(paste(dim(column_color))); testies <<- column_color }
    
    png(
      filename=my_output_png_filename, 
      width = 8, 
      height = 8, 
      units = "in", 
      res = 300) # 300 dpi
      scatterplot3d(
        x = pcoa_data[["eigen_vectors"]][,selected_eigen_vectors[1]], 
        y = pcoa_data[["eigen_vectors"]][,selected_eigen_vectors[2]], 
        z = pcoa_data[["eigen_vectors"]][,selected_eigen_vectors[3]], 
        color = column_color[,1], 
        pch = 19, 
        main = metadata_column_name,
        xlab = paste("xPC_",selected_eigen_vectors[1]), 
        ylab = paste("yPC_",selected_eigen_vectors[2]), 
        zlab = paste("zPC_",selected_eigen_vectors[3]))
    dev.off()
  
  }
  
}
  
   
  
    
  
}

#pdf("plot_high_res.pdf", width = 8, height = 8)
#scatterplot3d(x, y, z, color = "blue", pch = 19, xlab = "X Axis", ylab = "Y Axis", zlab = "Z Axis")
#dev.off()




# select a single column of metadata for color generation
column_name = "env_package.data.body_site"
#column_name = "mixs.seq_method"
#column_name = "library.data.seq_center"
metadata_column = metadata_matrix[,column_name, drop=FALSE]
# generate colors for selected column
metadata_colors <- create_colors(metadata_column)












create_plot <- function(
    PCoA_in,
    ncol.color_matrix,
    eigen_values, eigen_vectors, 
    components, 
    plot_samples,
    column_levels, num_levels, color_levels, 
    pcoa_colors, 
    plot_pch, pch_labels, pch_levels,
    image_out,
    figure_main,
    image_width_in, image_height_in, image_res_dpi,
    width_legend, 
    width_figure,
    title_cex, legend_cex, figure_cex, figure_symbol_cex, bar_cex, 
    bar_vert_adjust, label_points, vert_line, 
    debug
){
  
  if(debug==TRUE){print("creating figure")}
  
  png( # initialize the png 
    filename = image_out,
    width = image_width_in,
    height = image_height_in,
    res = image_res_dpi,
    units = 'in'
  )
  
  # LAYOUT CREATION HAS TO BE DICTATED BY PCH TO A DEGREE _ NUM LEVELS (1 or more)
  # Determine num levels for pch
  num_pch <- length(levels(as.factor(plot_pch)))
  # CREATE THE LAYOUT
  if ( num_pch > 1 ){
    my_layout <- layout( matrix(c(1,1,2,3,4,3,5,5), 4, 2, byrow=TRUE ), widths=c(0.5,0.5), heights=c(0.1,0.8,0.3,0.1) )
  }else{
    my_layout <- layout(  matrix(c(1,1,2,3,4,4), 3, 2, byrow=TRUE ), widths=c(width_legend,width_figure), heights=c(0.1,0.8,0.1) )
    # requires an extra plot.new() to skip over pch legend (frame 4 or none )
  }
  # my_layout <- layout(  matrix(c(1,1,2,3,4,3,5,5), 4, 2, byrow=TRUE ), widths=c(width_legend,width_figure), heights=c(0.1,0.4,0.8,0.4,0.1) ) # for auto pch legend
  layout.show(my_layout)
  
  # PLOT THE TITLE (layout frame 1)
  par( mai = c(0,0,0,0) )
  par( oma = c(0,0,0,0) )
  plot.new()
  if ( identical(title_cex, "default") ){ # automatically scale cex for the legend
    if(debug==TRUE){print("autoscaling the title cex")}
    title_par <- par_fetch()
    title_cex <- calculate_cex(figure_main, title_par$my_pin, title_par$my_mai, reduce_by=0.10)
  }
  text(x=0.5, y=0.5, figure_main, cex=title_cex)
  
  # PLOT THE LEGEND (layout frame 2)
  plot.new()
  if ( identical(legend_cex, "default") ){ # automatically scale cex for the legend
    if(debug==TRUE){print("autoscaling the legend cex")}
    legend_par <- par_fetch()
    legend_cex <- calculate_cex(column_levels, legend_par$my_pin, legend_par$my_mai, reduce_by=0.40)
  }
  legend( x="center", y="center", legend=column_levels, pch=15, col=color_levels, cex=legend_cex)
  
  # PLOT THE PCoA FIGURE (layout frame 3)
  # set par options (Most of the code in this section is copied/adapted from Dan Braithwaite's pco plotting in matR)
  
  #par(op)
  par <- list ()
  #par$mar <- par()['mar']
  #par$oma <- par()['oma']
  #par$mar <- c(4,4,4,4)
  #par$mar <- par(op)['mar']
  #par$oma <- par(op)['oma']
  #par$oma <- c(1,1,1,1)
  #par$mai <- c(1,1,1,1)
  par$main <- ""#figure_main
  #par$labels <- if (length (names (x)) != 0) names (x) else samples (x)
  if ( label_points==TRUE ){
    #par$labels <-  rownames(eigen_vectors)
    if ( identical(plot_samples, "all") ){
      par$labels <-  rownames(eigen_vectors)
    }else{
      par$labels <- plot_samples
    }
  } else {
    par$labels <- NA
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
