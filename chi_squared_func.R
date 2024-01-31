# Create a contingency table 
data <- matrix(c(120, 90, 40, 110, 95, 45, 80, 70, 30), ncol = 3, byrow = TRUE) 
colnames(data) <- c("Group A", "Group B", "Group C") 
rownames(data) <- c("Male", "Female", "Other") 
# Perform the chi-square test 
result <- chisq.test(data) 
# Print the results 
cat("Chi-Square Test Statistic: ", result$statistic, "\n") 
cat("Degrees of Freedom: ", result$parameter, "\n") 
cat("p-value: ", result$p.value, "\n")





import_metadata <- function(group_table){ #, group_column, sample_names){
  metadata_matrix <- as.matrix( # Load the metadata table (same if you use one or all columns)
    read.table(
      file=group_table,row.names=1,header=TRUE,sep="\t",
      colClasses = "character", check.names=FALSE,
      comment.char = "",quote="",fill=TRUE,blank.lines.skip=FALSE
    )
  )
}

import_data <- function(file_name)
{
  data.matrix(read.table(file_name, row.names=1, header=TRUE, sep="\t", comment.char="", quote="", check.names=FALSE))
} 


data_file="filtered_counts.txt" 
metadata_file="filtered_counts.metadata.txt" 
metadata_column="env_package.data.env_package"
stat_test ="Chi-Square"
p_adjust_method = "BH"

results <- perform_chi(
    data_file, #="filtered_counts.txt", 
    metadata_file, #="filtered_counts.metadata.txt", 
    metadata_column, #="env_package.data.body_site", 
    stat_test, #="Chi-Square",
    p_adjust_method # = "BH"
)






# load some data and metadata to play with
my_data <- import_data(data_file)
my_metadata <- import_metadata(metadata_file)
my_metadata <<- my_metadata

# Here, make sure that the data are sorted COLUMNWISE by id
my_data <-  my_data[,order(colnames(my_data))]

# make sure that the metadata matrix is sorted (ROWWISE) by id
my_metadata <-  my_metadata[order(rownames(my_metadata)),]

# check to make sure that the two listings of ids are identical
if( identical(colnames(my_data), rownames(my_metadata)) ){
  print("IDs are EQUAL")
}else{
  stop("IDs are NOT EQUAL, check your data and metadata files")
}






perform_chi <- function(
    data_file, #="filtered_counts.txt", 
    metadata_file, #="filtered_counts.metadata.txt", 
    metadata_column, #="env_package.data.body_site", 
    stat_test, #="Chi-Square",
    p_adjust_method # = "BH"
){
  
  # prep someplace to store the stats for all of the rows in my_data
  my_stats <- matrix(nrow = nrow(my_data), ncol=7)
  rownames(my_stats) <- rownames(my_data)
  
  # parts that need to be part of the switch
  colnames(my_stats) <- c("median", "mean", "sd", "X-squared_stat", "p", "bonferroni_p", paste(p_adjust_method,"_p",sep=""))
  
  # iterate through each row 
  for (i in 1:nrow(my_data)){
    
    # first calculate some simple state
    my_stats[i,"median"] <- median(my_data[i,])
    my_stats[i,"mean"] <- mean(my_data[i,])
    my_stats[i,"sd"] <- sd(my_data[i,])
    
    # prep data for chi_squared
    stat_input <- matrix(nrow=ncol(my_data), ncol=2) # a 
    colnames(stat_input) <- c("values","groups")
    stat_input[,"values"] <- my_data[i,]
    stat_input[,"groups"] <- my_metadata[,metadata_column]
    stat_input <- as.data.frame(stat_input)
    
    # Perform chi_squared using the ~ operator
    #aov_result <- aov(values ~ groups, data = stat_input)
    chi_result <- chisq.test(stat_input)
    
    # Summary of ANOVA results
    #aov_result_summary <- summary(chi_result)
    
    # Get Chi results into my_stats (note FDR and adjusted p have to be calculated later)
    my_stats[i,"X-squared_stat"] <- chi_result$statistic[[1]]
    my_stats[i,"p"]      <- chi_result$p.value[[1]]
    
  }
  
  # Calculate the Bonferroni adjusted p
  my_stats[,"bonferroni_p"] <- p.adjust(p=my_stats[,"p"], method = "bonferroni")
  
  # Calculate the Benjamini & Hochberg adjusted p
  my_stats[,paste(p_adjust_method,"_p",sep="")] <- p.adjust(p=my_stats[,"p"], method = p_adjust_method)
  
  # combine my_data and my_stats to create a single output object
  my_output_data <- cbind(my_data,my_stats)
  
  # sort the data by p value  
  my_output_data <- my_output_data[order(my_output_data[, "p"]), ]
  
  # output the object
  export_data(
    data_object = my_output_data, 
    file_name = paste(tools::file_path_sans_ext(data_file),".",stat_test,".", metadata_column, ".STAT_RESULTS.txt", sep="")
  )
  
}