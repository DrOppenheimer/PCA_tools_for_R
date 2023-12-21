sigtest <- function(data_matrix, groups.list, stat_test){
  
  "ANOVA-one-way", # c("Kruskal-Wallis", "t-test-paired", "Wilcoxon-paired", "t-test-unpaired", "Mann-Whitney-unpaired-Wilcoxon", "ANOVA-one-way"

  switch(stat_test, 
         "Kruskal-Wallis" = kruskal.test(data=data_matrix, g=groups.list,),
         "t-test-paired" = expr2,
         "Wilcoxon-paired" = ,
         "t-test-unpaired" = ,
         "Mann-Whitney-unpaired-Wilcoxon" - , 
         "ANOVA-one-way"...
         
         default = default_expr)
  
  , , , , , 
  
    
}


groups.list <- (metadata_matrix[,metadata_column])
my_stats.summary <- cbind(data_matrix, my_stats$mean, my_stats$sd, my_stats.statistic, my_stats.p, my_stats.fdr)


# retrieve the selected grouping
#groups.list <- as.list(metadata_matrix[,metadata_column])
groups.list <- (metadata_matrix[,metadata_column])
names(groups.list) <- rownames(metadata_matrix) 

mean
sd
stat
p
fdr
adjusted p

# set the working directory
setwd("/Users/kevinkeegan/Documents/GitHub/PCA_tools_for_R")

# load some data and metadat to play with
my_data <- import_data("filtered_counts.txt")
my_metadata <- import_data("filtered_counts.metadata.txt")

# Here, make sure that the data are sorted COLUMNWISE by id
my_data <-  my_data[,order(colnames(my_data))]

# make sure that the color matrix is sorted (ROWWISE) by id
my_metadata <-  my_metadata[order(rownames(my_metadata)),]

# check to make sure that the two listings of ids are identical
if( identical(colnames(my_data), rownames(my_metadata)) ){
  print("EQUAL")
}else{
  print("NOT EQUAL")
}

# run an anova for each row in my_data, i.e. for each function in the orginal mgrast based implementation

# prep someplace to store the stats for all of the rows in my_data
my_stats <- matrix(nrow = nrow(my_data), ncol=7)
rownames(my_stats) <- rownames(my_data)
colnames(my_stats) <- c("median", "mean", "sd", "F_stat", "p", "bonferroni_p", "BH_p")

# iterate through each row 
for (i in 1:nrow(my_data)){
  
  # first calculate some simple state
  my_stats[i,"median"] <- median(my_data[i,])
  my_stats[i,"mean"] <- mean(my_data[i,])
  my_stats[i,"sd"] <- sd(my_data[i,])
  
  # prep data for anova
  stat_input <- matrix(nrow=ncol(my_data), ncol=2) # a 
  colnames(stat_input) <- c("values","groups")
  stat_input[,"values"] <- my_data[i,]
  stat_input[,"groups"] <- my_metadata[,"env_package.data.body_site"]
  stat_input <- as.data.frame(stat_input)
  
  # Perform ANOVA using the ~ operator
  aov_result <- aov(values ~ groups, data = stat_input)
  
  # Summary of ANOVA results
  aov_result_summary <- summary(aov_result)
  
  # Get ANOVA results into my_stats (note FDR and adjusted p have to be calculated later)
  my_stats[i,"F_stat"] <- aov_result_summary[[1]]$`F value`[1]
  my_stats[i,"p"]    <- aov_result_summary[[1]]$`Pr(>F)`[1]
  
}
  
# Calculate the Bonferroni adjusted p
my_stats[,"bonferroni_p"] <- p.adjust(p=my_stats[,"p"], method = "bonferroni")


# Calculate the Benjamini & Hochberg adjusted p
my_stats[,"BH_p"] <- p.adjust(p=my_stats[,"p"], method = "BH")



# combine my_data and my_stats to create a single output object
my_output_data <- cbind(my_data,my_stats)

# output the object
export_data(data_object = my_output_data, file_name = "my_stat_output.txt")

anova_stats <- 
  
    # Example data
    set.seed(123)
  groupA <- rnorm(20, mean = 5)
  groupB <- rnorm(20, mean = 7)
  groupC <- rnorm(20, mean = 6)
  
  # Combine data into a data frame
  my_data <- data.frame(
    Value = c(groupA, groupB, groupC),
    Group = rep(c("A", "B", "C"), each = 20)
  )
  
  # Perform ANOVA using the ~ operator
  aov_result <- aov(Value ~ Group, data = my_data)
  
 
    
    
    
    
    
}




switch(EXPR, 
       value1 = expr1,
       value2 = expr2,
       ...
       valueN = exprN,
       default = default_expr)




  
  #my_stats <- sigtest(data_matrix, groups.list, stat_test)