sigtest <- function(data_file="filtered_counts.txt", 
                    metadata_file="filtered_counts.metadata.txt", 
                    metadata_column="env_package.data.body_site", 
                    stat_test="ANOVA-one-way", # c("Kruskal-Wallis","t-test-paired","Wilcoxon-paired","ANOVA-one-way","t-test-unpaired" = ,"Mann-Whitney-unpaired-Wilcoxon") 
                    p_adjust_method = "BH" # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")
                    ){
  
  # load some data and metadata to play with
  my_data <- import_data(data_file)
  my_metadata <- import_data(metadata_file)
  
  # Here, make sure that the data are sorted COLUMNWISE by id
  my_data <-  my_data[,order(colnames(my_data))]
  
  # make sure that the color matrix is sorted (ROWWISE) by id
  my_metadata <-  my_metadata[order(rownames(my_metadata)),]
  
  # check to make sure that the two listings of ids are identical
  if( identical(colnames(my_data), rownames(my_metadata)) ){
    print("EQUAL")
  }else{
    stop("NOT EQUAL")
  }

  switch(stat_test, 
         "Kruskal-Wallis" = print("test"),
         "t-test-paired" = print("test"),
         "Wilcoxon-paired" = print("test"),
         "t-test-unpaired" = print("test"),
         "Mann-Whitney-unpaired-Wilcoxon" = print("test"), 
         "ANOVA-one-way" = perform_anova(data_file,metadata_file,metadata_column,stat_test,p_adjust_method)
         )
}
         
  
    



## set the working directory
#setwd("/Users/kevinkeegan/Documents/GitHub/PCA_tools_for_R")



# run stat test for each row in my_data, i.e. for each function in the original mgrast based implementation


switch(EXPR, 
       value1 = expr1,
       value2 = expr2,
       ...
       valueN = exprN,
       default = default_expr)




perform_anova <- function(
    data_file="filtered_counts.txt", 
    metadata_file="filtered_counts.metadata.txt", 
    metadata_column="env_package.data.body_site", 
    stat_test="anova",
    p_adjust_method = "BH"
  ){

  # prep someplace to store the stats for all of the rows in my_data
  my_stats <- matrix(nrow = nrow(my_data), ncol=7)
  rownames(my_stats) <- rownames(my_data)
  
  # parts that need to be part of the switch
  colnames(my_stats) <- c("median", "mean", "sd", "F_stat", "p", "bonferroni_p", paste(p_adjust_method,"_p",sep=""))

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
    my_stats[i,"p"]      <- aov_result_summary[[1]]$`Pr(>F)`[1]
  
  }
  
  # Calculate the Bonferroni adjusted p
  my_stats[,"bonferroni_p"] <- p.adjust(p=my_stats[,"p"], method = "bonferroni")

  # Calculate the Benjamini & Hochberg adjusted p
  my_stats[,"BH_p"] <- p.adjust(p=my_stats[,"p"], method = p_adjust_method)

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






  
  #my_stats <- sigtest(data_matrix, groups.list, stat_test)