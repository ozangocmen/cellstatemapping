  ###this function is related to common factor analysis 
  
  cfa <- function(GEM, Class, GS) {
    
      ##Arguments##
    
      #  GEM: input data which is gene expression matrix (csv, tsv and xls file formats)
      #  class: class is represent to responder and non-responder samples
      #  GS: GeneSets are provides from curated signatures(eg. Pro inflammatory, Cytolytics effector pathway)which is must be a dataframe.  
      #  selectdata: Select only the specified genes from the data matrix
      #  averagescore: Calculate the average gene expression level for each row in the curated signature data
      #  results: Perform factor analysis on the entire gene expression data matrix
      #  correlation: Calculating the correlation between the average expression scores of GS and the factor analysis results
 
  cfa <- setClass(Class = c("responder","non-responder",  
                  slots = c(field1 = "data.frame",
                            field2 = "data.frame")
    
      selectedData <- data[GS, ]
      
      averagescore <- rowMeans(GS)
      
      results <- list(averageLevels = apply(GS, 1, mean),
                      factorResults = factanal(GEM, factors = 3, rotation = "varimax")
      
      correlation <- cor(averagescore, factorResults$loadings)
    
      return(f1cor = x, f2cor = y, f3cor = z)
    }
    

    

##### calling function #####
## Load the gene expression data matrix
  #data <- read.csv("GEM.csv")
  
## Load the selected gene set
  #selectedGenes <- read.csv("GS.csv")
  
## Calculate the correlation between the average gene expression levels and the factor analysis results
  #correlation <- cor(GEM, GS)
  
## Print the correlation
  #print(correlation)
  
