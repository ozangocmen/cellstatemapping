cfa <- function(GEM, Class, GS) {
  
  ##Arguments##
  
  #  GEM: input data which is gene expression matrix (csv, tsv, and xls file formats)
  #  class: class is represent responder and non-responder samples
  #  GS: GeneSets are provided from curated signatures(eg. Pro-inflammatory, Cytolytics effector pathway)which must be a data frame.  
  #  selectdata: Select only the specified genes from the data matrix
  #  averageScore: Calculate the average gene expression level for each row in the curated signature data
  #  FactorResults: Perform factor analysis on the entire gene expression data matrix
  #  correlation: Calculate the correlation between the average expression scores of GS and the factor analysis results
  CFA <- setClass("CFA", slots = list(GEM = "data.frame",
                                      Class = "character",
                                      GS = "data.frame"))
  
  cfa <- function(CFA) {
    
    selectedData <- CFA@GEM[CFA@GS, ]
    
    averageScore <- colMeans(selectedData)
    
    factorResults <- factanal(CFA@GEM, factors = 3, rotation = "varimax")
    
    correlation <- cor(averagescore, factorResults$loadings)
    
    return(list(f1cor = correlation[1],
                f2cor = correlation[2],
                f3cor = correlation[3]))
  }
  
  
  
  #################################################
  ############ call the function ##################
  #################################################
  ## Load the gene expression data matrix
  #data <- read.csv("GEM.csv")
  
  ## Load the selected gene set
  #selectedGenes <- read.csv("GS.csv")
  
  ## Calculate the correlation between the average gene expression levels and the factor analysis results
  #correlation <- cor(GEM, GS)
  
  ## Print the correlation
  #print(correlation)
  