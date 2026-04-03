
####################### DEFINE THE FUNCTIONS  #######################

#' This function returns a list of two arrays. The first is the Pearson correlation coefficients of 
#' mean gene expression of genes in a given gene set (Gs). The second is the df of Factor loading values per cell.
#' GEM: gene expression matrix. Has to be genes (rows) by cells (cols).
#' Gs: Gene set. Gene symbols have to be present the GEM row names.
#' nF: Number of factors to be run
RunCFA_GScors <- function( GEM, Gs, nF=3){
  
  #1st part: Compute the common factors
  fct <- factanal(GEM, nF, rotation="varimax")
  factors <- as.data.frame(head(fct$loadings, nrow(GEM)))
  #2nd part: Compute the mean Gs expression per cell
  GS.Cell.mean <- apply(GEM[which(rownames(GEM) %in% Gs), ], 2, mean, na.rm=T)#this requires tidyverse. update.
  
  GS.F.cors <- c()
  for(i in 1:ncol(factors)){
    GS.F.cors[i] <- cor(GS.Cell.mean, factors[[i]])
  }
  return(list(FactCors=GS.F.cors, Factors=factors))
}

#' GEM: gene expression matrix. Has to be genes (rows) by cells (cols).
#' Obs.cor: Observed Pearson correlation coefficient between one Factor and mean gene expression of a GS
#' Gs: Gene set. Needed for only size
#' Factor: An array of factor loadings per cell corresponing to the GEM.
#' ntry: Number of permutation
permutationTest <- function(Obs.cor, Gs, GEM, Factor, ntry=500){
  set.seed(878712)
  nGns <- nrow(GEM)#number of genes.
  Gs.size <- length(Gs)
  #ncd: Null correlations distribution
  ncd <- c()
  #p_val <- 0
  for(i in 1:ntry){
    ncd[i] <- cor(apply(GEM[sample(nGns, Gs.size), ], 2, mean), Factor)
  }
  p <- sum(ncd >= Obs.cor[1])/(ntry)
  if(p > 0.5){p_val <- 1-p}else{p_val <- p}
  results <- list(p_val=p_val, Ncd=ncd)
  return(results)
}
