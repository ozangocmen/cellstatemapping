#' Weighted t-test for cell ratios of patients from two groups.
#' This function utilizes R library 'weights'
#' @param Pid patient ids vector for each data point, i.e. cell.
#' @param Groups groups vector for each data point, i.e. cell.
#' @param Clusters cluster ids vector for each data point, i.e. cell.
#' @param nTry number of re-sampling for bootstrapping
#' @param N_r number of patient samples in one group
#' @param N_nr number of patient samples in the second group
#' @param GroupNames
WeiTtest <- function(Pid, Groups, Clusters, nTry=3000, N_r, N_nr, GroupNames){
  
  seu.meta.data <- data.frame(Pid= Pid, Groups = Groups, Clusters = Clusters)
  
  seu.meta.data %>%
    group_by(Clusters) %>%
    mutate(Uk=n()) %>% #Uk: total number of T cells in the cluster k
    group_by(Groups) %>% 
    mutate(Np=n()) %>% #Np: total number of T cells in group
    group_by(Pid) %>% 
    mutate(Ni=n()) %>% #Ni: total number of T cells of patient i
    group_by(Clusters, Pid) %>% 
    mutate(Nk=n()) %>% #Nk: total number of T cells of patient i in cluster k
    mutate(koi=Nk/Ni, iop=Ni/Np) %>% #koi: proportion of cells, iop: T cell contribution to the group pool by patient i
    mutate(Wi=ifelse(Groups==GroupNames[1], N_r*iop, N_nr*iop)) %>%
    unique() %>% as.data.frame() -> meta.ttest
  
  cls.exclude <- droplevels(as.data.frame(table(meta.ttest$Groups, meta.ttest$Clusters))$Var2[as.data.frame(table(meta.ttest$Groups, meta.ttest$Clusters))$Freq == 0])
  cls_list <- unique(Clusters)[!unique(Clusters) %in% cls.exclude]
  cls_stats <- c()
  examined.cls.list <- c()
  st.list <- list()
  rt.df <- data.frame()
  sm.plist <- list()
  for(k in cls_list){
    
    A <- meta.ttest[which(meta.ttest$Clusters == k & meta.ttest$Groups == GroupNames[1]), "koi"]
    B <- meta.ttest[which(meta.ttest$Clusters == k & meta.ttest$Groups == GroupNames[2]), "koi"]
    wA <- meta.ttest[which(meta.ttest$Clusters == k & meta.ttest$Groups == GroupNames[1]), "Wi"]
    wB <- meta.ttest[which(meta.ttest$Clusters == k & meta.ttest$Groups == GroupNames[2]), "Wi"]
    
    if(length(A) == 1 & length(B) == 1 | length(A) == 1){print(k); cls.exclude <- c(cls.exclude, k); next;}
    
    wtd.t.test(
      A,
      B,
      weight = wA,
      weighty = wB,
      samedata = F,
      alternative = "greater") -> tt
    cls_p_val <- tt$coefficients['p.value']
    
    U_k <- unique(meta.ttest[meta.ttest$Clusters == k,]$Uk)#The size of the cluster k
    pdist <- c()
    for(b in 1:nTry){
      rand.meta <- seu.meta.data
      rand.meta$Clusters <- as.vector(rand.meta$Clusters)
      #Randomly select U_k number of cells from the pool of all samples and compute the p-value using the same weighted t.test.
      cls = 'rand'
      rand.meta[sample(1:nrow(rand.meta), size = U_k, replace = T), 'Clusters'] <- cls
      
      rand.meta %>%
        group_by(Clusters) %>%
        mutate(Uk=n()) %>% #Uk: total number of T cells in the cluster k
        group_by(Groups) %>% 
        mutate(Np=n()) %>% #Np: total number of T cells in group
        group_by(Pid) %>% 
        mutate(Ni=n()) %>% #Ni: total number of T cells of patient i
        group_by(Clusters, Pid) %>% 
        mutate(Nk=n()) %>% #Nk: total number of T cells of patient i in cluster k
        mutate(koi=Nk/Ni, iop=Ni/Np) %>% #koi: proportion of cells, iop: T cell contribution to the group pool by patient i
        mutate(Wi=ifelse(Groups=="Responder", N_r*iop, N_nr*iop)) %>%
        unique() %>% 
        filter(Clusters == cls) %>%
        as.data.frame() -> rand.meta.ttest
      
      a <- rand.meta.ttest[which(rand.meta.ttest$Clusters == cls & rand.meta.ttest$Groups == GroupNames[1]), "koi"]
      b <- rand.meta.ttest[which(rand.meta.ttest$Clusters == cls & rand.meta.ttest$Groups == GroupNames[2]), "koi"]
      wa <- rand.meta.ttest[which(rand.meta.ttest$Clusters == cls & rand.meta.ttest$Groups == GroupNames[1]), "Wi"]
      wb <- rand.meta.ttest[which(rand.meta.ttest$Clusters == cls & rand.meta.ttest$Groups == GroupNames[2]), "Wi"]
      
      if(length(a) == 1 & length(b) == 1){next}
      
      #Compute the p-value on randomly selected cells.
      try(
        wtd.t.test(
          a,
          b,
          weight = wa,
          weighty = wb,
          samedata = F,
          alternative = "greater") -> tt
      )
      pdist <- c(pdist, tt$coefficients['p.value'])
      
    }
    
    cls_Q_val <- table(pdist < cls_p_val)['TRUE']/nTry
    names(cls_Q_val) <- "class_Q_val"
    
    print(paste("cluster:", k, "pval:", cls_p_val, "Qval:", cls_Q_val))
    
    sm.plist[[k]] <- pdist
    st.list[k] <- cls_Q_val
    
    examined.cls.list <- c(examined.cls.list, k)
    
    df <- data.frame(cbind(Cluster=k,
                           rbind(meta.ttest[which(meta.ttest$Clusters == k & meta.ttest$Groups == GroupNames[1]), ],
                                 meta.ttest[which(meta.ttest$Clusters == k & meta.ttest$Groups == GroupNames[2]), ])
    ))
    rt.df <- rbind(rt.df, df)
  }
  
  rt.df$Groups <- factor(rt.df$Groups, c(GroupNames[1], GroupNames[2]))
  
  object <- new(Class = "ClusCompStats",
                stats = st.list,
                ratios= rt.df,
                plist = sm.plist)
  
  return(object)
}
