install.packages("Seurat")
library(Seurat)
library(tidyverse)
# List of file paths
file_paths <- c("~/workspace/data/cmlnf_B1_dense.csv",
                "~/workspace/data/cmlnf_B2_dense.csv",
                "~/workspace/data/cmlnf_B5_dense.csv",
                "~/workspace/data/cmlnf_B6_dense.csv",
                "~/workspace/data/cmlnf_B9_dense.csv",
                "~/workspace/data/cmlnf_B12_dense.csv",
                "~/workspace/data/cmlnf_B21_dense.csv",
                "~/workspace/data/cmlnf_B22_dense.csv",
                "~/workspace/data/cmlnf_B23_dense.csv",
                "~/workspace/data/cmlnf_B24_dense.csv",
                "~/workspace/data/cmlnf_B25_dense.csv",
                "~/workspace/data/cmlnf_B26_dense.csv")

# Loop through the file paths
for (i in 1:length(file_paths)) {
  
  # Read in the CSV file
  cnt <- read.csv(file_paths[i], header=T, row.names = 1)
  
  # Transpose the data and randomly select 1000 rows
  cnt <- t(cnt[sample(nrow(cnt), 1000), ])
  
  # Create a Seurat object
  seu_r <- CreateSeuratObject(counts = cnt, project = paste0("B", i), min.cells = 3, min.features = 200)
  
  # Assign the Outcome variable in the meta.data slot
  seu_r@meta.data$Outcome <- c("Responder")
  
  # Print the head of the meta.data
  head(seu_r@meta.data)
  
  # Remove the count object from memory
  rm(cnt)
  
  # Run garbage collection
  gc()
}


# List of file paths
file_paths <- c("~/workspace/data/cmlnf_B3_dense.csv",
                "~/workspace/data/cmlnf_B4_dense.csv",
                "~/workspace/data/cmlnf_B7_dense.csv",
                "~/workspace/data/cmlnf_B8_dense.csv",
                "~/workspace/data/cmlnf_B17_dense.csv",
                "~/workspace/data/cmlnf_B18_dense.csv",
                "~/workspace/data/cmlnf_B19_dense.csv",
                "~/workspace/data/cmlnf_B20_dense.csv",
                "~/workspace/data/cmlnf_B27_dense.csv",
                "~/workspace/data/cmlnf_B28_dense.csv")

# Loop through the file paths
for (i in 1:length(file_paths)) {
  
  # Read in the CSV file
  cnt <- read.csv(file_paths[i], header=T, row.names = 1)
  
  # Transpose the data and randomly select 1000 rows
  cnt <- t(cnt[sample(nrow(cnt), 1000), ])
  
  # Create a Seurat object
  seu_nr <- CreateSeuratObject(counts = cnt, project = paste0("B", i), min.cells = 3, min.features = 200)
  
  # Assign the Outcome variable in the meta.data slot
  seu_nr@meta.data$Outcome <- c("Non-Responder")
  
  # Print the head of the meta.data
  head(seu_nr@meta.data)
  
  # Remove the count object from memory
  rm(cnt)
  
  # Run garbage collection
  gc()
}



#merge datasets
#seu <- merge(seu_r, seu_nr, add.cell.ids = c("B01", "B02", "B03", "B04", fill with sample names ), project = "DLIstudy")

