--
title: "Determine Direction of Information Flow FROM and TO MM"
author: "A. Divaraniya"
date: "February 24, 2015"
output: html_document
---
######################################################################
#Load Magic Module Genes
magic_module_genes <- read.delim("~/Google Drive/PhD Research/Joel Dudley/Research Files/Network Analysis/magic_module_genes_updated.txt", na.strings="", stringsAsFactors=FALSE)
magic_module_genes$EntrezID <- as.character(magic_module_genes$EntrezID)

######################################################################
#Generate omental network
omental_network_genes <- read.csv("~/Google Drive/PhD Research/Joel Dudley/Research Files/Network Analysis/omental_network_genes.csv", na.strings="", stringsAsFactors=FALSE)
omental_network_genes$Direction <- NULL

omental_fat_mapped_nodes <- read.delim("~/Google Drive/PhD Research/Joel Dudley/Research Files/Network Analysis/omental_fat_mapped_nodes.txt", stringsAsFactors=FALSE)

omental_network_genes$From_EntrezID <- omental_fat_mapped_nodes$EntrezID[match(omental_network_genes$From, omental_fat_mapped_nodes$GeneSymbol)]

omental_network_genes$To_EntrezID <- omental_fat_mapped_nodes$EntrezID[match(omental_network_genes$To, omental_fat_mapped_nodes$GeneSymbol)]

omental_network_genes <- na.omit(omental_network_genes) #removes all entries that have NA for either from or to
omental_Gene_IDs_Index <- omental_network_genes

omental_network_genes$From <- NULL
omental_network_genes$To <- NULL

KeyDrivers <- read.csv("~/Google Drive/PhD Research/Joel Dudley/Research Files/Network Analysis/KeyDrivers.csv")

relations_whole_omental <- data.frame(from = omental_network_genes$From_EntrezID, to = omental_network_genes$To_EntrezID)
g_whole_omental <- graph.data.frame(omental_network_genes, directed = TRUE)
vertices_whole_omental <- V(g_whole_omental)
network_genes_whole_omental <- vertices_whole_omental$name
length(which(!is.na(match(network_genes_whole_omental, magic_module_genes$EntrezID)))) #458 magic genes in omental network (with EntrezIds)

######################################################################
#Generate residual omental network
RO_network_genes <- V(g_whole_omental)$name[which(!(V(g_whole_omental)$name %in% magic_module_genes$EntrezID))]

RO_network_pos <- match(RO_network_genes, V(g_whole_omental)$name)

g_RO <- induced.subgraph(g_whole_omental, RO_network_pos, impl = "auto")

length(V(g_RO)$name) #4847

######################################################################
#Using Function Information_Flow_Btwn_Omental_MM to get X degree neighbors of each gene within each B subnetwork (Function is saved in Functions.Rmd file)
mode_magic_list <- c("in", "out")
mode_neighbors_list <- c("in", "out", "all")

for (num_neighbors in 1:5){
  for (mode_magic_index in 1:2){
    mode_magic <- mode_magic_list[mode_magic_index]
    for (mode_neighbors_index in 1:3){
      mode_neighbors <- mode_neighbors_list[mode_neighbors_index]
      temp_df <- Information_Flow_Btwn_Omental_MM(num_neighbors, mode_magic, mode_neighbors)
      
      #Assign relevant title for each data frame
      direction_magic <- ifelse(mode_magic == "out", "FROM", ifelse(mode_magic == "in", "TO", "ALL"))
      direction_neighbor <- ifelse(mode_neighbors == "out", "FROM", ifelse(mode_neighbors == "in", "TO", "ALL"))
      neighbors <- c("C", "D", "E", "F", "G")
      
      colnames(temp_df) <- paste(direction_neighbor, "_", colnames(temp_df), sep = "")
      
      assign(paste("B_", direction_neighbor, "_", neighbors[num_neighbors], "_", direction_magic, "_MM", sep = ""), temp_df)
      rm(temp_df)
    }
  }
}
######################################################################
#Identify Gene Symbol for gene of interest

omental_fat_mapped_nodes$GeneSymbol[which(omental_fat_mapped_nodes$EntrezID == "652106")]

######################################################################
#Using Function Information_Flow_Btwn_Omental_Permuted_MM to get X degree neighbors of each gene within each permuted B subnetwork (Function is saved in Functions.Rmd file)
numSamples <- 100
size_MM <- length(which(!is.na(match(magic_module_genes$EntrezID, V(g_whole_omental)$name))))
mode_magic_list <- c("in", "out")
mode_neighbors_list <- c("in", "out")

for (num_neighbors in 1:5){
  for (mode_magic_index in 1:2){
    mode_magic <- mode_magic_list[mode_magic_index]
    for (mode_neighbors_index in 1:2){
      mode_neighbors <- mode_neighbors_list[mode_neighbors_index]
      temp_df <- data.frame(matrix(NA, nrow = 1, ncol = 9))
      colnames(temp_df) <- c("# genes in SN", "SN to MM (OSN)",  "SN to RO (OSN)",  "SN to MM/RO (OSN)",  "SN to RO(ISN)", "SN from MM (OSN)", "SN from RO (OSN)",  "SN from MM/RO(OSN)", "SN from RO (ISN)")
      
      sample_table <- Information_Flow_Btwn_Omental_Permuted_MM(temp_df, numSamples, num_neighbors, mode_magic, mode_neighbors)
      
      #Assign relevant title for each data frame
      direction_magic <- ifelse(mode_magic == "out", "FROM", ifelse(mode_magic == "in", "TO", "ALL"))
      direction_neighbor <- ifelse(mode_neighbors == "out", "FROM", ifelse(mode_neighbors == "in", "TO", "ALL"))
      neighbors <- c("C", "D", "E", "F", "G")
      sample_table <- sample_table[-1,]
      assign(paste("Perm_B_", direction_neighbor, "_", neighbors[num_neighbors], "_", direction_magic, "_MM", sep = ""), sample_table)
      rm(temp_df)
    }
  }
}
######################################################################
#cbind all columns
aB_C_FROM_MM <- cbind(B_FROM_C_FROM_MM, B_TO_C_FROM_MM, B_ALL_C_FROM_MM)
aB_D_FROM_MM <- cbind(B_FROM_D_FROM_MM, B_TO_D_FROM_MM, B_ALL_D_FROM_MM)
aB_E_FROM_MM <- cbind(B_FROM_E_FROM_MM, B_TO_E_FROM_MM, B_ALL_E_FROM_MM)
aB_F_FROM_MM <- cbind(B_FROM_F_FROM_MM, B_TO_F_FROM_MM, B_ALL_F_FROM_MM)
aB_G_FROM_MM <- cbind(B_FROM_G_FROM_MM, B_TO_G_FROM_MM, B_ALL_G_FROM_MM)

aB_C_TO_MM <- cbind(B_FROM_C_TO_MM, B_TO_C_TO_MM, B_ALL_C_TO_MM)
aB_D_TO_MM <- cbind(B_FROM_D_TO_MM, B_TO_D_TO_MM, B_ALL_D_TO_MM)
aB_E_TO_MM <- cbind(B_FROM_E_TO_MM, B_TO_E_TO_MM, B_ALL_E_TO_MM)
aB_F_TO_MM <- cbind(B_FROM_F_TO_MM, B_TO_F_TO_MM, B_ALL_F_TO_MM)
aB_G_TO_MM <- cbind(B_FROM_G_TO_MM, B_TO_G_TO_MM, B_ALL_G_TO_MM)





```{r}
sample_output1 <- sample_output

a <- rbind(sample_output, sample_output1)
```

