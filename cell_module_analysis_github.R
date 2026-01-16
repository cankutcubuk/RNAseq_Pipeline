############################################################################
#             Relative quantification of cell populations                  # 
#           in bulk RNA-seq and spatial GeoMx datasets of                  #
#   https://www.biorxiv.org/content/10.64898/2026.01.10.698810v1.full.pdf  #
#                                                                          #
#                cankutcubuk [at] {gmail} [dot] {com}                      #
#                               2025				       	               #  
#                https://www.qmul.ac.uk/whri/emr/	                       #
#                     @ EMR-QMUL, London, UK		                       #
############################################################################

rm(list=ls())

library(Seurat)
library(edgeR)

load("./norm.vals.RData") # loads norm.vals = normalized gene expression 
load("./fantom5.RData") # loads modgen = cell-specific gene sets

OA_seurat <- CreateSeuratObject(counts = norm.vals, min.cells = 0, min.features = 0, project = "bulkRNAseq")

for(ct in names(modgen)){
  mysubset <- modgen[ct]
  OA_seurat <- AddModuleScore(OA_seurat, features = mysubset, k=F, ctrl.size = 5, name = paste0("SC_",ct), slot = "counts")  
}

names(OA_seurat@meta.data)[grep("SC_",names(OA_seurat@meta.data))] <- sapply(grep("SC_",names(OA_seurat@meta.data), value = T), function(x) substr(x,1,nchar(x)-1))

# normalize the scores between 0-1
for(sset in paste0("SC_",names(modgen))){
  myscore <- as.vector(OA_seurat[[sset]][,1])
  mymin <- min(myscore)
  mymin <- ifelse(mymin > 0, -1 * mymin, abs(mymin))
  myscore <- myscore + mymin
  mymax <- max(myscore)
  myscore <- myscore/mymax
  OA_seurat[[sset]] <- myscore
}

CellModuleScores <- t(OA_seurat@meta.data[,grep("SC_",names(OA_seurat@meta.data))])
rownames(CellModuleScores) <- gsub("SC_","",rownames(CellModuleScores))

write.table(CellModuleScores, file="./CellModuleScores_OA_Tarmac_July2025.csv", sep = ",")
