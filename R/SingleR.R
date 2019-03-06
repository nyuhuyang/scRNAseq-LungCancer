########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(sapply(c("Seurat","magrittr","SingleR","dplyr"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
args <- commandArgs(trailingOnly = TRUE)
#(load(file = "data/LungCancer_4_20190305.Rda"))
args[1] = as.character(args[1])
(load(file = args[1]))
print(species <- CheckSpecies(object))
singler = CreateSinglerObject(as.matrix(object@data), annot = NULL, 
                              project.name = "LungCancer",
                              min.genes = 500,technology = "10X", species = species, citation = "",
                              ref.list = list(),normalize.gene.length = F, variable.genes = "de",
                              fine.tune = T, do.signatures = F, clusters = NULL,
                              numCores = SingleR.numCores/7*3)

# if singler didn't find all cell labels
if(length(singler$singler[[1]]$SingleR.single$labels) != ncol(object@data)){
        all.cell = object@cell.names
        print(length(all.cell))
        know.cell = rownames(singler$singler[[1]]$SingleR.single$labels)
        print(length(know.cell))
        object = SubsetData(object, cells.use = know.cell)
}
singler$meta.data$orig.ident = object@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = object@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = object@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file=paste0("output/singler_T_",basename(args[1])))