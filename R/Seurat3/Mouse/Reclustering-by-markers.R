invisible(lapply(c("Seurat","dplyr","kableExtra","ggplot2","cowplot","sctransform",
                   "harmony"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("data")) dir.create("data")
if(!dir.exists("doc")) dir.create("doc")
TitleCenter
# 3.1.1 load data
# Rename ident
(load(file = "data/Lung_5_CCA_20200313.Rda"))
DefaultAssay(object) = "SCT"

# read Classifier genes 
df_genes <- readxl::read_excel("doc/2020_3_26_Single_Cell_Groups.xlsx")
genes_list <- df2list(df_genes)
genes <- FilterGenes(object, unlist(genes_list), unique = T)

VariableFeatures(object) = genes
object@assays$SCT@scale.data = matrix(0,0,0)
object %<>% ScaleData()
object@reductions$pca = NULL
object %<>% RunPCA(npcs = 100, verbose = FALSE)
npcs = 100
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(reduction = "pca",resolution = 0.6,
                         dims.use = 1:npcs,print.output = FALSE)
object@reductions$tsne = NULL
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object@reductions$umap = NULL
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
TSNEPlot.1(object, label = T, label.repel = T,do.print = T, no.legend = T)
