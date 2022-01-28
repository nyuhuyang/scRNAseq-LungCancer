library(Seurat)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

(load("data/Lung_5_20200313.Rda"))
object = UpdateSeuratObject(object)
TSNEPlot.1(object, group.by = "singler1sub",do.print = F)
