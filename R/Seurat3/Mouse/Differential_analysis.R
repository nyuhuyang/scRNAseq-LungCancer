########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(gplots)
library(MAST)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file = "data/Lung_harmony_6_2019728.Rda"))
Idents(object) <-"singler1main"
object <- sortIdent(object)
#TSNEPlot(object)
Cancer.markers <- FindAllMarkers.UMI(object = object, only.pos = F, logfc.threshold = 0.5,
                                        test.use = "MAST")
write.csv(Cancer.markers,paste0(path,"Cancer_cell_type_markers.csv"))

#DoHeatmap.1======
Top_n = 10
top <-  Cancer.markers %>% group_by(cluster) %>% 
        top_n(Top_n, avg_logFC) %>% as.data.frame()
add.genes = unique(c(as.character(top$gene)))
object %<>% ScaleData(features= add.genes)
DoHeatmap.1(object, add.genes = add.genes,
            Top_n = Top_n, do.print=T, angle = 45,
            group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0.02, label=T, cex.row= 400/length(add.genes), 
            #legend.size = 0,
            width=10, height=6.5,unique.name = T, pal_gsea = T,
            title = paste("Top",Top_n,"DE genes in all cell types"))
GC()
