library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/Lung_5_CCA_20200313.Rda"))
(load(file = "output/singler_T_Lung_5_20200313.Rda"))
# if singler didn't find all cell labels`
length(singler$singler[[1]]$SingleR.single$labels) == ncol(object@assays$RNA@data)
if(length(singler$singler[[1]]$SingleR.single$labels) < ncol(object@assays$RNA@data)){
        all.cell = colnames(object);length(all.cell)
        know.cell = rownames(singler$singler[[1]]$SingleR.single$labels);length(know.cell)
        object = object[,know.cell]
}
table(rownames(singler$singler[[1]]$SingleR.single$labels) == colnames(object))
all(rownames(singler$singler[[1]]$SingleR.single$labels) %in% colnames(object))

##############################
# add singleR label to Seurat
##############################
singlerDF1 = data.frame("singler1sub" = singler$singler[[1]]$SingleR.single$labels,
                       "singler1main" = singler$singler[[1]]$SingleR.single.main$labels,
                       "orig.ident" = gsub("\\_.*","",rownames(singler$singler[[1]]$SingleR.single$labels)),
                       row.names = rownames(singler$singler[[1]]$SingleR.single$labels))
singlerDF = apply(singlerDF1,2, as.character) %>% as.data.frame()
rownames(singlerDF) = rownames(singlerDF1)
dim(singlerDF)
head(singlerDF)
table(rownames(singlerDF) %in% colnames(object))
singlerDF = singlerDF[colnames(object),]
table(rownames(singlerDF) == colnames(object))
apply(singlerDF,2,function(x) length(unique(x)))

###############################
# check the spearman correlation
###############################
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"DrawHeatmap_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,normalize = F))
dev.off()
jpeg(paste0(path,"DrawHeatmap_sub1_N.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single,top.n = 50,normalize = T))
dev.off()

# adjust cell label
table(rownames(singlerDF) == rownames(object@reductions$umap@cell.embeddings))
#singlerDF = cbind(singlerDF,object@reductions$umap@cell.embeddings)
#T_cells <- singlerDF$UMAP_1 < 0 & singlerDF$UMAP_2 >-5
#other <- paste0("Endothelial_cells|Erythrocytes|Fibroblasts_senescent|Granulocytes|",
#       "Hepatocytes|Macrophages_activated|Mast_Cells|",
#       "aNSCs|B_cells|Stem_Cells|Stromal_Cells|test")
#singlerDF[grepl(other,singlerDF$singler1sub) & T_cells,
#           "singler1sub"] = "T_cells:CD4+"
#singlerDF = singlerDF[,-c(4:5)]
#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singlerDF$singler1sub, singlerDF$orig.ident)) %>%
        kable_styling()
singlerDF$orig.ident %>% table() %>% kable() %>% kable_styling()
singlerDF$singler1main %>% table() %>% kable() %>% kable_styling()
singlerDF$singler1sub %>% table() %>% prop.table %>% kable() %>% kable_styling()


##############################
# process color scheme
##############################
set.seed(101)
(l = length(unique(singlerDF$singler1sub)))
object <- AddMetaData(object = object,metadata = singlerDF)
Idents(object) = "singler1sub"
object %<>% sortIdent()
object <- AddMetaColor(object = object, label= "singler1sub", 
                       colors = Singler.colors[1:l])
lapply(c(UMAPPlot.1, TSNEPlot.1), function(fun) {
    fun(object, cols = ExtractMetaColor(object),
         label = T,pt.size = 1,no.legend = T,label.repel = T,
         label.size = 4, repel = T,do.return= T,do.print = T,alpha = 0.9,
         title = "All cell types identified by Mouse RNA-seq reference database")
})


save(object,file="data/Lung_5_CCA_20200313.Rda")

##############################
# draw tsne plot
##############################
Idents(object) <- "singler1sub"
object %<>% sortIdent()
UMAPPlot.1(object[,T_cells], group.by="singler1sub",
           cols = ExtractMetaColor(object[,T_cells]),
           split.by = "orig.ident", ncol = 3,border = T,
           label = F,pt.size = 1,no.legend = T,label.repel = F,
           label.size = 4, repel = T,do.return= T,do.print = T,alpha = 0.9,
           title = "Compare T cell types in each sample")
table(as.character(object[,T_cells]$singler1sub), 
      as.character(object[,T_cells]$orig.ident)) %>% kable() %>% kable_styling()
