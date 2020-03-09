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
(load(file = "data/Lung_6_20190817.Rda"))
(load(file = "output/singler_T_Lung_6_20190817.Rda"))
# if singler didn't find all cell labels`
length(singler$singler[[1]]$SingleR.single$labels) == ncol(object@assays$RNA@data)
if(length(singler$singler[[1]]$SingleR.single$labels) < ncol(object@assays$RNA@data)){
        all.cell = colnames(object);length(all.cell)
        know.cell = rownames(singler$singler[[1]]$SingleR.single$labels);length(know.cell)
        object = object[,know.cell]
}


table(names(singler$singler[[1]]$SingleR.single$labels) == colnames(object))
singler$meta.data$orig.ident = object@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = object@reductions$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = Idents(object) # the Seurat clusters (if 'clusters' not provided)
save(singler,file="output/singler_T_Lung_6_20190817.Rda")

##############################
# add singleR label to Seurat
##############################
singlerDF = data.frame("singler1sub" = singler$singler[[1]]$SingleR.single$labels,
                       "singler1main" = singler$singler[[1]]$SingleR.single.main$labels,
                       "orig.ident" = object@meta.data$orig.ident,
                       row.names = names(singler$singler[[1]]$SingleR.single$labels))

table(rownames(singlerDF) %in% colnames(object))
head(singlerDF)
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
singlerDF = cbind(singlerDF,object@reductions$umap@cell.embeddings)
T_cells <- singlerDF$UMAP_1 < 0 & singlerDF$UMAP_2 >-5
other <- paste0("Endothelial_cells|Erythrocytes|Fibroblasts_senescent|Granulocytes|",
       "Hepatocytes|Macrophages_activated|Mast_Cells|",
       "aNSCs|B_cells|Stem_Cells|Stromal_Cells|test")
singlerDF[grepl(other,singlerDF$singler1sub) & T_cells,
           "singler1sub"] = "T_cells:CD4+"
singlerDF = singlerDF[,-c(4:5)]
#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singlerDF$singler1main, singlerDF$orig.ident)) %>%
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
object <- AddMetaColor(object = object, label= "singler1sub", 
                       colors = base::sample(Singler.colors,size = l))
object <- AddMetaColor(object = object, label= "singler1main", colors = singler.colors)
Idents(object) <- "singler1sub"
object %<>% sortIdent()
UMAPPlot.1(object, cols = ExtractMetaColor(object),
           label = T,pt.size = 1,no.legend = T,label.repel = T,
         label.size = 4, repel = T,do.return= T,do.print = T,alpha = 0.9,
         title = "All cell types identified by Mouse RNA-seq reference database")

save(object,file="data/Lung_6_20190817.Rda")

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
