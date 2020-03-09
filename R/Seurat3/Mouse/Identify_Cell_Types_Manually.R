library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(readxl)
library(cowplot)
source("../R/Seurat3_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
marker_path <- paste0(path,"markers/")
if(!dir.exists(marker_path))dir.create(marker_path, recursive = T)
#====== 2.1 Identify cell types ==========================================
(load(file="data/Lung_6_20190817.Rda"))
DefaultAssay(object) = "RNA"

df_markers <- readxl::read_excel("../seurat_resources/bio-rad-markers.xlsx",sheet = "Human.sub")
colnames(df_markers) = gsub(" ","_",colnames(df_markers))
colnames(df_markers) = gsub(":|\\/","_",colnames(df_markers))
colnames(df_markers) = gsub("\\+","",colnames(df_markers))
df_markers = df_markers[,grep("Alias",colnames(df_markers),invert = T)]
marker.list <- df2list(df_markers)

marker.list %<>% lapply(function(x) x[1:18]) %>% 
    lapply(function(x) FilterGenes(object,x)) %>% 
    lapply(function(x) x[!is.na(x)]) %>% 
    lapply(function(x) x[1:min(length(x),12)])
marker.list <- marker.list[!is.na(marker.list)]
marker.list %>% list2df %>% t %>% kable() %>% kable_styling()

marker.list <- marker.list[sapply(marker.list,length)>0]
marker.list %>% list2df %>% t %>% kable() %>% kable_styling()
Idents(object) <- "integrated_snn_res.0.6"

for(i in 1:length(marker.list)){
    p <- lapply(marker.list[[i]], function(marker) {
        FeaturePlot(object = object, features = marker,pt.size = 0.5, label=T)+
            NoLegend()+
            ggtitle(paste0(marker,Alias(df = df_markers,gene = marker)))+
            theme(plot.title = element_text(hjust = 0.5,size = 15,face = "plain"))
    })
    jpeg(paste0(marker_path,names(marker.list)[i],".jpeg"),units="in", width=10, height=7,res=600)
    print(do.call(plot_grid, p)+ ggtitle(paste(names(marker.list)[i],"markers"))+
              theme(plot.title = element_text(hjust = 0.5,size = 20)))
    dev.off()
    print(paste0(i,":",length(marker.list)))
}

#======== rename ident =================
object$cell.type <- plyr::mapvalues(object$RNA_snn_res.0.6,
                                            from = 0:19,
                                            to = c("Alveolar macrophages",
                                                   "Alveolar type I&II cells \nDistal secretory cells",
                                                   "Endothelial cells",
                                                   "Macrophages",
                                                   "T cells",
                                                   "Stromal fibroblasts",
                                                   "Monocytes",
                                                   "Ciliated cells",
                                                   "Smooth muscle cells",
                                                   "Secretory cells",
                                                   "Endothelial cells",
                                                   "Alveolar type I&II cells \nDistal secretory cells",
                                                   "B cells",
                                                   "Red blood cells",
                                                   "Basal cells",
                                                   "HSC/progenitor cells",
                                                   "Lymphatic endothelial cells",
                                                   "Endothelial cells",
                                                   "Chondrocytes & Fibroblast",
                                                   "Alveolar type I&II cells \nDistal secretory cells"))
object$cell.type <- as.character(object$cell.type)

# rename NK cells
(load(file="output/singlerF_Lung_12_20190614.Rda"))
singlerDF = data.frame("singler1sub" = singler$singler[[1]]$SingleR.single$labels,
                       row.names = rownames(singler$singler[[1]]$SingleR.single$labels))
head(singlerDF)
NK_cells <-  rownames(singlerDF)[singlerDF$singler1sub %in% "NK_cells"]
# PECAM1_neg
SELE <- WhichCells(object, expression = SELE >0 )


object@meta.data[NK_cells,"cell.type"] = "NK cells"
object@meta.data[SELE,"cell.type"] = "Endothelial PECAM1-neg"

Idents(object) <- "cell.type"
TSNEPlot.1(object, cols = ExtractMetaColor(object),label = T,pt.size = 1,
           label.size = 5, repel = T,no.legend = F,do.print = F,
           title = "Cell types")
##############################
# process color scheme
##############################
object <- sortIdent(object)
table(Idents(object)) %>% kable %>% kable_styling()
singler_colors <- readxl::read_excel("doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors1[duplicated(singler_colors1)]
length(singler_colors1)
length(unique(object$cell.type))
object <- AddMetaColor(object = object, label= "cell.type", colors = singler_colors1)
Idents(object) <- "cell.type"

object <- sortIdent(object)

TSNEPlot.1(object, cols = ExtractMetaColor(object),label = T,pt.size = 1,
           label.size = 5, repel = T,no.legend = F,do.print = F,
           title = "Cell types")
save(object,file="data/Lung_harmony_12_20190614.Rda")

meta.data = cbind.data.frame(object@meta.data,object@reductions$tsne@cell.embeddings)
colnames(meta.data)
meta.data = meta.data[,c("tSNE_1","tSNE_2","RNA_snn_res.1.2","cell.type")]
meta.data$RNA_snn_res.1.2 = as.numeric(as.character(meta.data$RNA_snn_res.1.2))

meta.data$cell.type = gsub(" \n",", ",meta.data$cell.type)
table(meta.data$cell.type)
meta.data = meta.data[order(meta.data$RNA_snn_res.1.2),]
write.csv(meta.data, paste0(path,"tSNE_coordinates.csv"))

data = as.matrix(t(object@assays$RNA@data))
tsne_data = cbind(meta.data[,3:4], data[match(rownames(meta.data),
                                        rownames(data)), ])

tsne_data[1:4,1:5]
write.csv(t(tsne_data), paste0(path,"Lung_exp.csv"))

#====== 2.2 Identify Epi cell types ==========================================
(load(file="data/Epi_harmony_12_20190627.Rda"))

df_markers <- readxl::read_excel("doc/Renat.markers.xlsx",sheet = "20190613")

#df_markers <- readxl::read_excel("../seurat_resources/bio-rad-markers.xlsx",sheet = "Human.sub")
colnames(df_markers) = gsub(" ","_",colnames(df_markers))
colnames(df_markers) = gsub(":|\\/","_",colnames(df_markers))
colnames(df_markers) = gsub("\\+","",colnames(df_markers))
#markers = df_markers[,-grep("Alias",colnames(df_markers))]
marker.list <- df2list(df_markers)

marker.list %<>% lapply(function(x) x) %>% 
    lapply(function(x) FilterGenes(Epi,x)) %>% 
    lapply(function(x) x[!is.na(x)])
#marker.list %>% list2df %>% t %>% kable() %>% kable_styling()
Idents(Epi) <- "RNA_snn_res.0.3"

for(i in 1:length(marker.list)){
    p <- lapply(marker.list[[i]], function(marker) {
        FeaturePlot(object = Epi, features = marker,pt.size = 0.5, label=T)+
            NoLegend()+
            ggtitle(paste0(marker,Alias(df = df_markers,gene = marker)))+
            theme(plot.title = element_text(hjust = 0.5,size = 15,face = "plain"))
    })
    jpeg(paste0(path,names(marker.list)[i],".jpeg"),units="in", width=10, height=7,res=600)
    print(do.call(cowplot::plot_grid, p)+ ggtitle(paste(names(marker.list)[i],"markers"))+
              theme(plot.title = element_text(hjust = 0.5,size = 20)))
    dev.off()
    print(paste0(i,":",length(marker.list)))
}

#======== rename ident =================
Epi$cell.type <- plyr::mapvalues(Epi$RNA_snn_res.0.8,
                                    from = 1:9,
                                    to = c("Alveolar type II cells",
                                           "M-ciliated cells",
                                           "Airway secretory cells",
                                           "D-ciliated cells",
                                           "Alveolar type I cells",
                                           "Distal airway secretory cells",
                                           "Airway basal stem cells",
                                           "Unknown Alveolar cells",
                                           "Alveolar type II cells"))
Idents(Epi) = "cell.type"
Epi %<>% sortIdent()
table(Idents(Epi)) %>% kable() %>% kable_styling()
# process color scheme
singler_colors <- readxl::read_excel("doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors1[duplicated(singler_colors1)]
length(singler_colors1)
Epi <- AddMetaColor(object = Epi, label= "cell.type", colors = singler_colors1)
Epi %<>% sortIdent()

TSNEPlot.1(Epi, cols = ExtractMetaColor(Epi),label = F,pt.size = 1,
           label.size = 5, repel = T,no.legend = F,do.print = T,
           title = "Epithelial cell types")

p5 <- UMAPPlot(Epi, group.by="cell.type",pt.size = 1,label = F,
               cols = ExtractMetaColor(Epi),
               label.size = 4, repel = T)+ggtitle("Epithelial cell types")+
    theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))

jpeg(paste0(path,"Epi_umap_cell.type.jpeg"), units="in", width=10, height=7,res=600)
print(p5)
dev.off()
save(Epi, file = "data/Epi_harmony_12_20190627.Rda")
