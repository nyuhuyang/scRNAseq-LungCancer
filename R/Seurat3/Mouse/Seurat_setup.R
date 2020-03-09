########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(cowplot)
library(kableExtra)
library(magrittr)
library(harmony)
library(scran)
source("../R/Seurat3_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the mouse.eyes dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
df_samples <- readxl::read_excel("doc/20190529_scRNAseq_info.xlsx")
print(df_samples)
df_samples = as.data.frame(df_samples)
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c(paste0("test",3)))
df_samples <- df_samples[sample_n,]
df_samples
samples = df_samples$sample


#======1.2 load  SingleCellExperiment =========================
(load(file = "data/sce_6_20190817.Rda"))
names(sce_list)
object_list <- lapply(sce_list, as.Seurat)
(samples <- sapply(object_list,function(x) unique(x@meta.data$orig.ident)))

for(i in 1:length(samples)){
    object_list[[i]]@meta.data$tests <- df_samples$tests[i]
    object_list[[i]]@meta.data$conditions <- df_samples$conditions[i]
    object_list[[i]]@meta.data$groups <- df_samples$group[i]
    object_list[[i]]@meta.data$tsne <- df_samples$tsne[i]
    
}
#========1.3 merge ===================================
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), object_list)
object@assays$RNA@data = object@assays$RNA@data *log(2)
remove(sce_list,object_list);GC()

(remove <- which(colnames(object@meta.data) %in% "ident"))
meta.data = object@meta.data[,-remove]
object@meta.data = meta.data 
Idents(object) = "orig.ident"

#======1.2 QC, pre-processing and normalizing the data=========================
Idents(object) = factor(Idents(object),levels = df_samples$sample)
(load(file = "output/20190617/g1_6_20190817.Rda"))

#object %<>% subset(subset = nFeature_RNA > 100 & nCount_RNA > 200 & percent.mt < 15)
# FilterCellsgenerate Vlnplot before and after filteration
g2 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
    VlnPlot(object = object, features = features, ncol = 3, pt.size = 0.01)+
        theme(axis.text.x = element_text(size=15),legend.position="none")
})

save(g2,file= paste0(path,"g2_6_20190817.Rda"))
jpeg(paste0(path,"S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nFeature_RNA before filteration")+
                    scale_y_log10(limits = c(100,10000)),
                g2[[1]]+ggtitle("nFeature_RNA after filteration")+
                    scale_y_log10(limits = c(100,10000))))
dev.off()
jpeg(paste0(path,"S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nCount_RNA before filteration")+
                    scale_y_log10(limits = c(500,100000)),
                g2[[2]]+ggtitle("nCount_RNA after filteration")+ 
                    scale_y_log10(limits = c(500,100000))))
dev.off()
jpeg(paste0(path,"S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % before filteration")+
                    ylim(c(0,50)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                    ylim(c(0,50))))
dev.off()

######################################

# After removing unwanted cells from the dataset, the next step is to normalize the data.
#object <- NormalizeData(object = object, normalization.method = "LogNormalize", 
#                      scale.factor = 10000)
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(object), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
jpeg(paste0(path,"VariableFeaturePlot.jpeg"), units="in", width=10, height=7,res=600)
print(plot2)
dev.off()
hvf.info <- HVFInfo(object = object)
hvf.info = hvf.info[VariableFeatures(object),]
write.csv(hvf.info, file = paste0(path,"high_variable_genes.csv"))
#======1.3 1st run of pca-tsne  =========================
object %<>% ScaleData(features = VariableFeatures(object))
object %<>% RunPCA(features = VariableFeatures(object),verbose =F,npcs = 100)
#object %<>% JackStraw(num.replicate = 20,dims = 100)
#object %<>% ScoreJackStraw(dims = 1:85)

jpeg(paste0(path,"ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = 85)+
    ggtitle("ElbowPlot")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18)) 
dev.off()

jpeg(paste0(path,"JackStrawPlot.jpeg"), units="in", width=10, height=7,res=600)
JackStrawPlot(object, dims = 65:80)+
    ggtitle("JackStrawPlot")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18)) 
dev.off()

npcs =50
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(reduction = "pca",resolution = 0.6,
                       dims.use = 1:npcs,print.output = FALSE)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)

p0 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 no.legend = T,label.size = 4, repel = T, title = "Original")
p1 <- UMAPPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 no.legend = T,label.size = 4, repel = T, title = "Original")

#======1.4 Performing SCTransform and integration =========================
set.seed(100)
object_list <- SplitObject(object, split.by = "orig.ident")
object_list %<>% lapply(SCTransform)
object.features <- SelectIntegrationFeatures(object_list, nfeatures = 3000)
options(future.globals.maxSize= object.size(object_list)*1.5)
object_list <- PrepSCTIntegration(object.list = object_list, anchor.features = object.features, 
                                  verbose = FALSE)
anchors <- FindIntegrationAnchors(object_list, normalization.method = "SCT", 
                                  anchor.features = object.features)
object1 <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

remove(object.anchors,object_list);GC()
object %<>% RunPCA(npcs = 100, verbose = FALSE)
object <- JackStraw(object, num.replicate = 20,dims = 100)
object <- ScoreJackStraw(object, dims = 1:100)
jpeg(paste0(path,"JackStrawPlot_SCT.jpeg"), units="in", width=10, height=7,res=600)
JackStrawPlot(object, dims = 90:100)
dev.off()
npcs = 100
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(reduction = "pca",resolution = 0.6,
                         dims.use = 1:npcs,print.output = FALSE)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)

p2 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 label.size = 4, repel = T,title = "Intergrated tSNE plot")

p3 <- UMAPPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 label.size = 4, repel = T,title = "Intergrated UMAP plot")

#=======1.9 summary =======================================
jpeg(paste0(path,"S1_TSNEPlot_batch.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p0+ggtitle("Clustering without integration")+
              theme(plot.title = element_text(hjust = 0.5,size = 18)),
          p2+ggtitle("Clustering with integration")+
              theme(plot.title = element_text(hjust = 0.5,size = 18)))
dev.off()

jpeg(paste0(path,"S1_UMAPPlot_batch.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1+ggtitle("Clustering without integration")+
              theme(plot.title = element_text(hjust = 0.5,size = 18)),
          p3+ggtitle("Clustering with integration")+
              theme(plot.title = element_text(hjust = 0.5,size = 18)))
dev.off()

UMAPPlot.1(object, group.by="integrated_snn_res.0.6",pt.size = 1,label = T,no.legend = T,
           label.repel = T, alpha = 1,border = T,do.print = T,
           label.size = 4, repel = T,title = "All cluster in UMAP plot resolution = 0.6")

UMAPPlot.1(object, group.by="integrated_snn_res.0.6",split.by = "conditions",
           pt.size = 1,label = T,no.legend = T,
           label.repel = T, alpha = 1,border = T,do.print = T,
           label.size = 4, repel = T,title = "All cluster in UMAP plot resolution = 0.6")

p3 <- TSNEPlot.1(object, group.by="integrated_snn_res.0.6",pt.size = 1,label = T,no.legend = T,
                 label.repel = T, alpha = 1,border = T,
                 label.size = 4, repel = T,title = "Total Clusters",do.print = T)
p4 <- TSNEPlot.1(object, group.by="integrated_snn_res.0.6",pt.size = 1,label = T,no.legend = T,
                 label.repel = T, alpha = 1,border = T, split.by = "conditions",
                 label.size = 4, repel = T,title = NULL, do.print = T)

jpeg(paste0(path,"S1_split_TSNEPlot_all.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p3, p4, align = "h")
dev.off()
object@assays$RNA@scale.data = matrix(0,0,0)
object@assays$integrated@scale.data = matrix(0,0,0)
save(object, file = "data/Lung_6_20190817.Rda")

object_data = object@assays$SCT@data
save(object_data, file = "data/object.data_Lung_6_20190817.Rda")
