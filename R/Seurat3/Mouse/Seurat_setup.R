########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

invisible(lapply(c("Seurat","dplyr","kableExtra","ggplot2","cowplot","sctransform",
                   "harmony"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("data")) dir.create("data")
if(!dir.exists("doc")) dir.create("doc")
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the mouse.eyes dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
df_samples <- readxl::read_excel("doc/20200313_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% paste0("test",3))
df_samples = df_samples[sample_n,]
dim(df_samples)
head(df_samples)
(attach(df_samples))
(samples = df_samples$sample)

#======1.2 load  SingleCellExperiment =========================
(load(file = "data/sce_5_20200313.Rda"))
names(sce_list)
object_list <- lapply(sce_list, as.Seurat)
(samples <- sapply(object_list,function(x) unique(x@meta.data$orig.ident)))

for(i in 1:length(samples)){
    object_list[[i]]@meta.data$conditions <- df_samples$conditions[i]
    object_list[[i]]@meta.data$groups <- df_samples$group[i]
}
#========1.3 merge ===================================
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), object_list)
object@assays$RNA@data = object@assays$RNA@data *log(2)
remove(sce_list,object_list);GC()

(remove <- which(colnames(object@meta.data) %in% "ident"))
meta.data = object@meta.data[,-remove]
object@meta.data = meta.data 
Idents(object) = "orig.ident"

######################################

# After removing unwanted cells from the dataset, the next step is to normalize the data.
#object <- NormalizeData(object = object, normalization.method = "LogNormalize", 
#                      scale.factor = 10000)
DefaultAssay(object) = "SCT"
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
object <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

remove(object.anchors,object_list);GC()
object %<>% ScaleData()
object %<>% RunPCA(npcs = 100, verbose = FALSE)
object <- JackStraw(object, num.replicate = 20,dims = 100)
object <- ScoreJackStraw(object, dims = 1:100)
a <- seq(1,100, by = 10)
b <- a+9
for(i in seq_along(a)){
        jpeg(paste0(path,"JackStrawPlot_cca_",i,"_", a[i],"_",min(b[i],100),".jpeg"),
             units="in", width=10, height=7,res=600)
        print(JackStrawPlot(object, dims = a[i]:min(b[i],100)))
        Progress(i,length(a))
        dev.off()
}
npcs = 100
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(reduction = "pca",resolution = 0.6,
                         dims.use = 1:npcs,print.output = FALSE)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)

Idents(object) = "orig.ident"
p1 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 do.print = T,
                 label.size = 4, repel = T,title = "CCA Intergrated tSNE plot")

p2 <- UMAPPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 do.print = T,
                 label.size = 4, repel = T,title = "CCA Intergrated UMAP plot")
format(object.size(object), unit = "GB")
save(object, file = "data/Lung_5_CCA_20200313.Rda")

#======1.8  Harmony =========================
DefaultAssay(object)  = "SCT"
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20, nfeatures = 2000,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
object %<>% ScaleData
object %<>% RunPCA(verbose = T,npcs = 100)

object %<>% JackStraw(num.replicate = 20,dims = 100)
object %<>% ScoreJackStraw(dims = 1:100)
a <- seq(1,100, by = 10)
b <- a+9
for(i in seq_along(a)){
        jpeg(paste0(path,"JackStrawPlot_SCT_",i,"_",a[i],"_",min(b[i],100),".jpeg"), units="in", width=10, height=7,res=600)
        print(JackStrawPlot(object, dims = a[i]:min(b[i],100)))
        Progress(i,length(a))
        dev.off()
}

npcs = 74
jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                     theta = 2, plot_convergence = TRUE,
                                     nclust = 50, max.iter.cluster = 100))
dev.off()

object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
system.time(object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs))
object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
Idents(object) = "orig.ident"
p3 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,legend.size = 15,
                 do.print = T,
                 label.size = 4, repel = T,title = "Harmony intergrated tSNE plot")
p4 <- UMAPPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,legend.size = 15,
                 do.print = T,
                 label.size = 4, repel = T,title = "Harmony intergrated UMAP plot")
# =====

object@assays$integrated@scale.data = matrix(0,0,0)
save(object, file = "data/Lung_5_Harmony_20200313.Rda")

object_data <- object@assays$SCT@data
save(object_data, file = "data/Lung.data_5_20200313.Rda")
