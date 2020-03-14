########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("R.utils","Seurat","dplyr","kableExtra","ggplot2","scater",
                   "scran","BiocSingular","Matrix","cowplot"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("data")) dir.create("data")
if(!dir.exists("doc")) dir.create("doc")
########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======1.1 Load the data files and Set up Seurat object =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/20200313_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% paste0("test",3))
df_samples = df_samples[sample_n,]
dim(df_samples)
head(df_samples)
(attach(df_samples))
(samples = df_samples$sample)
# check missing data
current <- list.files("data")
(current <- current[!grepl(".Rda|RData",current)])
(missing_data <- df_samples$sample.id[!(df_samples$sample.id %in% current)])

# select species
if(unique(df_samples$species) %in%  c("Homo_sapiens","Human")) species <- "hg19"
if(unique(df_samples$species) %in%  c("Mus_musculus","Mouse")) species <- "mm10"
if(unique(df_samples$species) %in%  c("Danio_rerio")) species <- "danRer10"
#if(species == "hg19") suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
#if(species == "mm10") suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))

message("Copying the datasets")
if(length(missing_data)>0){
        # Move files from Download to ./data and rename them
        for(missing_dat in missing_data){
                old.pth  <- paste("~/Downloads", missing_dat,"outs",
                                  "filtered_feature_bc_matrix",sep = "/")
                list.of.files <- list.files(old.pth)
                new.folder <- paste("./data/scRNA-seq", missing_dat,"outs",
                                    "filtered_gene_bc_matrices",species,sep = "/")
                if(!dir.exists(new.folder)) dir.create(new.folder, recursive = T)
                # copy the files to the new folder
                file.copy(paste(old.pth, list.of.files, sep = "/"), new.folder)
                print(list_files <- list.files(new.folder))
        }
}

message("Loading the datasets")
## Load the dataset
Seurat_raw <- list()
Seurat_list <- list()
for(i in seq_along(df_samples$sample)){
        Seurat_raw[[i]] <- Read10X(data.dir = paste0("data/",df_samples$sample.id[i],
                                   "/outs/filtered_feature_bc_matrix"))
        colnames(Seurat_raw[[i]]) = paste0(df_samples$sample[i],"_",colnames(Seurat_raw[[i]]))
        Seurat_list[[i]] <- CreateSeuratObject(Seurat_raw[[i]],
                                                 min.cells = 0,
                                               min.features = 0)
        Seurat_list[[i]]@meta.data$tests <- df_samples$tests[i]
        Progress(i, length(df_samples$sample))
}
remove(Seurat_raw);GC()

#========1.1.3 g1 QC plots before filteration=================================
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), Seurat_list)
remove(Seurat_list);GC()

# read and select mitochondial genes
if(species == "hg19") mito = "^MT-"
if(species == "mm10") mito = "^mt-" # not Mt-
if(species == "danRer10") mito = "^mt-"
message("mito.genes:")

(mito.features <- grep(pattern = mito, x = rownames(object), value = TRUE))
object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = mito)
Idents(object) = factor(Idents(object),levels = df_samples$sample)
g1 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
        VlnPlot(object = object, features = features, ncol = 3, pt.size = 0.01)+
                theme(axis.text.x = element_text(size=10),legend.position="none")
})
save(g1,file= paste0(path,"g1","_",length(df_samples$sample),"_",gsub("-","",Sys.Date()),".Rda"))

#============1.2 scatter ======================
Seurat_list <- SplitObject(object, split.by = "orig.ident")
remove(object);GC()

for(i in 1:length(df_samples$sample)){
        high.mito <- isOutlier(Seurat_list[[i]]$percent.mt, nmads=3, type="higher")
        low.lib <- isOutlier(log10(Seurat_list[[i]]$nCount_RNA), type="lower", nmad=3)
        low.genes <- isOutlier(log10(Seurat_list[[i]]$nFeature_RNA), type="lower", nmad=3)
        discard <- high.mito | low.lib | low.genes
        print(data.frame(HighMito= sum(high.mito),LowLib=sum(low.lib), 
                         LowNgenes=sum(low.genes),Discard=sum(discard)))
        Seurat_list[[i]] <- Seurat_list[[i]][,!discard]
        #print(summary(!discard))
        print(i)
}

object <- Reduce(function(x, y) merge(x, y, do.normalize = F), Seurat_list)
remove(Seurat_list);GC()

object %<>% subset(subset = nFeature_RNA > 400 & nCount_RNA > 700 & percent.mt < 10)
# FilterCellsgenerate Vlnplot before and after filteration
Idents(object) = factor(Idents(object),levels = df_samples$sample)

g2 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
        VlnPlot(object = object, features = features, ncol = 3, pt.size = 0.01)+
                theme(axis.text.x = element_text(size=10),legend.position="none")
})
save(g2,file= paste0(path,"g2","_",length(df_samples$sample),"_",gsub("-","",Sys.Date()),".Rda"))
jpeg(paste0(path,"S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nFeature_RNA before filteration")+
                        scale_y_log10(limits = c(100,10000))+
                        theme(axis.text.x = element_text(size=13),
                              plot.title = element_text(hjust = 0.5)),
                g2[[1]]+ggtitle("nFeature_RNA after filteration")+
                        scale_y_log10(limits = c(100,10000))+
                        theme(axis.text.x = element_text(size=13),
                              plot.title = element_text(hjust = 0.5))))
dev.off()
jpeg(paste0(path,"S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nCount_RNA before filteration")+
                        scale_y_log10(limits = c(500,100000))+
                        theme(axis.text.x = element_text(size=13),
                              plot.title = element_text(hjust = 0.5)),
                g2[[2]]+ggtitle("nCount_RNA after filteration")+ 
                        scale_y_log10(limits = c(500,100000))+
                        theme(axis.text.x = element_text(size=13),
                              plot.title = element_text(hjust = 0.5))))
dev.off()
jpeg(paste0(path,"S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % before filteration")+
                        ylim(c(0,50))+
                        theme(axis.text.x = element_text(size=13),
                              plot.title = element_text(hjust = 0.5)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                        ylim(c(0,50))+
                        theme(axis.text.x = element_text(size=13),
                              plot.title = element_text(hjust = 0.5))))
dev.off()

#====
Seurat_list <- SplitObject(object, split.by = "orig.ident")
remove(object);GC()
#======1.1.2 record data quality before removing low quanlity cells =========================
# if args 2 is passed
message("QC")
cell.number <- sapply(Seurat_list, function(x) length(colnames(x)))
nCount_RNA <- sapply(Seurat_list, function(x) mean(x$nCount_RNA))
nFeature_RNA <- sapply(Seurat_list, function(x) mean(x$nFeature_RNA))
percent.mt <- sapply(Seurat_list, function(x) mean(x$percent.mt))
QC.list <- cbind(df_samples,cell.number, nCount_RNA, nFeature_RNA,percent.mt,
                 row.names = df_samples$sample)
write.csv(QC.list,paste0(path,"QC_list_",gsub("-","",Sys.Date()),".csv"))
#QC.list %>% kable() %>% kable_styling()
remove(cell.number,nCount_RNA,nFeature_RNA,percent.mt,QC.list);GC()

#========1.4 scran ===============================

sce_list <- lapply(Seurat_list, as.SingleCellExperiment)
remove(Seurat_list);GC()

# cluster
set.seed(1000)
clusters_list <- lapply(sce_list,function(x){
        scran::quickCluster(x, use.ranks=FALSE, BSPARAM=BiocSingular::IrlbaParam())
}) 
sapply(clusters_list, table)

# computeSumFactors and normalize
sce_list <- mapply(function(x,y){
        computeSumFactors(x, min.mean=0.1, cluster=y)},
        x=sce_list,y=clusters_list)
lapply(sce_list,function(x) summary(sizeFactors(x)))
remove(clusters_list);GC()
par(mfrow=c(1,1))
plot(sce_list[[1]]$nCount_RNA, sizeFactors(sce_list[[1]]), log="xy")
sce_list <- lapply(sce_list, scater::normalize)

save(sce_list, file = paste0("data/","sce_",length(df_samples$sample),"_",gsub("-","",Sys.Date()),".Rda"))
