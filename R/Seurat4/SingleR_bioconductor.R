#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.0.3 linux
invisible(lapply(c("Seurat","SingleR","SingleCellExperiment",
                   "magrittr","data.table","Matrix"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

# Need 64GB ?
set.seed(101)


# ====== load single cell =============
(load("data/Lung_5_20200313.Rda"))
object = UpdateSeuratObject(object)
sce <- SingleCellExperiment(list(logcounts=object[["SCT"]]@data),
                                colData=DataFrame(object@meta.data))
rm(object);GC()

# ======= load azimuth PBMC data ==============================
path = "../seurat_resources/azimuth/PBMC/"
counts <- Read10X(paste0(path, "GSE164378/GSM5008740_RNA_5P"))
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
meta.data = read.csv(paste0(path,"GSE164378/GSE164378_sc.meta.data_5P.csv"),row.names =1)
table(rownames(meta.data) == colnames(counts))
rownames(counts) %<>% tolower() %>% Hmisc::capitalize()

PBMC <- SingleCellExperiment(list(logcounts=log1p(t(t(counts)/size.factors))),
                             colData=DataFrame(meta.data))
rm(counts,meta.data,libsizes,size.factors);GC()



# ====== conbime data =============

common <- Reduce(intersect, list(rownames(sce),
                                 rownames(PBMC)
))
length(common)
table(PBMC$celltype.l3)
system.time(trained <- trainSingleR(ref = PBMC[common,],
                                    labels=PBMC$celltype.l3))
system.time(pred <- classifySingleR(sce[common,], trained))
# elapsed 4872.846 sec
saveRDS(object = pred, file = paste0("output/Lung_5_20200313_singleR_azimuth_PBMC.rds"))
