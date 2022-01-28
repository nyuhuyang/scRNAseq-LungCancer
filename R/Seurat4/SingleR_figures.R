# conda activate r4.0.3
library(Seurat)
library(magrittr)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggpubr)
library(S4Vectors)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================

##############################
# create singleR data frame
###############################
pred = readRDS("output/Lung_5_20200313_singleR_azimuth_PBMC.rds")
(load("data/Lung_5_20200313.Rda"))
object = UpdateSeuratObject(object)

singlerDF = data.frame("celltype.l3" = pred$pruned.labels,
                       row.names = rownames(pred))
table(is.na(singlerDF$celltype.l3))
singlerDF$celltype.l3[is.na(singlerDF$celltype.l3)]= "unknown"

##############################
# adjust cell label
##############################
path = "../seurat_resources/azimuth/PBMC/"
meta.data = read.csv(paste0(path,"GSE164378/GSE164378_sc.meta.data_5P.csv"),row.names =1)
meta.data = meta.data[!duplicated(meta.data$celltype.l3),]

singlerDF$celltype.l2 = plyr::mapvalues(singlerDF$celltype.l3,
                                        from = meta.data$celltype.l3,
                                        to = meta.data$celltype.l2)
singlerDF$celltype.l1 = plyr::mapvalues(singlerDF$celltype.l3,
                                        from = meta.data$celltype.l3,
                                        to = meta.data$celltype.l1)
##############################
# process color scheme
##############################
table(colnames(object) == rownames(singlerDF))
object@meta.data %<>% cbind(singlerDF)

saveRDS(object, file = "data/Lung_5_20200313.rds")


# by barplot
cell_Freq <- table(object$label.fine) %>% as.data.frame
cell_Freq$Percent <- prop.table(cell_Freq$Freq) %>% round(digits = 2) %>% scales::percent()
cols = ExtractMetaColor(object)
cell_Freq$cols = cols[cell_Freq$Var1]
cell_Freq = cell_Freq[order(cell_Freq$Var1),]

cell_Freq = cell_Freq[order(cell_Freq$Freq,decreasing = T),]
cell_Freq$Var1 %<>% factor(levels = as.character(cell_Freq$Var1))
colnames(cell_Freq)[1:2] = c("Cell_Type", "Cell_Number")

jpeg(paste0(path,"cell_type_numbers.jpeg"), units="in", width=6, height=6,res=600)
ggbarplot(cell_Freq, "Cell_Type", "Cell_Number",
          fill = "Cell_Type", color = "black",xlab = "",
          palette = cell_Freq$col,x.text.angle = 90,
          ylab = "Cell Number",
          label = cell_Freq$Percent,
          lab.size = 3,
          sort.val = "desc",
          width = 1, size = 0.5,
          title = "Numbers of cell types in total 30 samples")+NoLegend()+
    theme(plot.title = element_text(hjust = 0.5,size=15),
          axis.text.x = element_text(vjust = 0.5))
dev.off()