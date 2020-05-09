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
library(ggpubr)
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
(load(file = "data/Lung_5_CCA_20200313.Rda"))
DefaultAssay(object) = "SCT"
Idents(object) = "integrated_snn_res.0.6"
T_CD4.markers <- FindMarkers.UMI(object = object, ident.1 = 5, only.pos = F, 
                                 logfc.threshold = 0.1, test.use = "MAST")
write.csv(T_CD4.markers, paste0(path, "T_CD4-vs-everything_else.csv"))
T_CD4_Treg.markers <- FindMarkers.UMI(object = object, ident.1 = 5, ident.2 = 14, 
                                 only.pos = F, 
                                 logfc.threshold = 0.1, test.use = "MAST")
write.csv(T_CD4_Treg.markers, paste0(path, "T_CD4-vs-Treg.csv"))

Idents(object) = "orig.ident"
exp <- AverageExpression(object, assays = "SCT")
write.csv(exp$SCT, paste0(path, "Average_Expression_by_sample.csv"))

#
Signatures1 <- readxl::read_excel("doc/Signatures.xlsx", sheet = "Primary_Signatures")
Signatures2 <- readxl::read_excel("doc/Signatures.xlsx", sheet = "Additional_Signatures")
Signatures1 %<>% df2list
Signatures2 %<>% df2list
Signatures <- c(Signatures1, Signatures2)
Signatures %<>% lapply(function(sig) FilterGenes(object, marker.genes = sig))
object %<>% AddModuleScore(features = Signatures, assay = "SCT", name = names(Signatures))
Colnames =  colnames(object@meta.data)
ind = grep(paste(names(Signatures), collapse = "|"),Colnames)
colnames(object@meta.data)[ind] = names(Signatures)

object$orig.ident %<>% factor(levels = c("IRE-CD45-plus-D7","IRE-CD45-plus-D9",
                                         "SCR-CD45-plus-D7", "SCR-CD45-plus-D9",
                                         "Naive-CD45-plus-D7"))

for(i in seq_along(Signatures)){
        jpeg(paste0(path,"FeaturePlot_",names(Signatures)[i],".jpeg"),
             units="in", width=8, height=7,res=600)
        print(FeaturePlot.1(object, features = names(Signatures)[i], reduction = "tsne",
                            split.by = "orig.ident", ncol = 2,
                            title = names(Signatures)[i]))
        dev.off()
        Progress(i, length(Signatures))
}
features <- c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4")
features <- FilterGenes(object, features)
DotPlot(object, features = features, cols = c("blue","red"),
        split.by = "orig.ident") + RotatedAxis()

# subset signature feature plot
jpeg(paste0(path,"TSNEPlot_subset_label.jpeg"), units="in", width=10, height=10,res=600)
TSNEPlot.1(object, group.by="singler1sub",
           cols = ExtractMetaColor(object),
           label = T,pt.size = 1,no.legend = T,label.repel = T,
           label.size = 4, repel = T,do.return= T,do.print = F,alpha = 0.9,
           title = "Compare T cell types in each sample") +
        #rectangle(-32, -24, -25, -16,colour = "black") + # alveolar_macrophages
        #rectangle(-26, -19, -35, -29,colour = "black") + # pDc
        #rectangle(-20, -10, 37, 43,colour = "black") + # cDc2
        #rectangle(25, 35, -22, -12,colour = "black") + # cDC1
        #rectangle(2.5, 14, 25, 30,colour = "black") + # Cd8_Effector
        #rectangle(14, 20.5, 28, 38,colour = "black")+ #Cd8_Exhausted
        #rectangle(21, 27, 9, 17,colour = "black")+ #Th
        #rectangle(27, 36, 18.5, 29,colour = "black")+ #Treg
        rectangle(13, 24, -32, -26,colour = "black") #MDSC
dev.off()

cell.embeddings = object@reductions$tsne@cell.embeddings
# alveolar_macrophages
alveolar_macrophages = cell.embeddings[,"tSNE_1"] > -33 & cell.embeddings[,"tSNE_1"] < -24 & 
                        cell.embeddings[,"tSNE_2"] > -25 & cell.embeddings[,"tSNE_2"] < -19
alveolar_macrophages = names(alveolar_macrophages)[alveolar_macrophages]
# pDc
pDc = cell.embeddings[,"tSNE_1"] > -26 & cell.embeddings[,"tSNE_1"] < -15 & 
        cell.embeddings[,"tSNE_2"] > -35 & cell.embeddings[,"tSNE_2"] < -29
pDc = names(pDc)[pDc]
# cDc2
cDc2 = cell.embeddings[,"tSNE_1"] > -20 & cell.embeddings[,"tSNE_1"] < -10 & 
        cell.embeddings[,"tSNE_2"] > 37 & cell.embeddings[,"tSNE_2"] < 43
cDc2 = names(cDc2)[cDc2]
# cDC1
cDC1 = cell.embeddings[,"tSNE_1"] > 25 & cell.embeddings[,"tSNE_1"] < 35 & 
        cell.embeddings[,"tSNE_2"] > -22 & cell.embeddings[,"tSNE_2"] < -12
cDC1 = names(cDC1)[cDC1]
# Cd8_Effector
Cd8_Effector = cell.embeddings[,"tSNE_1"] > 2.5 & cell.embeddings[,"tSNE_1"] < 14.5 & 
        cell.embeddings[,"tSNE_2"] > 25 & cell.embeddings[,"tSNE_2"] < 30
Cd8_Effector = names(Cd8_Effector)[Cd8_Effector]
# Cd8_Exhausted
Cd8_Exhausted = cell.embeddings[,"tSNE_1"] > 14 & cell.embeddings[,"tSNE_1"] < 20.5 & 
        cell.embeddings[,"tSNE_2"] > 28 & cell.embeddings[,"tSNE_2"] < 38
Cd8_Exhausted = names(Cd8_Exhausted)[Cd8_Exhausted]
deparse(substitute(Cd8_Exhausted))

sub_cells <- list(alveolar_macrophages, pDc, cDc2, cDC1, Cd8_Effector, Cd8_Exhausted)
sub_cells_names <- c("alveolar macrophages", "pDc", "cDc2", "cDC1", "Cd8 Effector", "Cd8 Exhausted")

# subset tSNE plots and Signatures FeaturePlot
for(i in seq_along(sub_cells)){
        sub_object <- subset(object, cells = sub_cells[[i]])
        TSNEPlot.1(sub_object, group.by="singler1sub",
                   cols = ExtractMetaColor(sub_object),
                   label = T,pt.size = 1,no.legend = T,label.repel = T,
                   save.path = paste0(path,sub_cells_names[i],"/"),
                   label.size = 4, repel = T,do.return= F,do.print = T,alpha = 0.9,
                   title = paste("tSNE Plot of", sub_cells_names[i]))
        for(m in seq_along(Signatures)){
                jpeg(paste0(path,sub_cells_names[i],"/FeaturePlot_",sub_cells_names[i],
                            "_",names(Signatures)[m],".jpeg"), units="in", width=8, height=7,res=600)
                print(FeaturePlot.1(sub_object, features = names(Signatures)[m], reduction = "tsne",
                                    split.by = "orig.ident", ncol = 2,
                                   cols = c("lightgrey", "blue"),
                                    title = names(Signatures)[m]))
                dev.off()
                Progress(m, length(Signatures))
        }
        df = sub_object@meta.data[,names(Signatures)]
        write.csv(df, paste0(path,sub_cells_names[i],"/signatures_expression.csv"))
        Progress(i, length(sub_cells))
}

# subset Signatures boxplot
for(i in seq_along(sub_cells)){
        sub_object <- subset(object, cells = sub_cells[[i]])
        for(m in seq_along(Signatures)){
                jpeg(paste0(path,sub_cells_names[i],"/FeaturePlot_",sub_cells_names[i],
                            "_",names(Signatures)[m],".jpeg"), units="in", width=8, height=7,res=600)
                print(FeaturePlot.1(sub_object, features = names(Signatures)[m], reduction = "tsne",
                                    split.by = "orig.ident", ncol = 2,
                                    cols = c("lightgrey", "blue"),
                                    title = names(Signatures)[m]))
                dev.off()
                Progress(m, length(Signatures))
        }
        df = sub_object@meta.data[,names(Signatures)]
        write.csv(df, paste0(path,sub_cells_names[i],"/signatures_expression.csv"))
        Progress(i, length(sub_cells))
}


features <- names(Signatures)       
g <- list()
sub_cells_names <- c("alveolar macrophages", "pDc", "cDc2", "cDC1", "Cd8 Effector", "Cd8 Exhausted")
for(i in seq_along(sub_cells_names)){
        sub_object <- subset(object, cells = sub_cells[[i]])
        df = sub_object@meta.data[,c("orig.ident",features)]
        colnames(df) %<>% gsub(".*-","",.)
        df$gene = rownames(df)
        df1 <- reshape2::melt(df)
        colnames(df1) = c("gene","sample", "UMI")
        g <- ggbarplot(df1, x = "variable", y = "value", color = "orig.ident",
                            add = "mean_se", #palette = c("#00AFBB", "#E7B800"),
                            legend = "right",
                            position = position_dodge())+
                ggtitle(sub_cells_names[i])+
                theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))
}
jpeg(paste0(path,"ggbarplot.jpeg"), units="in", width=10, height=7,res=600)
print(cowplot::plot_grid(g[[1]]+theme(axis.title.x=element_blank()),
                         g[[2]]+theme(axis.title.x=element_blank()),
                         g[[3]]+theme(axis.title.x=element_blank()),
                         ncol = 1,align = "hv"))
dev.off()
