########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
library(ggplot2)
library(fgsea)
library(tibble)
library(ggpubr)
library(ggsci)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# 3.1.1 load data
# Load some results from Seurat
#============ read pathways  =====================
(load(file="data/Lung_6_20190817.Rda"))

res = read.csv(file="output/20190817/DEG/IRE1a WT_IRE1aKO_Naive_SCT.csv",
                        row.names = 1, stringsAsFactors=F)
table(res$cluster)
head(res)
res = res[order(res["p_val_adj"]),]
head(res, 20)
(clusters <- unique(res$cluster))
hallmark <- gmtPathways("../seurat_resources/msigdb/h.all.v6.2.symbols.gmt")
biocarta <- gmtPathways("../seurat_resources/msigdb/c2.cp.biocarta.v6.2.symbols.gmt")
kegg <- gmtPathways("../seurat_resources/msigdb/c2.cp.kegg.v6.2.symbols.gmt")
tft <- gmtPathways("../seurat_resources/msigdb/c3.tft.v6.2.symbols.gmt")
c6 <- gmtPathways("../seurat_resources/msigdb/c6.all.v6.2.symbols.gmt")
GO <- gmtPathways("../seurat_resources/msigdb/c5.all.v6.2.symbols.gmt")

hallmark_mouse = readRDS(file ="../seurat_resources/msigdb/mouse_h.all.v6.2.symbols.rds")
allpathways_mouse = readRDS(file ="../seurat_resources/msigdb/mouse_allpathways.v6.2.symbols.rds")

hallmark %>% head() %>% lapply(head)
hallmark_mouse <- vector("list",length(hallmark))
names(hallmark_mouse) = gsub("HALLMARK_","",names(hallmark))

for (i in 1:length(hallmark)){
        hallmark_mouse[[i]] <- FilterGenes(object,hallmark_mouse[[i]],verbose = F)
        #print(paste0(i,":",length(hallmark)))
}
lapply(hallmark_mouse,length)
list.save(hallmark, "../seurat_resources/msigdb/mouse_h.all.v6.2.symbols.rds")
allpathways_mouse <- vector("list",length(allpathways))
names(allpathways_mouse) = names(allpathways)
for (i in 1:length(allpathways)){
        allpathways_mouse[[i]] <- Human2Mouse(allpathways[[i]])
        print(paste0(i,":",length(allpathways)))
}
list.save(allpathways_mouse, "../seurat_resources/msigdb/mouse_allpathways.v6.2.symbols.rds")

# Now, run the fgsea algorithm with 1000 permutations:
(titles <- paste(ident.2,"vs.", ident.1))
(clusters <- unique(gde.pair$cluster1.vs.cluster2))

for(i in 1:length(clusters)) FgseaBarplot(pathways=hallmark_mouse, stats=res, nperm=1000,
                               cluster = clusters[i],no.legend = F,
                               cut.off = "padj",cut.off.value = 0.25,
                               sample="Lynch Syndrome",pathway.name = "Hallmark", hjust=0.5,
                               width=10, height=7)
res$cluster1.vs.cluster2 %<>% plyr::mapvalues(from = clusters,
                                              to = celltypes)

for(i in 1:length(clusters)) FgseaBarplot(pathways=allpathways_mouse, stats=res, nperm=1000,
                                          cluster = clusters[i],no.legend = F,
                                          cut.off = "padj",cut.off.value = 0.25,
                                          sample="Lynch Syndrome",pathway.name = "All_pathways", hjust=0.5,
                                          width=10, height=7)

FgseaDotPlot(stats=res, pathways=hallmark_mouse, nperm=1000,pval = 0.01,
             order.by = c(4,"NES"),decreasing = F,
             size = "-log10(pval)", fill = "NES",sample = "each cell type", 
             pathway.name = "Hallmark",rotate.x.text = T)

fgseaRes <- FgseaDotPlot(stats=res, pathways=allpathways_mouse, nperm=1000,padj = 0.1,pval = 0.005,
             order.by = c(4,"NES"),decreasing = F,do.return = T,
             size = "-log10(pval)", fill = "NES",sample = "each B_MCL clusters", 
             rotate.x.text = F, pathway.name = "Hallmark, biocarta,and KEGG")
df_fgseaRes <- data.table::rbindlist(fgseaRes) %>% as.data.frame()
write.csv(df_fgseaRes, paste0(path,"GSEA_all_pathways.csv"))



df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
groups = c("Untreated","Pt-17","Pt-25")
for(i in 1:length(groups)){
        res_B = read.csv(file = paste0("output/20190622/B/B_MCL_DE/B_",groups[i],".csv"))
        res_T = read.csv(file = paste0("output/20190622/T/T_NK_DE/T_",groups[i],".csv"))
        
        (samples = df_samples$sample[df_samples$sample %in% unique(res_B$cluster)])
        #res_B$cluster %<>% factor(levels = samples)
        #res_T$cluster %<>% factor(levels = samples)
        
        FgseaDotPlot(stats=res_B, pathways=hallmark, nperm=1000,padj = 0.25,pval = 0.05,
                     order.by = c(4,"NES"),decreasing = F,
                     size = "-log10(pval)", fill = "NES",
                     sample = paste(groups[i],"B_MCL clusters"), 
                     pathway.name = "Hallmark",rotate.x.text = F)
        (samples = df_samples$sample[df_samples$sample %in% unique(res_T$cluster)])
        FgseaDotPlot(stats=res_T, pathways=hallmark, nperm=1000,padj = 0.25,pval = 0.05,
                     order.by = c(4,"NES"),decreasing = F,
                     size = "-log10(pval)", fill = "NES",
                     sample = paste(groups[i],"T_NK clusters"), 
                     pathway.name = "Hallmark",rotate.x.text = F)

}

#============ T cells  =====================
df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
tests <- paste0("test",4)
sample_n = which(df_samples$tests %in% tests)
df <- as.data.frame(df_samples[sample_n,])
(samples <- c("Normal",unique(df$sample)))
cell.type <- c("T_cells:CD4+","T_cells:CD8+","NK_cells")

res_T_list <- list()
for(i in 2:length(samples)){
        res_T_list[[i-1]] = read.csv(file = paste0(path,samples[i],"_vs_Normal",".csv"))
}
res_T <- do.call("rbind.data.frame", res_T_list)
res_T = res_T[grep("CD8+",res_T$cluster1.vs.cluster2),]
res_T$cluster1.vs.cluster2 %<>% as.character %>% gsub('\\..*',"",.)

res_T$cluster1.vs.cluster2 %<>% factor(levels = samples[2:5])
        
for(i in 2:length(samples)) FgseaBarplot(pathways=hallmark, stats=res_T, nperm=1000,
                                          cluster = samples[i],no.legend = F,
                                          cut.off = "padj",cut.off.value = 0.25,
                                          sample="CD8+ T cells of",pathway.name = "Hallmark", hjust=0.5,
                                          width=10, height=7)

FgseaDotPlot(stats=res_T, pathways=hallmark, nperm=1000,padj = 0.25,pval = 0.05,
                     order.by = c(4,"NES"),decreasing = F,
                     size = "-log10(pval)", fill = "NES",
                     sample = paste("Pt-17's CD8+ T cells"), 
                     pathway.name = "Hallmark",rotate.x.text = F)
