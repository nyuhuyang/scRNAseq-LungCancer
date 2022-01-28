########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
library(Seurat)
library(dplyr)
library(kableExtra)
library(magrittr)
library(ggplot2)
library(cowplot)
library(fgsea)
library(tibble)
library(ggsci)
library(progress)

source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

csv_list <- list.files(pattern = "FC0.01",path = "output/20210511",full.names = T)
deg_list <- pbapply::pblapply(csv_list, function(x){
    tmp = read.csv(x,row.names = 1)
    tmp = tmp[tmp$cluster == "KO", ]
    tmp$cluster %<>% paste0("/WT")
    tmp = tmp[order(tmp$avg_log2FC,decreasing = T), ]
    #tmp = tmp[tmp$avg_log2FC > 0, ]
    tmp

})


head(deg_list[[1]], 20)
# read pathway

hallmark <- fgsea::gmtPathways("../seurat_resources/msigdb/h.all.v7.4.symbols.gmt")
names(hallmark) = gsub("HALLMARK_","",names(hallmark))
names(hallmark) = gsub("\\_"," ",names(hallmark))


hallmark <- gmtPathways("../seurat_resources/msigdb/h.all.v7.4.symbols.gmt")
Biocarta <- gmtPathways("../seurat_resources/msigdb/c2.cp.biocarta.v7.4.symbols.gmt")
kegg <- gmtPathways("../seurat_resources/msigdb/c2.cp.kegg.v7.4.symbols.gmt")
pid <- gmtPathways("../seurat_resources/msigdb/c2.cp.pid.v7.4.symbols.gmt")
reactome <- gmtPathways("../seurat_resources/msigdb/c2.cp.reactome.v7.4.symbols.gmt")
tft <- gmtPathways("../seurat_resources/msigdb/c3.tft.v7.4.symbols.gmt")
go_cc <- gmtPathways("../seurat_resources/msigdb/c5.go.cc.v7.4.symbols.gmt")
go_bp <- gmtPathways("../seurat_resources/msigdb/c5.go.bp.v7.4.symbols.gmt")
go_mf <- gmtPathways("../seurat_resources/msigdb/c5.go.mf.v7.4.symbols.gmt")
go <- do.call(c, list(go_cc, go_bp, go_mf))

Development <- gmtPathways("../seurat_resources/msigdb/Development.gmt")
Physiology <- gmtPathways("../seurat_resources/msigdb/Physiology.gmt")
Disease <- gmtPathways("../seurat_resources/msigdb/Disease.gmt")
Cancer <- gmtPathways("../seurat_resources/msigdb/Cancer.gmt")

msigdb_list <- list("hallmark" = hallmark,
                    "Biocarta" = Biocarta,
                    "kegg" = kegg,
                    "reactome" = reactome,
                    "GO Biological Process" = go_bp,
                    "GO Cellular Component" = go_cc,
                    "GO Molecular Function" = go_mf)
msigdb_list %<>% pbapply::pblapply(function(x) {
    x %<>% lapply(function(x1) Hmisc::capitalize(tolower((x1))))
    x
})

hallmark %<>% lapply(function(x) Hmisc::capitalize(tolower((x))))
#=========
res =  bind_rows(deg_list)
res$cluster = res$cell.type
(clusters = unique(as.character(res$cluster)))

# hallmark
Fgsea_res <- FgseaDotPlot(stats=res, pathways=hallmark,Rowv = T,
                 title = "enriched hallmark pathways in KO",
                 plot.title = element_text(hjust = 1,size = 15),
                 axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
                 width = 6,do.return = T)
colnames(Fgsea_res)[5] = "cell.type"
openxlsx::write.xlsx(split(Fgsea_res,f = Fgsea_res$cell.type ),
                     file =  paste0(path,"hallmark_gsea.xlsx"),
                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))

# Biocarta
Fgsea_res <- FgseaDotPlot(stats=res, pathways=msigdb_list[["Biocarta"]],Rowv = T,
                          title = "enriched Biocarta pathways in KO",
                          plot.title = element_text(hjust = 1,size = 15),
                          axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
                          height= 10,width = 5,do.return = T)
colnames(Fgsea_res)[5] = "cell.type"
openxlsx::write.xlsx(split(Fgsea_res,f = Fgsea_res$cell.type ),
                     file =  paste0(path,"Biocarta_gsea.xlsx"),
                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))

# kegg
Fgsea_res <- FgseaDotPlot(stats=res, pathways=msigdb_list[["kegg"]],Rowv = T,
                          title = "enriched kegg pathways in KO",
                          plot.title = element_text(hjust = 1,size = 15),
                          axis.text.x = element_text(angle = 45, hjust = 1,size = 11),
                          height= 7.5,width = 6,do.return = T)
colnames(Fgsea_res)[5] = "cell.type"
openxlsx::write.xlsx(split(Fgsea_res,f = Fgsea_res$cell.type ),
                     file =  paste0(path,"kegg_gsea.xlsx"),
                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))


# reactome
Fgsea_res <- FgseaDotPlot(stats=res, pathways=msigdb_list[["reactome"]],
                          Rowv = T,
                          padj = 0.01, pval=0.001,
                          title = "enriched reactome pathways in KO",
                          plot.title = element_text(hjust = 1,size = 15),
                          axis.text.x = element_text(angle = 45, hjust = 1,size = 11),
                          height= 8,width = 12,do.return = T)
colnames(Fgsea_res)[5] = "cell.type"
openxlsx::write.xlsx(split(Fgsea_res,f = Fgsea_res$cell.type ),
                     file =  paste0(path,"reactome_gsea.xlsx"),
                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))

# go
Fgsea_res <- FgseaDotPlot(stats=res, pathways=msigdb_list[["GO Biological Process"]],
                          Rowv = T,
                          padj = 0.001, pval=0.0001,
                          title = "enriched GO Biological Process pathways in KO",
                          plot.title = element_text(hjust = 1,size = 15),
                          axis.text.x = element_text(angle = 45, hjust = 1,size = 11),
                          height= 10,width =8,do.return = T)
colnames(Fgsea_res)[5] = "cell.type"
openxlsx::write.xlsx(split(Fgsea_res,f = Fgsea_res$cell.type ),
                     file =  paste0(path,"go_bp_gsea.xlsx"),
                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))


Fgsea_res <- FgseaDotPlot(stats=res, pathways=msigdb_list[["GO Cellular Component"]],
                          Rowv = T,
                          padj = 0.05, pval=0.01,
                          title = "enriched GO Cellular Component pathways in KO",
                          plot.title = element_text(hjust = 1,size = 15),
                          axis.text.x = element_text(angle = 45, hjust = 1,size = 11),
                          height= 10,width =8,do.return = T)
colnames(Fgsea_res)[5] = "cell.type"
openxlsx::write.xlsx(split(Fgsea_res,f = Fgsea_res$cell.type ),
                     file =  paste0(path,"go_cc_gsea.xlsx"),
                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))

Fgsea_res <- FgseaDotPlot(stats=res, pathways=msigdb_list[["GO Molecular Function"]],
                          Rowv = T,
                          padj = 0.05, pval=0.01,
                          title = "enriched GO Molecular Function pathways in KO",
                          plot.title = element_text(hjust = 1,size = 15),
                          axis.text.x = element_text(angle = 45, hjust = 1,size = 11),
                          height= 8,width =9,do.return = T)
colnames(Fgsea_res)[5] = "cell.type"
openxlsx::write.xlsx(split(Fgsea_res,f = Fgsea_res$cell.type ),
                     file =  paste0(path,"go_mf_gsea.xlsx"),
                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))
