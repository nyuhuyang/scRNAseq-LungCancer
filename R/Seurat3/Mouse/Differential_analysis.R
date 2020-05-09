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
Idents(object) <-"singler1sub"
object <- sortIdent(object)
DefaultAssay(object) = "RNA"
#TSNEPlot(object)
Cancer.markers <- FindAllMarkers.UMI(object = object, only.pos = F, logfc.threshold = 0.5,
                                        test.use = "MAST")
write.csv(Cancer.markers,paste0(path,"Cancer_cell_type_markers.csv"))

#DoHeatmap.1======
Top_n = 10
top <-  Cancer.markers %>% group_by(cluster) %>% 
        top_n(Top_n, avg_logFC) %>% as.data.frame()
add.genes = unique(c(as.character(top$gene)))
object %<>% ScaleData(features= add.genes)
DoHeatmap.1(object, add.genes = add.genes,
            Top_n = Top_n, do.print=T, angle = 45,
            group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0.02, label=T, cex.row= 400/length(add.genes), 
            #legend.size = 0,
            width=10, height=6.5,unique.name = T, pal_gsea = T,
            title = paste("Top",Top_n,"DE genes in all cell types"))
GC()

# T cells ======================
meta.data = object@meta.data
meta.data = cbind(meta.data,object@reductions$umap@cell.embeddings)
T_cells <- meta.data$UMAP_1 < 0 & meta.data$UMAP_2 >-5
T_cells <- grepl("T_cells",object$singler1sub) & T_cells
UMAPPlot.1(object[,T_cells], group.by="singler1sub",
           cols = ExtractMetaColor(object[,T_cells]), ncol = 3,border = T,
           label = T,pt.size = 1,no.legend = T,label.repel = T,
           label.size = 4, repel = T,do.return= T,do.print = F,alpha = 0.9,
           title = "Compare T cell types in each sample")

ident.1 <- c("IRE1a WT", "IRE1aKO","IRE1aKO")
ident.2 <- c("Naive","Naive","IRE1a WT")

Idents(object) = "conditions"
table(Idents(object))
subfolder <- paste0(path,"DEG/")
DefaultAssay(object) = "SCT"
gde.pair <- FindPairMarkers(object[,T_cells], ident.1 = ident.1, ident.2 = ident.2,slot = "data",
                            logfc.threshold = 0.01, min.cells.group =3,assay.type = "SCT",
                            return.thresh = 1, only.pos = F, save.path = subfolder)
# Volcano plot=========
(titles <- paste(ident.2,"vs.", ident.1))
(clusters <- unique(gde.pair$cluster1.vs.cluster2))
for(i in 1:length(clusters)){
        df <- gde.pair[gde.pair$cluster1.vs.cluster2 %in% clusters[i],]
        df$log10_p_val_adj = -log10(df$p_val_adj)
        df$log10_p_val_adj[df$log10_p_val_adj == "Inf"] = 400
        left = df[df$avg_logFC < -0.1,]
        right = df[df$avg_logFC > 0.1,]
        left = rownames(left)[left$log10_p_val_adj >= head(sort(left$log10_p_val_adj,decreasing = T),15) %>% tail(1)]
        right = rownames(right)[right$log10_p_val_adj >= head(sort(right$log10_p_val_adj,decreasing = T),15) %>% tail(1)]
        g <- ggplot(df,aes(avg_logFC,log10_p_val_adj)) + 
                geom_point() + 
                ggtitle(titles[i]) + 
                ylab("-log10(p_value_adj)")+
                theme_minimal()+
                theme(plot.title = element_text(size=20,hjust = 0.5))+
                ggrepel::geom_text_repel(aes(label = gene), 
                                         data=df[c(left,right),]) +
                geom_point(color = ifelse((df$avg_logFC > 0.1  & df$p_val_adj < 0.05) , "red",
                                          ifelse((df$avg_logFC < -0.1 & df$p_val_adj < 0.05), "blue","gray")))
        jpeg(paste0(path,"Volcano_plot",clusters[i],".jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
}


# ggplot =========
Idents(object) = "orig.ident"

T_cells_exp <- AverageExpression(object[,T_cells],assays = "SCT")
T_cells_exp = T_cells_exp$SCT/log(2)
write.csv(T_cells_exp,paste0(path,"T_cells_exp.csv"))
T_cells_exp = T_cells_exp[unique(gde.pair$gene),]
T_cells_exp = T_cells_exp[,c("Naive", "IRE1aKO", "IRE1a WT")]
T_cells_delta <- sweep(T_cells_exp, 1, T_cells_exp$Naive,"-")
T_cells_delta = T_cells_delta[,-grep("Naive",colnames(T_cells_exp))]
colnames(T_cells_delta) = c("IRE1a_KO","IRE1a_WT")
T_cells_delta$genes = rownames(T_cells_delta)
dim(T_cells_delta);dim(T_cells_delta1)

data = T_cells_delta
ClevelandDotPlots <- function(data, x, y, top =20, group = NULL,color = "black",
                              sorting = c("ascending", "descending"), sort.by = "IRE1a_KO",
                              ggtheme = theme_bw(),palette = NULL,
                              y.text.col = FALSE, rotate = T, key = "strains",
                              value = "UMI_Increment", do.print = T,...){
        v <- deparse(substitute(data))
        v <- paste0(v, "_top", top, "by_",sort.by)

        not.sort.by = colnames(data)[!colnames(data) %in% c(sort.by,"genes")]
        data = data[data[,sort.by] > data[,not.sort.by[1]],]
        data %<>% gather(key = "strains", value = "UMI_Increment",-genes)
        genes <- data[,key] %in% sort.by %>% 
                data[., value] %>% 
                order(decreasing = ifelse(sorting =="descending",T,F)) %>%
                .[1:top] %>% data[.,"genes"]
        
        g <- ggdotchart(data[data$genes %in% genes,], x = x, y = y,
                   group = group, color = group,
                   palette = palette,
                   rotate = rotate,
                   sorting = sorting,
                   ggtheme = ggtheme,
                   y.text.col = y.text.col,...)+
                geom_line(aes(group = genes))+
                theme(text = element_text(size=12),
                      plot.title = element_text(size = 16,hjust = 0.5))
        
        if(do.print) {
                path <- paste0("output/",gsub("-","",Sys.Date()),"/")
                if(!dir.exists(path)) dir.create(path, recursive = T)
                jpeg(paste0(path,"ClevelandDotPlots_",v,".jpeg"),
                     units="in", width=10, height=7,res=600)
                print(g)
                dev.off()
        }
        return(g)
}


ClevelandDotPlots(T_cells_delta, x = "genes", y = "UMI_Increment",rotate = T,
                  top = 20, group = "strains", color = "strains",palette = c('#1F78B4','#E6AB02'),
                  sorting = "descending",sort.by = "IRE1a_WT",dot.size = 2,
                  title = "Top 20 upregulated genes in IRE1a_WT")

ClevelandDotPlots(T_cells_delta, x = "genes", y = "UMI_Increment",rotate = T,
                  top = 20, group = "strains", color = "strains",palette = c('#1F78B4','#E6AB02'),
                  sorting = "descending",sort.by = "IRE1a_KO",dot.size = 2,
                  title = "Top 20 upregulated genes in IRE1a_KO")

# T cell expression
T_cells <- object[,T_cells]
DefaultAssay(T_cells) = "SCT"
Idents(T_cells) = "orig.ident"
(samples <- as.character(unique(T_cells$orig.ident)))
subset_T_cells <- lapply(samples,function(x) subset(T_cells, idents = x))
for(i in 1:length(samples)){
        write.csv(subset_T_cells[[i]]@assays$SCT@data,
                  paste0(path,samples[i],"_T_cells.csv"))
}
