########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(sapply(c("Seurat","tidyr","magrittr","dplyr","gplots"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 2.1 identify phenotype for each cluster  ==========================================
args <- commandArgs(trailingOnly = TRUE)
args[1] = as.character(args[1])
(load(file = args[1]))
df_markers <- readxl::read_excel("../seurat_resources/bio-rad-markers.xlsx")
#df_markers <- readxl::read_excel("../seurat_resources/bio-rad-markers.xlsx",sheet = "Human.sub")
colnames(df_markers) = gsub(" ","_",colnames(df_markers))
colnames(df_markers) = gsub(":|\\/","_",colnames(df_markers))
colnames(df_markers) = gsub("\\+","",colnames(df_markers))
markers = df_markers[,-grep("Alias",colnames(df_markers))]
marker.list <- df2list(markers)

marker.list %<>% lapply(function(x) x[1:16]) %>% 
        lapply(function(x) FilterGenes(object,x)) %>% 
        lapply(function(x) x[!is.na(x)])
#marker.list %>% list2df %>% t %>% kable() %>% kable_styling()

for(i in 1:length(marker.list)){
        p <- lapply(marker.list[[i]], function(marker) {
                SingleFeaturePlot.1(object = object, feature = marker,pt.size = 0.5,
                                    gradient.use = c("lightblue", "blue3"),threshold=0.1)+
                        ggtitle(paste0(marker,Alias(df = df_markers,gene = marker)))+
                        theme(plot.title = element_text(hjust = 0.5,size = 15, face = "bold"))
                })

        jpeg(paste0(path,names(marker.list)[i],".jpeg"),
             units="in", width=10, height=7,res=600)
        print(do.call(plot_grid, p)+ ggtitle(paste(names(marker.list)[i],"markers"))+
                      theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold")))
        dev.off()
        print(paste0(i,":",length(marker.list)))
}