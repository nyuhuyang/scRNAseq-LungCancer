########################################################################
#
#  setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(sapply(c("Seurat","magrittr","dplyr","reshape2",
                   "pheatmap","tidyr","kableExtra"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
args <- commandArgs(trailingOnly = TRUE)
args = sapply(args,as.character)
#(load(file="data/LungCancer_4_20190305.Rda"))
(load(file =args[1]))
#args="doc/190212_scRNAseq_info.xlsx"
df_samples <- readxl::read_excel(args[2])
print(df_samples)
df_samples = as.data.frame(df_samples)
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",1:10)))

##############################
# subset Seurat
###############################
table(object@meta.data$orig.ident)
object %<>% SetAllIdent(id = "orig.ident")
table(object@ident)

tests <- paste0("test",1:2)
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        df <- as.data.frame(df_samples[sample_n,])
        samples <- unique(df$sample)
        rownames(df) = samples
        
        #print(samples <- df$sample[order(df$tsne)])
        
        g <- lapply(samples,function(sample) {
                SubsetData(object, ident.use = sample) %>%
                        SetAllIdent(id = "singler2main") %>%
                        TSNEPlot.1(no.legend = T,do.label =T,label.size=3,size=20,
                                   colors.use = ExtractMetaColor(.),
                                   return.plots =T, label.repel = T,force=2)+
                        ggtitle(sample)+theme(text = element_text(size=20),
                                              plot.title = element_text(hjust = 0.5))
        })
        jpeg(paste0(path,test,"_main_TSNEPlot.jpeg"), units="in", width=10, height=7,
             res=600)
        print(do.call(plot_grid, c(g, ncol = 2)))
        dev.off()
}