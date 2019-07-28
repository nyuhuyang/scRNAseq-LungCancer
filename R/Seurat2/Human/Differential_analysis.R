########################################################################
#
#  0 setup environment, install libraries if nLynchessary, load libraries
# 
# ######################################################################

invisible(sapply(c("Seurat","magrittr","tidyr","dplyr","kableExtra","gplots"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("./data/")) dir.create("data")

#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
args <- commandArgs(trailingOnly = TRUE)
#(load(file = "./data/MouseTumor_2_20190219.Rda"))
args[1] = as.character(args[1])

# Manually RenameIdent cell type markers ================
object %<>% SetAllIdent(id = "res.0.6")
new.ident <- c("Monocytes/Macrophages",
               "Monocytes/Macrophages",
               "Monocytes/Macrophages",
               "Monocytes/Macrophages",
               "Monocytes/Macrophages", 
               "Monocytes/Macrophages",
               "Monocytes/Macrophages",
               "Monocytes/Macrophages",
               "T/NK cells",
               "T/NK cells",
               "Fibroblasts",
               "Monocytes/Macrophages",
               "T/NK cells",
               "Monocytes/Macrophages",
               "B cells")
for (i in 0:14) {
        object <- RenameIdent(object = object, old.ident.name = i, 
                              new.ident.name = new.ident[i + 1])
}

p4 <- TSNEPlot.1(object = object, do.label = T, group.by = "ident",
                 do.return = TRUE, no.legend = T,
                 #colors.use = ExtractMetaColor(object),
                 pt.size = 1,label.size = 4,force = 2)+
        ggtitle("Major cell types")+
        theme(text = element_text(size=10),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"PlotTsne_manual.jpeg"), units="in", width=10, height=7,res=600)
print(p4)
dev.off()
object %<>% StashIdent(save.name = "manual")

# Identify cell type markers ================
gde.all <- FindAllMarkers.UMI(object,logfc.threshold = 0.25, only.pos = F,return.thresh = 0.05, 
                                 test.use = "MAST")
write.csv(gde.all, paste0(path,"all_markers.csv"))
head(gde.all,10) %>% kable %>% kable_styling
# FindPairMarkers
(ident <- unique(object@meta.data$manual) %>% sort)
print(ident.1 <- paste0(ident, "_NAM"))
print(ident.2 <- paste0(ident, "_Control"))
object@meta.data$manual = paste0(object@meta.data$manual, "_", object@meta.data$orig.ident)
object %<>% SetAllIdent(id = "manual")
table(object@ident)
gde.pair <- FindPairMarkers(object, ident.1 = dent.1, ident.2 = ident.2,
                               logfc.threshold = 0.25, min.cells.group =3,
                               return.thresh = 0.05, only.pos = FALSE)
write.csv(gde.pair, paste0(path,"pairwise_comparision.csv"))
head(gde.pair,10) %>% kable %>% kable_styling
