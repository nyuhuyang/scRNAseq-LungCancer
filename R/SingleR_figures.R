########################################################################
#
#  setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(sapply(c("Seurat","magrittr","SingleR","dplyr","reshape2",
                   "pheatmap","tidyr","kableExtra"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
args <- commandArgs(trailingOnly = TRUE)
args = sapply(args,as.character)

#(load(file="data/LungCancer_4_20190305.Rda"))
#(load(file="output/singler_T_LungCancer_4_20190305.Rda"))
(load(file =args[1]))
(load(file = args[2]))
 ##############################
# add singleR label to Seurat
###############################
singlerDF = data.frame("singler1sub" = FineTune(singler$singler[[1]]$SingleR.single$labels),
                       "singler1main" = FineTune(singler$singler[[1]]$SingleR.single.main$labels,
                                                 main.type = T),
                       "singler2sub" = FineTune(singler$singler[[2]]$SingleR.single$labels),
                       "singler2main" = FineTune(singler$singler[[2]]$SingleR.single.main$labels,
                                                 main.type = T),
                       "orig.ident" = object@meta.data$orig.ident,
                       row.names = rownames(singler$singler[[1]]$SingleR.single$labels))

table(rownames(singlerDF) %in% object@cell.names)

apply(singlerDF,2,function(x) length(unique(x)))
object <- AddMetaData(object = object,
                   metadata = singlerDF)
object <- SetAllIdent(object = object, id = "singler1sub")

##############################
# check the spearman correlation
###############################
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"/DrawHeatmap_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single, top.n = 50,normalize = F))
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

#kable(table(object@meta.data$singler2sub)) %>% kable_styling()
CellType <- table(object@meta.data$singler2sub, object@meta.data$orig.ident) %>% 
        as.data.frame() %>% spread(Var2,Freq)
colnames(CellType)[1] = "sub cell type"
write.csv(CellType,paste0(path,"SubCellType.csv"))

CellType <- table(object@meta.data$singler2main, object@meta.data$orig.ident) %>% 
        as.data.frame() %>% spread(Var2,Freq)
colnames(CellType)[1] = "cell type"
write.csv(CellType,paste0(path,"CellType.csv"))

##############################
# process color scheme
##############################

singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])

singler_colors1[duplicated(singler_colors1)]
singler_colors2[duplicated(singler_colors2)]
length(singler_colors1);length(singler_colors2)
apply(object@meta.data[,c("singler2sub","singler2main")],
      2,function(x) length(unique(x)))
#object@meta.data[,c("singler2sub")] %>% table() %>% kable() %>% kable_styling()
object <- AddMetaColor(object = object, label= "singler2sub", colors = singler_colors1)
object <- AddMetaColor(object = object, label= "singler2main", colors = singler_colors2)

#TSNEPlot.1(object, colors.use = ExtractMetaColor(object),no.legend = F)
##############################
# draw tsne plot
##############################
object <- SetAllIdent(object = object, id = "singler2main")
p3 <- TSNEPlot.1(object = object, do.label = T, group.by = "ident",
                 do.return = TRUE, no.legend = T,
                 colors.use = ExtractMetaColor(object),
                 pt.size = 1,label.size = 4,force = 2)+
  ggtitle("Supervised cell type labeling by Blueprint + Encode")+
  theme(text = element_text(size=10),							
        plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

object <- SetAllIdent(object = object, id = "singler2sub")
p4 <- TSNEPlot.1(object = object, do.label = T, group.by = "ident",
                 do.return = TRUE, no.legend = T,
                 colors.use = ExtractMetaColor(object),
                 pt.size = 1,label.size = 3,force = 2)+
        ggtitle("Supervised sub cell type labeling by Blueprint + Encode")+
        theme(text = element_text(size=10),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"PlotTsne_main.jpeg"), units="in", width=10, height=7,
     res=600)
print(p3)
dev.off()

jpeg(paste0(path,"PlotTsne_sub.jpeg"), units="in", width=10, height=7,
     res=600)
print(p4)
dev.off()

save(object,file=args[1])
