########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
####################################
invisible(lapply(c("Seurat","dplyr","magrittr","tidyr","openxlsx",
                   "MAST","future","gplots"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
arg <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",arg))

opts = data.frame(c(-32, -24, -25, -16, "alveolar_macrophages"),
                  c(-26, -19, -35, -29, "pDc"),
                  c(-20, -10, 37, 43, "cDc2"),
                  c(25, 35, -22, -12, "cDC1"),
                  c(2.5, 14, 25, 30, "Cd8_Effector"),
                  c(14, 20.5, 28, 38, "Cd8_Exhausted"),
                  c(21, 27, 9, 17, "Th"),
                  c(27, 36, 18.5, 29, "Treg"),
                  c(13, 24, -32, -26, "MDSC"),
                  stringsAsFactors = F) %>% t %>% as.data.frame()

colnames(opts) = c("x_left", "x_right", "y_bottom","y_top", "cell.types")
opts[,1:4] %<>% apply(2, as.numeric)
(opt = opts[arg,])
save.path = paste0(path, opt$cell.types)

# load data
(load(file = "data/Lung_5_CCA_20200313.Rda"))
DefaultAssay(object) = "SCT"

# subset
cell.embeddings = object@reductions$tsne@cell.embeddings
keep_cells = cell.embeddings[,"tSNE_1"] > opt$x_left & cell.embeddings[,"tSNE_1"] < opt$x_right & 
        cell.embeddings[,"tSNE_2"] > opt$y_bottom & cell.embeddings[,"tSNE_2"] < opt$y_top

keep_cells = names(keep_cells)[keep_cells]
object %<>% subset(cells = keep_cells)
Idents(object) = "singler1sub"
TSNEPlot.1(object, group.by="singler1sub",
           cols = ExtractMetaColor(object),
           label = T,pt.size = 1,no.legend = T,label.repel = T,
           save.path = paste0(save.path,"-"),
           label.size = 4, repel = T,do.return= F,do.print = T,alpha = 0.9,
           title = paste("tSNE Plot of", opt$cell.types))

# Differential analysis
object$conditions %<>% gsub("-CD45","",.)
Idents(object) = "conditions"
table(Idents(object))
res_list <- list()
pairs = data.frame(c("IRE",  "SCR"),
                   c("Naive","IRE"),
                   c("Naive", "SCR"),
                   stringsAsFactors = F)
for(i in 1:3) {
        (pair = pairs[,i])
        sub_object <- subset(object, idents = c(pair[1], pair[2]))
        res_list[[i]] = FindAllMarkers.UMI(sub_object,
                                           logfc.threshold = 0.2, 
                                           only.pos = T,
                                           return.thresh = 1,
                                           test.use = "MAST",
                                           latent.vars = "nCount_SCT")
        names(res_list)[i] = paste0(pair[1],"_", pair[2])
}
write.xlsx(res_list, file = paste0(save.path,"-DE_results.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
