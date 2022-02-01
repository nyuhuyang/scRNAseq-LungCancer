# conda activate r4.1.1
invisible(lapply(c("Seurat","dplyr","magrittr","tidyr",
                   "future"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# Need 64GB ?
set.seed(101)
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

object = readRDS(file = "data/Lung_5_20200313.rds")

step = c("cell_types")[1]

if(step == "cell_types"){# 32~64GB

    DefaultAssay(object) = "SCT"
    object %<>% subset(subset = celltype.l2 != "unknown")
    
    opts = data.frame(ident = c(rep("singler1sub",20),
                                rep("singler1main",12),
                                rep("celltype.l3",42),
                                rep("celltype.l2",22),
                                rep("celltype.l1",8)),
                      num = c(1:20,
                              1:12,
                              1:42,
                              1:22,
                              1:8)
                      )
    opt = opts[args,] # 1-104

    #==========================
    Idents(object) = opt$ident
    opt$type = sort(levels(object))[opt$num]
    print(opt)
    

    markers = FindMarkers_UMI(object, ident.1 = opt$type,
                              group.by = opt$ident,
                              assay = "SCT",
                              logfc.threshold = 0.1,
                             only.pos = F,
                             test.use = "wilcox")
    markers$cluster = as.character(opt$type)
    markers$Cell_category = opt$ident
    num = opt$num
    if(args < 10) num = paste0("0",num)
    if(args < 100) num = paste0("0",num)

    arg = args
    if(args < 10) arg = paste0("0",arg)
    if(args < 100) arg = paste0("0",arg)

    save_path <- paste0(path,step,"/")
    if(!dir.exists(save_path)) dir.create(save_path, recursive = T)
    
    write.csv(markers,paste0(save_path,arg,"-",opt$ident,"-",num,".",opt$type, ".csv"))
}