library(SingleR)
source("../R/SingleR_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/Lung.data_5_20200313.Rda"))
(load(file = "../SingleR/data/ref_GSE109125_mouse.rnaseq.RData"))
singler <- CreateBigSingleRObject.1(object_data, annot = NULL, project.name="Mouse_LungCancer",
                                    N = 5050, min.genes = 3, technology = "10X",
                                    species = "Mouse", citation = "", ref.list = list(ref),
                                    normalize.gene.length = F, variable.genes = "de", 
                                    fine.tune = T, do.signatures = F, do.main.types = T,
                                    temp.dir = getwd(), numCores = SingleR.numCores)
save(singler,file="output/singler_T_Lung_5_20200313.Rda")
