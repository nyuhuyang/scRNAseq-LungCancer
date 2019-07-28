library(SingleR)
source("../R/SingleR_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/Lung.data_harmony_6_2019728.Rda"))
attach(mouse.rnaseq)
singler <- CreateBigSingleRObject(object_data, annot = NULL, project.name="Mouse_LungCancer",
                                    N = 5000, min.genes = 200, technology = "10X",
                                    species = "Mouse", citation = "", ref.list = list(mouse.rnaseq),
                                    normalize.gene.length = F, variable.genes = "de", 
                                    fine.tune = T, do.signatures = F, clusters = NULL)
save(singler,file="output/singler_T_Lung_6_2019728.Rda")
  
