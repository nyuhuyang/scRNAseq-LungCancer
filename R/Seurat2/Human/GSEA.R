########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
library(Seurat)
library(SingleR)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# 5.1 load data ==============
(load(file = "./data/MouseTumor_2_20190219.Rda"))
##################################
# select Monocytes only
##################################
object <- SetAllIdent(object, id="manual")
table(object@ident)
TSNEPlot.1(object,do.label = T)
Monocytes <- SubsetData(object, ident.use = "Monocytes/Macrophages")
table(Monocytes@meta.data$singler1main)
TSNEPlot.1(Monocytes,do.label = T)

##############################
# PrepareGSEA
###############################
Monocytes <- SetAllIdent(Monocytes, id="orig.ident")
PrepareGSEA(Monocytes, k = 100)

##############################
# Run GSEA and generate reports
###############################
GSEA <- readr::read_delim("output/20190220/gsea_report_for_Control_go.xls",
                      "\t", escape_double = FALSE, trim_ws = TRUE)
GSEA %>% head(40) %>% kable() %>% kable_styling()
(gsea_path <- paste0("~/gsea_home/output/",tolower(format(Sys.Date(), "%b%d")), 
                     "/c5.all.Control_versus_NAM.Gsea.1550700096840"))
(c5.all.Control_versus_NAM <- sapply(GSEA$NAME[1:9], function(name) {
        paste0("enplot_",name, "_([0-9]+)*\\.png$")}) %>%
                sapply(function(x) list.files(path = gsea_path, pattern =x)) %>%
                .[sapply(.,length)>0] %>% #Remove empty elements from list with character(0)
                paste(gsea_path, ., sep = "/")) 
CombPngs(c5.all.Control_versus_NAM, ncol = 3)
