# scRNAseq-LungCancer
single cell RNA-Seq analysis of human lung cancer. Data from Vivek Mittal <vim2010@med.cornell.edu>  and Michael J. Crowley <mic2029@med.cornell.edu>

## How to run Rscript

### 1. Quality control and sample alignment
https://github.com/nyuhuyang/scRNAseq-LungCancer/blob/master/R/QC.R <br />
`Rscript R/QC.R doc/190212_scRNAseq_info.xlsx` <br />
This script generates ggplot figures for raw data.

https://github.com/nyuhuyang/scRNAseq-LungCancer/blob/master/R/scater.R <br />
`Rscript R/scater.R doc/190212_scRNAseq_info.xlsx` <br />
This script will run scater pipeline and generate one sce object stored in data folder.

https://github.com/nyuhuyang/scRNAseq-LungCancer/blob/master/R/Seurat_setup.R <br />
`Rscript R/Seurat_setup.R doc/190212_scRNAseq_info.xlsx` <br />
This script will run Seurat pipeline after scater and generate figures stored in output folder and Seurat object as LungCancer_{number}_{date}.Rda in data folder. The seruat object name might change

### 2. Identify cell types
https://github.com/nyuhuyang/scRNAseq-LungCancer/blob/master/R/Identify_Cell_Types_Manually.R <br />
`Rscript R/Identify_Cell_Types_Manually.R data/LungCancer_4_20190305.Rda` <br />
This script generates marker Feature plots stored in output folder.

https://github.com/nyuhuyang/scRNAseq-LungCancer/blob/master/R/SingleR.R <br />
`Rscript R/SingleR.R data/LungCancer_4_20190305.Rda` <br />
This script generates SingleR object as LungCancer_{number}_{date}.Rda and store in output folder. File name might change.

https://github.com/nyuhuyang/scRNAseq-LungCancer/blob/master/R/SingleR_figures.R <br />
`Rscript R/SingleR_figures.R data/LungCancer_4_20190305.Rda output/singler_T_LungCancer_4_20190305.Rda`<br />
This script generates TSNEplots stored in output folder. File name LungCancer_{number}_{date}.Rda might change.

https://github.com/nyuhuyang/scRNAseq-LungCancer/blob/master/R/Subset_SingleR_figures.R <br />
`Rscript R/Subset_SingleR_figures.R data/LungCancer_4_20190305.Rda  doc/190212_scRNAseq_info.xlsx` <br />
This script generates subset TSNEplots stored in output folder. File name LungCancer_{number}_{date}.Rda might change.
