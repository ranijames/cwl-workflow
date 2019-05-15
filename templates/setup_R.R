list.of.packages <- c("data.table", "RColorBrewer","ggplot2","ggsignif","pheatmap","plyr","Rtsne","dml","ROCR","lfda","kernlab","scales","VennDiagram","fastICA")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='https://stat.ethz.ch/CRAN/')

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager", repos='https://stat.ethz.ch/CRAN/')
if ("easyRNASeq" %in% rownames(installed.packages()) == FALSE) 
    BiocManager::install("easyRNASeq", version = "3.8")
if ("DESeq2" %in% rownames(installed.packages()) == FALSE) 
    BiocManager::install("DESeq2", version = "3.8")
if ("genefilter" %in% rownames(installed.packages()) == FALSE) 
    BiocManager::install("genefilter", version = "3.8")
if ("biomaRt" %in% rownames(installed.packages()) == FALSE) 
    BiocManager::install("biomaRt", version = "3.8")

