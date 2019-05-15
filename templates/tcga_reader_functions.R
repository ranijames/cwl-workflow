### reader functions for TCGA data

source("TUPRO_BASEDIR/hypoxia/R_scripts/DeSeq_normalization.R")
source("TUPRO_BASEDIR/hypoxia/R_scripts/DeSeq_DE.R")

library("data.table")
require("RColorBrewer")
library("DESeq2")
library("ggplot2")
library("ggsignif")
library("pheatmap")
library("genefilter")
library("plyr")
library("Rtsne")
library("dml")
library("ROCR")
library("lfda")
library("kernlab")
library("scales")

is.nan.data.frame <- function(x){ do.call(cbind, lapply(x, is.nan)) }

options(digits=6)

translate_ENS_to_HGNC <- function(in_table){
	print(head(in_table))
	in_table$id = remove_period(in_table$id)
  annot_table = data.frame(fread("/cluster/work/grlab/projects/Fagin/Project_06525_Update/annotation/GRCh38_p10_ens_to_trans.txt",
    header=T), check.names=F)
  colnames(annot_table) = c("id", "ensembl_trans_id", "hgnc_symbol")
  annot_table = unique(annot_table[,c(1,3)])

  in_table = merge(annot_table, in_table, by="id", all.y=T)
  in_table$hgnc_symbol[is.na(in_table$hgnc_symbol)] = ""
  in_table$hgnc_symbol[in_table$hgnc_symbol == ""] = in_table$id[in_table$hgnc_symbol == ""]
  return(in_table)

}


translate_HGNC_to_ENS <- function(in_table){
	print(head(in_table))
  annot_table = data.frame(fread("/cluster/work/grlab/projects/Fagin/Project_06525_Update/annotation/GRCh38_p10_ens_to_trans.txt",
    header=T), check.names=F)
  colnames(annot_table) = c("ensembl_gene_id", "ensembl_trans_id", "id")
  annot_table = unique(annot_table[,c(1,3)])

  in_table = merge(annot_table, in_table, by="id", all.y=T)
  in_table$ensembl_gene_id[is.na(in_table$ensembl_gene_id)] = ""
  in_table$ensembl_gene_id[in_table$ensembl_gene_id == ""] = in_table$id[in_table$ensembl_gene_id == ""]
  return(in_table)

}


get_expr_psi_data_total_gtex <- function(){

  # get the HIF gene expr values
  hif_gene_rand = data.frame(fread(HIF_GENE_SIG_GTEX_FILE))

  # get the GTEX PSI values
  psi_tcga_table = data.frame(fread(PSI_GTEX_FILE), check.names=F, stringsAsFactors=F)
  psi_tcga_table = psi_tcga_table[,c(1,13:ncol(psi_tcga_table))]

  psi_tcga_table_melt = melt(psi_tcga_table, id.vars=c("event_id"))
  psi_tcga_table_melt$variable = gsub("[.].*", "", psi_tcga_table_melt$variable)

  colnames(psi_tcga_table_melt) = c("event_id", "tcga_id", "psi")

  psi_tcga_table_cast = dcast(unique(psi_tcga_table_melt), tcga_id~event_id, value.var="psi", fun.aggregate=mean)

  hif_join_data = join(hif_gene_rand, psi_tcga_table_cast)

  # get metadata
  metadata_table = data.frame(fread(GTEX_META_DATA_FILE))
  metadata_table = metadata_table[,c("Run_s", "body_site_s", "sex_s")]
  metadata_table = unique(metadata_table)
  colnames(metadata_table) = c("tcga_id", "tissueSpecific", "gender")
  hif_join_data = join(metadata_table, hif_join_data, type="right")

  return(hif_join_data)
}


get_expr_psi_data_total <- function(){

  # get the HIF gene expr values
  hif_gene_sig = data.frame(fread(HIF_GENE_SIG_FILE))

  # get the TCGA PSI values
  psi_tcga_table = data.frame(fread(PSI_TCGA_FILE), check.names=F, stringsAsFactors=F)
  psi_gene_map = psi_tcga_table[,c(1,4)]
  psi_gene_map$gene_name = gsub("[.].*", "", psi_gene_map$gene_name)
  psi_tcga_table = psi_tcga_table[,c(1,13:ncol(psi_tcga_table))]

  # get the exon_skip 926 change
  psi_gene_map = rbind(psi_gene_map, c("exon_skip_926B", "ENSG00000117620"))
  psi_926_table = data.frame(fread(EXON_SKIP_926_PATH), check.names=F, stringsAsFactors=F)
  psi_926_table = psi_926_table[psi_926_table$exon_start == 100435992,8:ncol(psi_926_table)]
  psi_926_table = log10(apply(psi_926_table, 2, sum)+1)
  psi_926_table = data.frame(event_id="exon_skip_926B", t(psi_926_table), check.names=F, stringsAsFactors=F)
  psi_tcga_table = rbind(psi_tcga_table, psi_926_table)

  psi_tcga_table_melt = melt(psi_tcga_table, id.vars=c("event_id"))
  psi_tcga_table_melt$variable = gsub("[.].*", "", psi_tcga_table_melt$variable)

  colnames(psi_tcga_table_melt) = c("event_id", "tcga_id", "psi")

  psi_tcga_table_cast = dcast(unique(psi_tcga_table_melt), tcga_id~event_id, value.var="psi", fun.aggregate=mean)

  hif_join_data = join(hif_gene_sig, psi_tcga_table_cast)

	hif_join_data$tcga_id = substr(hif_join_data$tcga_id, 1, 16)



  # get the gender
  survival_table = read.table(GENDER_FILE, sep="\t", header=T, fill=T)
  col_interest = c("gender", "bcr_patient_barcode")
  survival_table = survival_table[,col_interest]
  colnames(survival_table)[2] = "tcga_id_mini"

  orig_colnames = colnames(hif_join_data)
  hif_join_data$tcga_id_mini = substr(hif_join_data$tcga_id, 1, 12)
  hif_join_data = join(survival_table, hif_join_data, type="right")

  hif_join_data = hif_join_data[,c("gender", orig_colnames)]

  hif_join_data = hif_join_data[!is.na(hif_join_data$tcga_id), ]

  return(list(hif_join_data, psi_gene_map))
}
