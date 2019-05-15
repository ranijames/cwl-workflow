## HIF experiment to run with new patients
suppressWarnings(suppressMessages(library(DESeq)))
suppressWarnings(suppressMessages(library(lattice)))
suppressWarnings(suppressMessages(library(BiocGenerics)))
suppressWarnings(suppressMessages(library(easyRNASeq)))
suppressWarnings(suppressMessages(library(Biobase)))
suppressWarnings(suppressMessages(library(VennDiagram)))

source("/cluster/home/aalva/Projects/Kjo_proct/projects2018_tumorprofiler/hypoxia/R_scripts/tcga_reader_functions.R")
source("/cluster/home/aalva/Projects/Kjo_proct/projects2018_tumorprofiler/hypoxia/R_scripts/DCA_methods.R")
source("/cluster/home/aalva/Projects/Kjo_proct/projects2018_tumorprofiler/hypoxia/R_scripts/set_r_paths.R")
source("/cluster/home/aalva/Projects/Kjo_proct/projects2018_tumorprofiler/hypoxia/R_scripts/DeSeq_normalization.R")
source("/cluster/home/aalva/Projects/Kjo_proct/projects2018_tumorprofiler/hypoxia/R_scripts/tcga_reader_functions.R")

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
library(stringr)

toChar <- function(vec){
  return(as.character(unlist(vec)))
}

get_normalization_factor <- function(){
  lib_size_table = data.frame(fread(HDF5_LIBRARY_SIZE_FILE))
  tcga_norm_factor = max(lib_size_table$libsize_75percent)

  lib_size_table = data.frame(fread(HDF5_LIBRARY_SIZE_GTEX_FILE))
  gtex_norm_factor = max(lib_size_table$libsize_75percent)

  return(c(tcga_norm_factor, gtex_norm_factor))
}

plot_boxplot_newsamp_expr <- function(plot_df, genes_interest, curr_ct, outdir, outstr){

  plot_df = plot_df[plot_df$study == curr_ct,]
  plot_df[,genes_interest] = t(scale(log10(t(plot_df[,genes_interest]+1))))


  plot_df_melt = melt(plot_df[,c(genes_interest, "is_new_sample", "is_normal")], id=c("is_normal", "is_new_sample"))
  outfile = paste(outdir, outstr, "_expr_boxplot.pdf", sep="")
  pdf(outfile)
    gg = ggplot(plot_df_melt, aes(x=is_normal, y=log10(value+1), fill=is_new_sample)) +
          geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(gg)

  dev.off()

  plot_df = plot_df[plot_df$study == curr_ct,]
  plot_df_melt = melt(plot_df[,c(genes_interest, "is_new_sample", "sample_id")], id=c("sample_id", "is_new_sample"))

  mean_val = by(plot_df_melt$value, plot_df_melt$sample_id, median)
  ranked_samps = names(mean_val)[order(mean_val)]
  plot_df_melt$sample_id = factor(plot_df_melt$sample_id, levels=ranked_samps)

  outfile = paste(outdir, outstr, "_expr_boxplot_samp.pdf", sep="")
  pdf(outfile, height=6, width=50)
    gg = ggplot(plot_df_melt, aes(x=sample_id, y=value, fill=is_new_sample)) +
          geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(gg)

  dev.off()

}

make_patients_table <- function(genes_interest, tcga_norm_factor){


  #patient_files = grep("counts.tsv", list.files(NEW_PATIENTS_DIR, recursive=T, full.names=T), value=T)
  #libsize_files = grep("libsize.tsv", list.files(NEW_PATIENTS_DIR, recursive=T, full.names=T), value=T)
  patient_files = count_expression
  libsize_files = lib_size
  all_patients = NA
  for(curr_file_idx in 1:length(patient_files)){
    curr_file = patient_files[curr_file_idx]
    libsize_file = libsize_files[curr_file_idx]

    curr_patient = data.frame(fread(curr_file))
    patient_id = basename(curr_file)
    patient_id = stringi::stri_extract(patient_id,regex = "[^__]+")

    libsize = data.frame(fread(libsize_file))
    patient_id2 = basename(libsize_file)
    patient_id2 = stringi::stri_extract(patient_id2,regex = "[^__]+")

    # make sure we are using the right library size file
    if(patient_id != patient_id2){
      print(patient_id)
      print(patient_id2)
      print("WRONG LIBSIZE")
      return(1)
    }

    colnames(curr_patient) = c("gene_id", patient_id)
    curr_patient$gene_id = remove_period(curr_patient$gene_id)
    curr_patient[,2] = curr_patient[,2]/(libsize$libsize_75percent/tcga_norm_factor)
    print(libsize)
    print(curr_file)



    if(is.na(all_patients)){
      all_patients = curr_patient
    }else{
      all_patients = join(all_patients, curr_patient)
    }
  }

  all_patients = all_patients[all_patients$gene_id %in% genes_interest,]
  all_patients_t = data.frame(t(all_patients[,2:ncol(all_patients)]))
  colnames(all_patients_t) = all_patients$gene_id

  ## Condition for sample_ID, if the input sample_id is eq=1 or greater than one
  if (row.names(all_patients_t)==1){
     all_patients_t$sample_id = sample_id
  }else{
    all_patients_t$sample_id = row.names(all_patients_t)
    all_patients_t$sample_id = str_split_fixed(all_patients_t$sample_id, "-", 2)[,1]
    }

  all_patients_t$is_normal = FALSE
  all_patients_t$study = cancertype

  return(all_patients_t)

}
read_in_new_patients <- function(hif_full_data, genes_interest, tcga_norm_factor){


  new_patients = make_patients_table(genes_interest, tcga_norm_factor)
  new_patients$is_new_sample = TRUE

  new_patients = new_patients[,c("sample_id", "is_normal", "study", "is_new_sample", genes_interest)]

  sample_id_col = which(colnames(hif_full_data) == "tcga_id")
  colnames(hif_full_data)[sample_id_col] = "sample_id"
  hif_full_data$is_new_sample = FALSE
  hif_full_data = hif_full_data[,c("sample_id", "is_normal", "study", "is_new_sample", genes_interest)]

  hif_data_new = rbind(hif_full_data, new_patients)

  return(hif_data_new)

}

get_prob_hypoxic <- function(hif_data_new, gtex_table, genes_interest, outdir, outstr){

    gtex_table$is_normal = TRUE
    gtex_table$is_new_sample = "gtex"

    sample_id_col = which(colnames(gtex_table) == "tcga_id")
    colnames(gtex_table)[sample_id_col] = "sample_id"

    print(gtex_table[1:10,1:10])
    gtex_full_data_sig = gtex_table[,c("sample_id", "study", "is_normal", "is_new_sample", "tissueSpecific", genes_interest)]
    hif_data_new$tissueSpecific = cancertype
    hif_full_data_sig = hif_data_new[,c("sample_id", "study", "is_normal", "is_new_sample", "tissueSpecific", genes_interest)]

    hif_full_data_sig = rbind(hif_full_data_sig, gtex_full_data_sig)
    hif_full_data_sig$sample_type = paste(hif_full_data_sig$is_normal, hif_full_data_sig$is_new_sample, sep="_")

    # get  w
    response_Y = as.numeric(hif_full_data_sig$is_normal)
    ct_try = c(cancertype, "HNSC", "PRAD", "BRCA", "KIRP", "THCA", "KIRC", "BLCA", "ESCA", "PAAD", "UCEC", "LIHC", "STAD", "SARC", "LUSC", "LUAD", "COAD", "READ", "KICH", "CHOL")
    train_idx = which(hif_full_data_sig$study %in% ct_try & hif_full_data_sig$is_new_sample == FALSE)
    test_idx =  which(hif_full_data_sig$study %in% cancertype)

    dca_test_res = do_dca_proj(hif_full_data_sig, genes_interest, train_idx=train_idx,
                      test_idx=test_idx, target_values=response_Y)

    dca_proj = dca_test_res[[2]]
    #perf_roc_test = dca_test_res[[5]]
    #perf_pr_test = dca_test_res[[6]]

    # merge projection estimate with all the other data
    hif_full_data_test = hif_full_data_sig[test_idx,]
    hif_full_data_test$row_id = row.names(hif_full_data_test)
    dca_proj = join(dca_proj[,c("dca_res", "row_id", "Y")], hif_full_data_test)
    dca_proj = dca_proj[order(dca_proj$dca_res),]
    row.names(dca_proj) = dca_proj$row_id


    # plot the histograms of GTEx vs ICGC
    dca_proj_plot = dca_proj[dca_proj$study == cancertype & dca_proj$is_new_sample != TRUE,]

    #filter out the one normal TCGA sample
    dca_proj_plot = dca_proj_plot[which(!(dca_proj_plot$is_new_sample == FALSE & dca_proj_plot$is_normal == TRUE)),]

    # get the scores for all the new samples
    new_samp_idx = dca_proj$dca_res[dca_proj$is_new_sample == TRUE]

    #outfile = paste(outdir, outstr, "_new_patient_distr.pdf", sep="")
    #pdf(outfile)
      #gg = ggplot(dca_proj_plot, aes(x=dca_res, fill=is_normal)) +
            #geom_histogram(alpha=0.5, position="identity") +
            #geom_vline(xintercept=new_samp_idx)
      #print(gg)
      #gg = ggplot(dca_proj_plot, aes(x=dca_res, fill=is_normal)) +
            #geom_density() +
            #geom_vline(xintercept=new_samp_idx)
      #print(gg)
      #gg = ggplot(dca_proj_plot, aes(x=dca_res, fill=tissueSpecific)) +
            #geom_density() +
            #geom_vline(xintercept=new_samp_idx)
      #print(gg)

    #dev.off()

    # calculate odds ratio
    mean_hyp = mean(dca_proj$dca_res[dca_proj$is_new_sample == FALSE])
    sd_hyp = sd(dca_proj$dca_res[dca_proj$is_new_sample == FALSE])

    mean_norm = mean(dca_proj$dca_res[dca_proj$is_new_sample == "gtex"])
    sd_norm = sd(dca_proj$dca_res[dca_proj$is_new_sample == "gtex"])

    liklihood_hyp = dnorm(new_samp_idx, mean_hyp, sd_hyp)
    liklihood_norm = dnorm(new_samp_idx, mean_norm, sd_norm)

    bayes_factor = liklihood_hyp/liklihood_norm
    2*log(bayes_factor)
    print(bayes_factor)
    print(2*log(bayes_factor))

    df_sample = data.frame(dca_res=c(rnorm(1000, mean_hyp, sd_hyp), rnorm(1000, mean_norm, sd_norm)),
                            is_normal = c(rep("hyp", 1000), rep("norm", 1000)))
    #outfile = paste(outdir, outstr, "_estimated_distr.pdf", sep="")
    #pdf(outfile)
      #gg = ggplot(df_sample, aes(x=dca_res, fill=is_normal)) +
            #geom_histogram(alpha=0.5, position="identity") +
            #geom_vline(xintercept=new_samp_idx)
      #print(gg)
      #gg = ggplot(df_sample, aes(x=dca_res, fill=is_normal)) +
            #geom_density() +
            #geom_vline(xintercept=new_samp_idx)
      #print(gg)

    #dev.off()


}


project_new_patients_plot_qc <- function(hif_data_new, gtex_table, genes_interest, outdir, outstr){

  gtex_table$is_normal = TRUE
  gtex_table$is_new_sample = "gtex"

  sample_id_col = which(colnames(gtex_table) == "tcga_id")
  colnames(gtex_table)[sample_id_col] = "sample_id"

  print(gtex_table[1:10,1:10])
  gtex_full_data_sig = gtex_table[,c("sample_id", "study", "is_normal", "is_new_sample", genes_interest)]
  hif_full_data_sig = hif_data_new[,c("sample_id", "study", "is_normal", "is_new_sample", genes_interest)]

  hif_full_data_sig = rbind(hif_full_data_sig, gtex_full_data_sig)
  hif_full_data_sig$sample_type = paste(hif_full_data_sig$is_normal, hif_full_data_sig$is_new_sample, sep="_")


  hif_data_plot = hif_full_data_sig[hif_full_data_sig$study == cancertype, ]
  color_labels = hif_data_plot$is_new_sample
  #plot_tsne(hif_data_plot[,genes_interest], color_labels, outdir, outstr=cancertype)
  #plot_pca_dca(hif_data_plot[,genes_interest], color_labels, outdir, outstr=cancertype)
  #plot_boxplot_newsamp_expr(hif_data_plot, genes_interest, curr_ct=cancertype, outdir, outstr=cancertype)

  hif_data_plot = hif_full_data_sig[hif_full_data_sig$study == cancertype & hif_full_data_sig$is_normal==FALSE, ]
  color_labels = hif_data_plot$is_new_sample
  #plot_pca_dca(hif_data_plot[,genes_interest], color_labels, outdir, outstr=cancertype"_tumor")


  # get  w
  response_Y = as.numeric(hif_full_data_sig$is_normal)
  ct_try = c(cancertype, "HNSC", "PRAD", "BRCA", "KIRP", "THCA", "KIRC", "BLCA", "ESCA", "PAAD", "UCEC", "LIHC", "STAD", "SARC", "LUSC", "LUAD", "COAD", "READ", "KICH", "CHOL")
  train_idx = which(hif_full_data_sig$study %in% ct_try & hif_full_data_sig$is_new_sample == FALSE)
  test_idx =  which(hif_full_data_sig$study %in% cancertype)

  dca_test_res = do_dca_proj(hif_full_data_sig, genes_interest, train_idx=train_idx,
                    test_idx=test_idx, target_values=response_Y)

  dca_proj = dca_test_res[[2]]
  #perf_roc_test = dca_test_res[[5]]
  #perf_pr_test = dca_test_res[[6]]

  dca_train_res = do_dca_proj(hif_full_data_sig, genes_interest, train_idx=train_idx,
                    test_idx=train_idx, target_values=response_Y)

  #perf_roc_train = dca_train_res[[5]]
  #perf_pr_train = dca_train_res[[6]]

  brca_test_res = do_dca_proj(hif_full_data_sig, genes_interest, train_idx=which(!hif_full_data_sig$study %in% c("BRCA")),
                    test_idx=which(hif_full_data_sig$study %in% c("BRCA")), target_values=response_Y)

  #perf_brca_roc_test = brca_test_res[[5]]
  #perf_brca_pr_test = brca_test_res[[6]]


  #outfile = paste(outdir, outstr, "_full_ROC.pdf", sep="")
  #pdf(outfile)
    #plot(perf_roc_train, col="purple", main=paste("train:", dca_train_res[[3]], "test:", dca_test_res[[3]]), ylim=c(0,1))
    #plot(perf_roc_test, col="orange", add=T)
    #plot(perf_pr_train, col="purple", main=paste("train:", dca_train_res[[4]], "test:", dca_test_res[[4]]), ylim=c(0,1))
    #plot(perf_pr_test, col="orange", add=T)

    #plot(perf_roc_train, col="purple", main=paste("noBRCA train:", dca_train_res[[3]], "BRCA test:", brca_test_res[[3]]), ylim=c(0,1))
    #plot(perf_brca_roc_test, col="orange", add=T)
    #plot(perf_pr_train, col="purple", main=paste("noBRCA train:", dca_train_res[[4]], "BRCA test:", brca_test_res[[4]]), ylim=c(0,1))
    #plot(perf_brca_pr_test, col="orange", add=T)

  #dev.off()

  hif_full_data_test = hif_full_data_sig[test_idx,]
  hif_full_data_test$row_id = row.names(hif_full_data_test)
  dca_proj = join(dca_proj[,c("dca_res", "row_id", "Y")], hif_full_data_test)
  dca_proj = dca_proj[order(dca_proj$dca_res),]
  row.names(dca_proj) = dca_proj$row_id


  # plot the projection
  #outfile = paste(outdir, outstr, "_new_patient.pdf", sep="")
  #pdf(outfile)
    #gg = ggplot(dca_proj, aes(x=dca_res, fill=is_normal)) +
          geom_histogram(alpha=0.5, position="identity")
    #print(gg)
    #gg = ggplot(dca_proj, aes(x=is_new_sample, y=dca_res, fill=is_normal)) +
        #  geom_boxplot() + geom_jitter()+
          #theme(axis.text.x = element_text(angle = 90, hjust = 1))
    #print(gg)

  #dev.off()


  # plot all cancer types
  dca_test_res = do_dca_proj(hif_full_data_sig, genes_interest, train_idx=train_idx,
                    test_idx=1:nrow(hif_full_data_sig), target_values=response_Y)
  dca_proj_full = dca_test_res[[2]]

  hif_full_data_test_all = hif_full_data_sig
  hif_full_data_test_all$row_id = row.names(hif_full_data_test_all)
  dca_proj_full = join(dca_proj_full[,c("dca_res", "row_id", "Y")], hif_full_data_test_all)
  dca_proj_full = dca_proj_full[order(dca_proj_full$dca_res),]
  row.names(dca_proj_full) = dca_proj_full$row_id


  # plot the projection with all samples
  dca_proj_full$is_new_sample[dca_proj_full$is_new_sample==TRUE] = dca_proj_full$sample_id[dca_proj_full$is_new_sample==TRUE]
  #outfile = paste(outdir, outstr, "_new_patient_all.pdf", sep="")
  #pdf(outfile)
    #gg = ggplot(dca_proj_full, aes(x=dca_res, fill=is_new_sample)) +
          #geom_histogram(alpha=0.5, position="identity")
    #print(gg)
    #gg = ggplot(dca_proj_full[dca_proj_full$study %in% c("BLCA", "BRCA", "HNSC", "LUAD", "LUSC", cancertype),], aes(x=study, y=dca_res, fill=is_new_sample)) +
          #geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    #print(gg)

  #dev.off()

  # make a heatmap
  compare_table = dca_proj_full[(dca_proj_full$study %in% c("BLCA", "BRCA", "HNSC", "LUAD", "LUSC") | dca_proj_full$is_new_sample==TRUE),]
  print(unique(compare_table$is_new_sample))
  compare_table = compare_table[order(compare_table$dca_res),]

  # plot the heatmap version
  annot_row = data.frame(is_new=as.factor(compare_table$is_new_sample), dca_score=compare_table$dca_res)
  row.names(annot_row) = row.names(compare_table)
  compare_table = compare_table[order(compare_table$dca_res),]

  #outfile = paste(outdir, "test_heatmap.pdf", sep="")
  #pdf(outfile)
  #pheatmap(log10(compare_table[,genes_interest]+1), show_rownames=F, show_colnames=T,
            #annotation_row=annot_row, fontsize_col=3, clustering_method="ward.D2",
            #cluster_rows=T)
  #pheatmap(log10(compare_table[,genes_interest]+1), show_rownames=F, show_colnames=T,
            #annotation_row=annot_row, fontsize_col=3, clustering_method="ward.D2",
            #cluster_rows=F)
  #dev.off()

  return(dca_proj_full)

}

do_cancertype_brca_dca_TN_experiment <- function(genes_interest, hif_full_data, gtex_table, outstr, outdir){

  gtex_table$is_normal = TRUE
  gtex_table$is_new_sample = "gtex"

  sample_id_col = which(colnames(gtex_table) == "tcga_id")
  colnames(gtex_table)[sample_id_col] = "sample_id"

  print(gtex_table[1:10,1:10])
  gtex_full_data_sig = gtex_table[,c("sample_id", "study", "is_normal", "is_new_sample", genes_interest)]
  hif_full_data_sig = hif_data_new[,c("sample_id", "study", "is_normal", "is_new_sample", genes_interest)]

  hif_full_data_sig = rbind(hif_full_data_sig, gtex_full_data_sig)
  hif_full_data_sig$sample_type = paste(hif_full_data_sig$is_normal, hif_full_data_sig$is_new_sample, sep="_")

  response_Y = as.numeric(hif_full_data_sig$is_normal)
  ct_try = c(cancertype, "HNSC", "PRAD", "BRCA", "KIRP", "THCA", "KIRC", "BLCA", "ESCA", "PAAD", "UCEC", "LIHC", "STAD", "SARC", "LUSC", "LUAD", "COAD", "READ", "KICH", "CHOL")
  train_idx = which(hif_full_data_sig$study %in% ct_try & hif_full_data_sig$is_new_sample == FALSE)
  test_idx =  1:nrow(hif_full_data_sig)

  dca_brca_res = do_dca_proj(hif_full_data_sig, genes_interest, train_idx=train_idx,
                    test_idx=test_idx, target_values=response_Y)

  dca_proj = dca_brca_res[[2]]



  # get metadata
  hif_full_data_BRCA = hif_full_data_sig
  hif_full_data_BRCA$row_id = row.names(hif_full_data_BRCA)
  dca_brca_res = join(dca_proj[,c("dca_res", "row_id", "Y")], hif_full_data_BRCA)

  # merge with the survival data to get TN score
	survival_table = read.delim(survival_file, sep="\t", header=T, fill=T)
	col_interest = c("days_to_last_known_alive", "bcr_patient_uuid", "bcr_patient_barcode",
									"acronym", "vital_status", "days_to_death", "project_code",
									"disease_code", "pathologic_stage", "clinical_stage",
									"patient_death_reason", "breast_carcinoma_estrogen_receptor_status",
									"breast_carcinoma_progesterone_receptor_status", "person_neoplasm_cancer_status")

	survival_table = survival_table[,col_interest]
	colnames(survival_table)[4] = "cancer_type"
	colnames(survival_table)[3] = "tcga_id_mini"
  survival_table$triple_neg = survival_table$breast_carcinoma_estrogen_receptor_status=="Negative" &
    survival_table$breast_carcinoma_progesterone_receptor_status=="Negative"

  dca_brca_res$tcga_id_mini = substr(dca_brca_res$sample_id, 1, 12)
  dca_brca_res = join(dca_brca_res, survival_table, type="left")
  dca_brca_res[is.na(dca_brca_res$triple_neg) | dca_brca_res$triple_neg==FALSE, "triple_neg"] = "not_triple_neg"
  dca_brca_res[dca_brca_res$triple_neg=="not_triple_neg" & dca_brca_res$study==cancertype, "triple_neg"] = paste(cancertype,"_samples")
  dca_brca_res[dca_brca_res$is_normal %in% c("gtex", TRUE), "triple_neg"] = "normal"
  dca_brca_res[dca_brca_res$triple_neg==TRUE, "triple_neg"] = "triple_neg"

  dca_brca_res$triple_neg = factor(dca_brca_res$triple_neg, levels=c("normal", "not_triple_neg", "triple_neg", cancertype,"_samples"))
  dca_full_res = dca_brca_res
  dca_brca_res = dca_brca_res[dca_brca_res$study=="BRCA" | dca_brca_res$is_new_sample==TRUE,]


  #outfile = paste(outdir, "new_dist_BRCA_in_", outstr, "_space_TN.pdf", sep="")
  #pdf(outfile)
    #gg = ggplot(dca_brca_res, aes(x=dca_res, fill=is_normal)) +
          #geom_histogram(alpha=0.5, position="identity")
    #print(gg)
    #gg = ggplot(dca_brca_res, aes(x=is_normal, y=dca_res, fill=is_normal)) +
          #geom_boxplot()
    #print(gg)
    #gg = ggplot(dca_brca_res, aes(x=dca_res, fill=triple_neg)) +
          #geom_histogram(alpha=0.5, position="identity")
    #print(gg)
    #gg = ggplot(dca_brca_res, aes(x=triple_neg, y=dca_res, fill=triple_neg)) +
          #geom_boxplot() +
          #geom_signif(comparisons = list(c("normal", "not_triple_neg"), c("not_triple_neg", "triple_neg")),
              #map_signif_level=TRUE)
    #print(gg)

  #dev.off()

  return(dca_full_res)

}

outdir = paste(RESULTS_FOLDER, sep="")

res = get_expr_psi_data_total()
psi_gene_map = res[[2]]
hif_full_data = res[[1]]

genes_interest = unique(c(HIF_GENES))
genes_interest = genes_interest[genes_interest %in% colnames(hif_full_data)]


# now only get the samples that have atleast 4 normals and SARC
normal_freq = table(hif_full_data[,c("study", "is_normal")])
normal_freq = data.frame(normal_freq)
normal_freq = normal_freq[(normal_freq$Freq > 3 & normal_freq$is_normal==TRUE) |
                (normal_freq$Freq > 3 & normal_freq$is_normal==FALSE),]

normal_freq_pass = normal_freq$study[duplicated(normal_freq$study)]
hif_full_data = hif_full_data[hif_full_data$study %in% c(cancertype, toChar(normal_freq_pass)),]


hif_full_data = unique(hif_full_data)


# get gtex data
gtex_res = get_expr_psi_data_total_gtex()
colnames(gtex_res)[4] = "tissues"

# remove the transformed fibroblasts from skin
gtex_res = gtex_res[gtex_res$tissueSpecific != "Cells - Transformed fibroblasts", ]
gtex_res = gtex_res[gtex_res$tissueSpecific != "Cells - Transformed fibroblasts", ]
gtex_res = gtex_res[gtex_res$tissueSpecific != "Cells - Leukemia cell line (CML)", ]



tissues = c("Esophagus", "Breast", "Thyroid", "Kidney", "Kidney", "Kidney", "Pancreas", "Bladder", "Brain", "Colon", "Liver", "Lung", "Lung", "Ovary", "Prostate", "Skin", "Uterus", "Stomach")
study = c("ESCA", "BRCA", "THCA", "KICH", "KIRC", "KIRP", "PAAD", "BLCA", "GBM", "COAD", "LIHC", "LUAD", "LUSC", "OV", "PRAD", cancertype, "UCEC", "STAD")
cancer_tissue_map = data.frame(tissues, study)

gtex_table = merge(cancer_tissue_map, gtex_res, by=c("tissues"))

print(gtex_res[1:10,1:10])

# get projection vector
norm_factors = get_normalization_factor()
tcga_norm_factor = norm_factors[1]
gtex_norm_factor = norm_factors[2]
hif_data_new = read_in_new_patients(hif_full_data, genes_interest, tcga_norm_factor)

# filter for gene expression cuttoff
genes_interest = colnames(hif_data_new)[which(colnames(hif_data_new) %in% genes_interest)]
mean_expr = apply(hif_data_new[,genes_interest], 2, quantile, 0.25)
mean_expr = mean_expr[mean_expr > 1000]
genes_interest = genes_interest[genes_interest %in% names(mean_expr)]


#plot_boxplot_expr(hif_data_new, genes_interest, outdir, outstr="test")

#plot_boxplot_newsamp_expr(hif_data_new, genes_interest, cancertype, outdir, outstr="new")

#color_labels = hif_data_new[hif_data_new$study==cancertype,]
#color_labels = paste(color_labels$is_normal, color_labels$is_new_sample, sep="_")
#plot_pca_dca(hif_data_new[hif_data_new$study==cancertype, genes_interest], color_labels, outdir, outstr="new_pca")

dca_all_samples = project_new_patients_plot_qc(hif_data_new, gtex_table, genes_interest, outdir, outstr="test")
dca_all_samples = do_cancertype_brca_dca_TN_experiment(genes_interest, hif_data_new, gtex_table, outstr="test", outdir)

formatted_res_table = dca_all_samples[,c("sample_id", "dca_res", "study", "is_normal", "is_new_sample", "triple_neg")]
colnames(formatted_res_table) = c("sample_id", "hypoxia_score", "cancer_type_match", "is_normal", "study", "BRCA_ONLY_is_triple_neg")
formatted_res_table$BRCA_ONLY_is_triple_neg[formatted_res_table$cancer_type_match != "BRCA"] = "NA"
formatted_res_table_non_new = formatted_res_table[formatted_res_table$study %in% c("gtex", FALSE),]
formatted_res_table_non_new$study[formatted_res_table_non_new$study==FALSE] = "tcga"

formatted_res_table_new = formatted_res_table[!formatted_res_table$study %in% c("gtex", FALSE),]
formatted_res_table_new$study = "tumor_profiler"

write.table(formatted_res_table_non_new, paste(outdir, sample_id,labkey_prefix,"__tcga_gtex_hypoxia_scores.tsv", sep=""), sep="\t", row.names=F, quote=F)
formatted_res_table_new= formatted_res_table_new[grepl(sample_id, formatted_res_table_new$sample_id)]
write.table(formatted_res_table_new, paste(outdir, sample_id,labkey_prefix,"__hypoxia_scores.tsv", sep=""), sep="\t", row.names=F, quote=F)
