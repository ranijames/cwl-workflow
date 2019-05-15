#source("http://bioconductor.org/biocLite.R")

library(easyRNASeq)
library(DESeq)
library(fastICA)

toChar <- function(vec){
  return(as.character(unlist(vec)))
}

getGeneLengths <- function(org="hsapiens_gene_ensembl"){

  # select mart and data set
  bm <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  bm <- useDataset(org, mart=bm)

  # Get ensembl gene ids and gene length
  # this is our translation table
  id_length <- getBM(mart=bm, attributes=c('ensembl_gene_id', 'cds_length'))
  id_length <- na.omit(id_length)

  #now take the average of the lengths
  mean_length = aggregate(cds_length~ensembl_gene_id, id_length, mean)

  return(mean_length)

}

ica_rank_importance_plot <- function(ica_table, num_of_comp=2, takeLog=T, RESULTS_FOLDER, plot_colors, plot_labels, plot_main, doAbs=TRUE, outname){
  ## whatever ica_table we get, we are going to project them to ICA space
  ## using the eigenvectors V, we will rank the "importance" of the gene on the PC of interest
  ## the rownames must be set, they are returned in the new ranked order

  # first we take the transpose, add one, then do log2
  log.samples <- t(ica_table)
  if(takeLog){
    log.samples <- log2(t(ica_table)+1)   
  } 

  # do projection pursuit
  a_pp <- fastICA(log.samples, num_of_comp, alg.typ = "parallel", fun = "logcosh", alpha = 1, method = "R", row.norm = FALSE, maxit = 200, tol = 0.0001, verbose = TRUE)
  pdf(paste(RESULTS_FOLDER, "/", outname, "_proj_pursuit_", num_of_comp, ".pdf", sep=""))
  par(mfrow = c(1, 3))
  plot(a_pp$X, main = "Pre-processed data", col=plot_colors)
  plot(a_pp$X %*% a_pp$K, main = "PCA components", col=plot_colors)
  plot(a_pp$S, main = "ICA components", col=plot_colors)
  dev.off()

  if(num_of_comp ==2){
    png(paste(RESULTS_FOLDER, "/", outname, "_proj_pursuit_", num_of_comp, ".png", sep=""))
    plot(a_pp$S, main = "ICA components", col=plot_colors)
    dev.off()
  }else{
    pdf(paste(RESULTS_FOLDER, "/", outname, "_proj_pursuit_", num_of_comp, ".pdf", sep=""))
    pairs(a_pp$S, panel = function(x,y) text(x,y, labels=plot_labels,cex=0.5,col=plot_colors), main=plot_main)
    dev.off()

  }


  return(a_pp)

}

ica_rank_importance_write_out_ucs <- function(a_pp, ica_table, num_of_comp=2, comp_of_interest=1, takeLog=T, RESULTS_FOLDER, plot_colors, plot_labels, plot_main, doAbs=TRUE, outname){

  log.samples <- t(ica_table)
  if(takeLog){
    log.samples <- log2(t(ica_table)+1)   
  } 

  if(doAbs){
    weight_PC = abs(a_pp$A[comp_of_interest,])
  }else{
    weight_PC = a_pp$A[comp_of_interest,]
  }

  #now rank
  names(weight_PC) = rownames(ica_table)
  weight_PC = weight_PC[rev(order(weight_PC))]

  #rank samples
  rank_samp = data.frame(rownames(log.samples), a_pp$S[,comp_of_interest])
  colnames(rank_samp) = c("bcr_patient_barcode", "score")
  
  # add in the clinical information
  clinical_ucs = read.table("/cbio/grlab/projects/ucs_davidson/annotation/ucs_clinical_patient.txt", header=T, sep="\t")
  clinical_ucec = read.table("/cbio/grlab/projects/ucs_davidson/annotation/ucec_clinical_patient.txt", header=T, sep="\t")
  clinical_sarc = read.table("/cbio/grlab/projects/ucs_davidson/annotation/sarc_clinical_patient.txt", header=T, sep="\t")

  clinical_ucs = merge(rank_samp, clinical_ucs, by="bcr_patient_barcode")
  write.table(clinical_ucs, paste(RESULTS_FOLDER, "/", outname, "_proj_pursuit_", num_of_comp, "_ic", comp_of_interest, "_ucs_clinical.txt", sep=""), quote=F, col.name=T, row.name=F, sep="\t")

  clinical_ucec = merge(rank_samp, clinical_ucec, by="bcr_patient_barcode")
  write.table(clinical_ucec, paste(RESULTS_FOLDER, "/", outname, "_proj_pursuit_", num_of_comp, "_ic", comp_of_interest, "_ucec_clinical.txt", sep=""), quote=F, col.name=T, row.name=F, sep="\t")

  filename_ucec_corr = paste(RESULTS_FOLDER, "/", outname, "_ucec_corr_", num_of_comp, "_ic", comp_of_interest, "_ucec_clinical.txt", sep="")
  capture.output(cor.test(scale(clinical_ucec$score), clinical_ucec$histological_type_num),file=filename_ucec_corr)

  filename_ucec_stage_corr = paste(RESULTS_FOLDER, "/", outname, "_ucec_stage_corr_", num_of_comp, "_ic", comp_of_interest, "_ucec_stage_clinical.txt", sep="")
  capture.output(cor.test(scale(clinical_ucec$score), clinical_ucec$clinical_stage_num),file=filename_ucec_stage_corr)

  #filename_ucs_stage_corr = paste(RESULTS_FOLDER, "/", outname, "_ucs_stage_corr_", num_of_comp, "_ic", comp_of_interest, "_ucs_stage_clinical.txt", sep="")
  #capture.output(cor.test(scale(clinical_ucs$score), clinical_ucs$clinical_stage_num),file=filename_ucs_stage_corr)


  clinical_sarc = merge(rank_samp, clinical_sarc, by="bcr_patient_barcode")
  write.table(clinical_sarc, paste(RESULTS_FOLDER, "/", outname, "_proj_pursuit_", num_of_comp, "_ic", comp_of_interest, "_sarc_clinical.txt", sep=""), quote=F, col.name=T, row.name=F, sep="\t")


  write.table(rank_samp, paste(RESULTS_FOLDER, "/", outname, "_proj_pursuit_", num_of_comp, "_ic", comp_of_interest, ".txt", sep=""), quote=F, col.name=F, row.name=F)


  return(weight_PC)

}

ica_rank_importance_write_out <- function(a_pp, ica_table, num_of_comp=2, comp_of_interest=1, takeLog=T, RESULTS_FOLDER, plot_colors, plot_labels, plot_main, doAbs=TRUE, outname){

  log.samples <- t(ica_table)
  if(takeLog){
    log.samples <- log2(t(ica_table)+1)   
  } 

  if(doAbs){
    weight_PC = abs(a_pp$A[comp_of_interest,])
  }else{
    weight_PC = a_pp$A[comp_of_interest,]
  }

  #now rank
  names(weight_PC) = rownames(ica_table)
  weight_PC = weight_PC[rev(order(weight_PC))]

  #rank samples
  rank_samp = data.frame(rownames(log.samples), a_pp$S[,comp_of_interest])
  colnames(rank_samp) = c("sample_id", "score")
  

  write.table(rank_samp, paste(RESULTS_FOLDER, "/", outname, "_proj_pursuit_", num_of_comp, "_ic", comp_of_interest, ".txt", sep=""), quote=F, col.name=F, row.name=F)


  return(weight_PC)

}


pca_rank_importance <- function(pca_table, PC_of_interest, takeLog=T, RESULTS_FOLDER, plot_colors, plot_labels, plot_main, doAbs=TRUE){
  ## whatever pca_table we get, we are going to project them to PCA space
  ## using the eigenvectors V, we will rank the "importance" of the gene on the PC of interest
  ## the rownames must be set, they are returned in the new ranked order

  # first we take the transpose, add one, then do log2
  log.samples <- t(pca_table)
  if(takeLog){
    log.samples <- log2(t(pca_table)+1)   
  } 

  

  #now we have the principal component projection
  ir.pca <- prcomp(log.samples, center=F) 

  #what we want is the matrix of eigenvalues, or the rotation matrix
  V = ir.pca$rotation

  #now get the eigenvectors on the PC of interest
  if(doAbs){
    eig_PC = abs(V[,PC_of_interest])
  }else{
    eig_PC = V[,PC_of_interest]
  }

  #now rank
  eig_PC = eig_PC[rev(order(eig_PC))]

  #plot just for checking purposes
  plot_pca(pca_table, paste(RESULTS_FOLDER, "/pca_for_rank.pdf", sep=""), plot_colors, plot_labels, plot_main)

  return(eig_PC)

}

pca_feature_fix_to_point <- function(pca_table, RESULTS_FOLDER, plot_colors, plot_labels, plot_main){

  ## whatever pca_table we get, we are going to project them to PCA space
  ## to get the matrix of the eigenvectors V
  ## once we have that matrix, we assume all of the points are at position 1,1,1,1,1,1,...,1
  ## then use V^(-1) or V transpose (since eigenvectors) to back project to get the original
  ## points
  ## I think this works best with a normalized countsTable

  # first we take the transpose, add one, then do log2
  log.samples <- t(pca_table)
  log.samples <- log2(t(pca_table)+1)   
  

  #now we have the principal component projection
  ir.pca <- prcomp(log.samples, center=F) 

  #what we want is the matrix of eigenvalues, or the rotation matrix
  V = ir.pca$rotation

  #the projected values is x
  x_pca = ir.pca$x
  
  # t(V) * x_pca should equal log.samples
  orig_val = x_pca %*% t(V)

  #we should get back the original matrix
  if(!isTRUE(all.equal(t((2^orig_val)-1), pca_table))){
    print("PCA FEATURE FIX FAILED DO NOT CONTINUE")
    return(NULL)
  }

  #now we take a single point
  #x_new_pca = matrix(c(100), nrow=nrow(x_pca), ncol=ncol(x_pca), byrow=T)

  #now we condense the variance on the second dimension
  x_new_pca = x_pca
  x_new_pca[,2] = 1
  x_new_pca[,3] = 1

  #ok, now we back project using V
  x_new_pca_v = x_new_pca %*% t(V)

  # all of our samples have the SAME expression values for the same gene, but the genes have different
  # expresion values in comparison with each other, which is what you expect IF there is ABSOLUTELY NO
  # transcriptional changes...

  # now just unlog it
  x_new_pca_v_t = t((2^x_new_pca_v)-1)
  #x_new_pca_v_t = t(x_new_pca_v)

  #then add some value to make them non-neg
  x_new_pca_v_t = x_new_pca_v_t + -1 * min(x_new_pca_v_t)

  # then plot it again
  plot_pca(x_new_pca_v_t, paste(RESULTS_FOLDER, "reduced_pca.pdf", sep=""), plot_colors, plot_labels, plot_main)
  summary(x_new_pca_v_t)

 
 
  return(x_new_pca_v_t)

}

get_genes_high_interrep_var <- function(sub_countsTable, quant, RESULTS_FOLDER, plot_name){
  ## this function will return a subset of a table that
  ## has the highest interreplicate variance as calculated 
  ## by top "quant" quantile of the coefficient of variance

  ## this assumes that countsTable is only a table of 1 sample
  ## and all of its replicates.

  ## quant is the quantile that is accepted

  ## RETURNS: two vectors, the row names (gene names) with 
  ## interreplicate variance > quant and < quant

  ## SIDE EFFECTS
  ## diagnostic plots will be put in RESULTS_FOLDER

  sub_counts_var = apply(sub_countsTable, 1, var)
  sub_counts_std = sqrt(sub_counts_var)

  sub_counts_mean = apply(sub_countsTable, 1, mean)

  sub_counts_CoV = na.omit(sub_counts_std/sub_counts_mean)

  pdf(paste(RESULTS_FOLDER, plot_name, "_CoV_hist_before.pdf", sep=""))
  hist(sub_counts_CoV, breaks=200)
  dev.off()


  #remove the top quant% most over dispersed genes
  use = (sub_counts_CoV > quantile(sub_counts_CoV, probs=quant))
  table(use)
  sub_counts_CoV_high = sub_counts_CoV[use]
  sub_counts_CoV_low = sub_counts_CoV[!use]


  pdf(paste(RESULTS_FOLDER, plot_name, "_CoV_hist_after.pdf", sep=""))
  hist(sub_counts_CoV_low, breaks=200)
  dev.off()

  return(list(names(sub_counts_CoV_high), names(sub_counts_CoV_low)))

}

do_deseq_normalization_and_plot_ontable <- function(countsTable, RESULTS_FOLDER, plot_colors, plot_labels, plot_main){
  

  ##normalize the data 
  norm_counts=countsTable

  iqr<-apply(countsTable,1,IQR)
  vars<-apply(countsTable,1,var)
  means<-apply(countsTable,1,mean)

  for ( i in 1:nrow(countsTable)){
    norm_counts[i,]<-(countsTable[i,]-means[i])/(iqr[i]+0.00001)
  }


  #if needed, center the data 
  #norm_counts <- scale(norm_counts,scale = TRUE)


  do_normalization_check_plotting(NULL, norm_counts, countsTable, RESULTS_FOLDER, plot_colors, plot_labels, plot_main, plot_disp=F, doLog=F)
  hierarchical_clust(countsTable, paste(RESULTS_FOLDER, "unnormalized_samp_clust.pdf", sep=""), plot_main)

  return(norm_counts)

}

#filtering out the bottom theta quantile
filter_lowest_counts <- function(countsTable, theta){
  #filter
  rs = rowSums(countsTable )
  use = (rs > quantile(rs, probs=theta))
  table(use)
  countsTable_filtered = countsTable[use,]

  return(countsTable_filtered)

}

filter_thresh_lowest_counts <- function(countsTable, val, theta){

  countsTableFiltered = countsTable[rowSums(countsTable>=val) == ncol(countsTable),]

  #filter
  rs = rowSums(countsTableFiltered )
  use = (rs > quantile(rs, probs=theta))
  table(use)
  countsTable_filtered = countsTable[use,]
  return(countsTable_filtered)
}

filter_threshhold_counts <- function(countsTable, val){
  #ensure that every row has atleast n counts in every column

  countsTableFiltered = countsTable[rowSums(countsTable>=val) == ncol(countsTable),]

  return(countsTableFiltered)
}

#for plotting dispersion
plotDispNorm <- function (cds, name = NULL, ymin,  int_norm_counts, linecol = "#ff000080", xlab = "mean of normalized counts",
ylab = "dispersion", log = "xy", cex = 0.45, ...)
{
	px = rowMeans(int_norm_counts)
	sel = (px > 0)
	px = px[sel]
	py = fitInfo(cds, name = name)$perGeneDispEsts[sel]
	if (missing(ymin))
	ymin = 10^floor(log10(min(py[py > 0], na.rm = TRUE)) -
	0.1)
	plot(px, pmax(py, ymin), xlab = xlab, ylab = ylab, log = log,
	pch = ifelse(py < ymin, 6, 16), cex = cex, ...)
	xg = 10^seq(-0.5, 5, length.out = 100)
	lines(xg, fitInfo(cds, name = name)$dispFun(xg), col = linecol,
	lwd = 4)
}    

#plotting PCA of samples
plot_pca <-function(counting_table, plot_name, plot_colors, plot_labels, plot_main, takeLog=TRUE){
	
  ### plot_name is the name of the PDF, include the path is working directory is not set
  ### plot_colors is a vector of colors consisting of how each sample's color (this is usually based on column names
  ### counting_table is the normalized gene counts in a data frame samples are the columns, genes are the rows
  log.samples <- t(counting_table)
  if(takeLog){
    log.samples <- log2(t(counting_table)+1)   
  } 
  ir.pca <- prcomp(log.samples) 
  vars <- apply(ir.pca$x, 2, var)  
  props <- vars / sum(vars)
  pc1_var_importance = as.character(format(props[1], digits=3))
  pc2_var_importance = as.character(format(props[2], digits=3))
  pc3_var_importance = as.character(format(props[3], digits=3))

  plot_main = paste(plot_main, 'PC1:', pc1_var_importance, 'PC2:', pc2_var_importance, 'PC3:', pc3_var_importance, sep = ' ')

  pdf(plot_name)
  pairs(ir.pca$x[,1:3], panel = function(x,y) text(x,y, labels=plot_labels,cex=0.5,col=plot_colors), main=plot_main)
  dev.off()

  print(ir.pca$x[,1:2])


  return(ir.pca)

}

#plotting PCA of samples
plot_concordance <-function(counting_table, plot_name, plot_colors, plot_labels, plot_main, takeLog=TRUE){
	
  ### plot_name is the name of the PDF, include the path is working directory is not set
  ### plot_colors is a vector of colors consisting of how each sample's color (this is usually based on column names
  ### counting_table is the normalized gene counts in a data frame samples are the columns, genes are the rows


  #also plost concordance of the columns in data.frame
  pdf(plot_name)
  pairs(counting_table, main=plot_main)
  dev.off()


}

#gsub("Sample_", "", colnames(counting_table))

#look at the sample clustering
hierarchical_clust <-function(counting_table, plot_name, plot_main){
	
  ### plot_name is the name of the PDF, include the path is working directory is not set
  ### plot_colors is a vector of colors consisting of how each sample's color (this is usually based on column names
  ### counting_table is the normalized gene counts in a data frame samples are the columns, genes are the rows


  #do ehirarcical clustering
  d <- dist(as.matrix(t(counting_table)))   # find distance matrix 
  hc <- hclust(d)                # apply hirarchical clustering 

  pdf(plot_name)
  plot(hc, main=plot_main)
  dev.off()


}

remove_dash <- function(GeneID_vector){
  myInterestingGenes = as.character(GeneID_vector)
  myInterestingGenes=gsub("[-]", ".", myInterestingGenes)
  return(myInterestingGenes)
}

#make plots to check normalization
do_normalization_check_plotting <- function(cds, norm_counts, counts_table, RESULTS_FOLDER, plot_colors, plot_labels, plot_main, plot_disp=T, doLog=T){

    ### cds is cds object for DeSeq
    ### norm_counts is the counts table after normalization
    ### counts_table is raw counts table
    ### RESULTS_FOLDER is where you want all the plots to go

    #do the original plots
    #make plots to show how the data looks after normalization

    num_genes = nrow(norm_counts)

    #plot the unnormalized boxplots
    pdf(paste(RESULTS_FOLDER, "raw_log2.pdf", sep=""))
    if(doLog){
      boxplot(log2(counts_table+1), las=2, main="Raw log2 Values of Expression", col=plot_colors)
    }else{
      boxplot((counts_table), las=2, main="Raw Values of Expression", col=plot_colors)

    }
    dev.off()

   #plot the normalized boxplots
    pdf(paste(RESULTS_FOLDER, "normalized_log2.pdf", sep=""))
    if(doLog){
      boxplot(log2(norm_counts+1), las=2, main=paste("Normalized log2 Values of Expression \n ngenes: ", num_genes, sep=""), col=plot_colors)
    }else{
      boxplot((norm_counts), las=2, main=paste("Normalized log2 Values of Expression \n ngenes: ", num_genes, sep=""), col=plot_colors)
    }
    dev.off()

    #plot the normalizaed pca of samples
    if(doLog){
      plot_pca(counts_table, paste(RESULTS_FOLDER, "raw_pca.pdf", sep=""), plot_colors, plot_labels, plot_main)
    }else{
      plot_pca(counts_table, paste(RESULTS_FOLDER, "raw_pca.pdf", sep=""), plot_colors, plot_labels, plot_main, takeLog=F)
    }

    
    #make plots to show how the data looks after
    #plot dispersion
    #if(plot_disp){
    #  pdf(paste(RESULTS_FOLDER, "normalized_disp_cancer.pdf", sep=""))
    #  plotDispNorm(cds, name="disp_cancer", int_norm_counts=norm_counts)
    #  dev.off()

    #  pdf(paste(RESULTS_FOLDER, "normalized_disp_normal.pdf", sep=""))
    #  plotDispNorm(cds, name="disp_normal", int_norm_counts=norm_counts)
    #  dev.off()
    #}

 
    #plot the normalizaed pca of samples
    if(doLog){
      plot_pca(norm_counts, paste(RESULTS_FOLDER, "normalized_pca.pdf", sep=""), plot_colors, plot_labels, plot_main)
    }else{
      plot_pca(norm_counts, paste(RESULTS_FOLDER, "normalized_pca.pdf", sep=""), plot_colors, plot_labels, plot_main, takeLog=F)
    }
    hierarchical_clust(norm_counts, paste(RESULTS_FOLDER, "normalized_samp_clust.pdf", sep=""), plot_main)


}

do_deseq_normalization_and_plot_onfile <-function(countTable_file, RESULTS_FOLDER, conditions, plot_colors, filter_type, filter_param, plot_labels, plot_main){
  ### countTable_file is the location of the raw counts file
  ### RESULTS_FOLDER is where you want all the plots to go
  ### conditions are the things you want to test like Treated vs Untreated = c(rep("treated", 10), rep("untreated", 10))
      ### they should be associated with a column in the countTable_file

  #read table
  countsTable <- read.delim(countTable_file,header=TRUE)
  rownames(countsTable) <- countsTable$gene
  countsTable <- countsTable[,-1]

  if(filter_type == "filter_threshhold_counts"){
    countsTableFiltered = filter_threshhold_counts(countsTable, filter_param)
  }else{
    countsTableFiltered = filter_lowest_counts(countsTable, filter_param)

  }

  conds=factor(conditions)
  cds <- newCountDataSet(countsTableFiltered, conds)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions( cds)

  norm_counts=counts(cds, normalized=TRUE)

  do_normalization_check_plotting(cds, norm_counts, countsTableFiltered, RESULTS_FOLDER, plot_colors, plot_labels, plot_main)
  
  return(list(cds, norm_counts, countsTableFiltered))

}

do_deseq_normalization_and_plot_ontable <-function(countsTable, RESULTS_FOLDER, conditions, plot_colors, filter_type, filter_param, plot_labels, plot_main, fitType="parametric"){
  ### countsTable is the raw counts
  ### RESULTS_FOLDER is where you want all the plots to go
  ### conditions are the things you want to test like Treated vs Untreated = c(rep("treated", 10), rep("untreated", 10))
      ### they should be associated with a column in the countTable_file


  #remove any genes where there is less than 1 counts across in any column
  if(filter_type == "filter_threshhold_counts"){
    countsTableFiltered = filter_threshhold_counts(countsTable, filter_param)
  }else{
    countsTableFiltered = filter_lowest_counts(countsTable, filter_param)

  }
  if(is.null(dim(conditions))){
    conds=factor(conditions)
  }else{
    conds = conditions
  }
  cds <- newCountDataSet(countsTableFiltered, conds)
  cds <- estimateSizeFactors(cds)
  if(ncol(countsTableFiltered) > 200){
    cds <- estimateDispersions(cds, fitType = "local")
  }else{
    cds <- estimateDispersions(cds, fitType=fitType)
  }

  norm_counts=counts(cds, normalized=TRUE)

  #just in case the plot labels were not set, use the conditions as the plot labels 
  if(length(plot_labels) == 0){
    plot_labels = conditions
  }

  do_normalization_check_plotting(cds, norm_counts, countsTableFiltered, RESULTS_FOLDER, plot_colors, plot_labels, plot_main)

  return(list(cds, norm_counts, countsTableFiltered))

}


do_deseq_spike_normalization_and_plot_ontable <-function(countsTable, RESULTS_FOLDER, conditions, plot_colors, filter_type, filter_param, plot_labels, plot_main, fitType="parametric", spike_idx){
  ### countsTable is the raw counts
  ### RESULTS_FOLDER is where you want all the plots to go
  ### conditions are the things you want to test like Treated vs Untreated = c(rep("treated", 10), rep("untreated", 10))
      ### they should be associated with a column in the countTable_file


  #remove any genes where there is less than 1 counts across in any column
  if(filter_type == "filter_threshhold_counts"){
    countsTableFiltered = filter_threshhold_counts(countsTable, filter_param)
  }else{
    countsTableFiltered = filter_lowest_counts(countsTable, filter_param)

  }
  if(is.null(dim(conditions))){
    conds=factor(conditions)
  }else{
    conds = conditions
  }
  cds <- newCountDataSet(countsTableFiltered, conds)
  cds <- estimateSizeFactors(cds)
  if(ncol(countsTableFiltered) > 200){
    cds <- estimateDispersions(cds, fitType = "local")
  }else{
    cds <- estimateDispersions(cds, fitType=fitType)
  }

  norm_counts=counts(cds, normalized=TRUE)

  #just in case the plot labels were not set, use the conditions as the plot labels 
  if(length(plot_labels) == 0){
    plot_labels = conditions
  }

  do_normalization_check_plotting(cds, norm_counts, countsTableFiltered, RESULTS_FOLDER, plot_colors, plot_labels, plot_main)

  return(list(cds, norm_counts, countsTableFiltered))

}
