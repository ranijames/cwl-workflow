###  basic DCA methods


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


# DCA takes in the prediction label, expression matr, genes of interest
# and returns the projected points, the projection vector, AUROC within and AUROC pan

plot_tsne <- function(plot_matr, color_labels, outdir, outstr){

  #PC_tot = Rtsne(scale(log10(plot_matr+1)))
  PC_tot = try(Rtsne((log10(plot_matr+1))), TRUE)
  if(inherits(PC_tot, "try-error")){
    return()
  }


  #PC_tot[is.nan(PC_tot)] <- 0

  PC_plot = data.frame(t1=PC_tot$Y[,1], t2=PC_tot$Y[,2])
  colnames(PC_plot) = c("t1", "t2")
  PC_plot$types = color_labels

  outfile = paste(outdir, outstr, "_tsne.pdf", sep="")
  pdf(outfile)
    gg = ggplot(PC_plot, aes(t1, t2, color=types)) +
          geom_point()
    print(gg)
  dev.off()

}

plot_boxplot_expr <- function(plot_df, genes_interest, outdir, outstr){

  plot_df_melt = melt(plot_df[,c(genes_interest, "study", "is_normal")], id=c("study", "is_normal"))
  outfile = paste(outdir, outstr, "_expr_boxplot.pdf", sep="")
  pdf(outfile)
    gg = ggplot(plot_df_melt, aes(x=study, y=log10(value+1), fill=is_normal)) +
          geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(gg)

  dev.off()

}

plot_pca_dca <- function(plot_matr, color_labels, outdir, outstr){

  PC_tot = prcomp(log10(plot_matr+1))
  PC_plot = data.frame(PC_tot$x[,1:3])
  colnames(PC_plot) = c("PC1", "PC2", "PC3")
  PC_plot$types = color_labels

  eigs <- PC_tot$sdev^2
  var_explained = eigs[1] / sum(eigs)

  outfile = paste(outdir, outstr, "_pca.pdf", sep="")
  pdf(outfile)
    gg = ggplot(PC_plot, aes(PC1, PC2, color=types)) +
          geom_point() + labs(x = eigs[1] / sum(eigs), y=eigs[2] / sum(eigs))
    print(gg)
  dev.off()

}


get_pred_perf <- function(dca_proj_compare, col_interest){

  pred_dca = prediction(dca_proj_compare[,col_interest], dca_proj_compare$Y == 0)
  perf_dca_roc = performance(pred_dca, "tpr", "fpr")
  auc_dca = performance(pred_dca, "auc")
  if(unlist(auc_dca@y.values) < 0.5){
    pred_dca = prediction(dca_proj_compare[,col_interest], dca_proj_compare$Y == 1)
    perf_dca_roc = performance(pred_dca, "tpr", "fpr")
  }
  perf_dca_pr = performance(pred_dca, "prec", "rec")

  require(MESS)
  x <- perf_dca_pr@x.values[[1]] # Recall values
  y <- perf_dca_pr@y.values[[1]]

  aucpr = try(auc(x,y, type = 'spline'), TRUE)
  if(inherits(aucpr, "try-error")){
    aucpr = NA
  }

  return(list(perf_dca_roc, perf_dca_pr, auc_dca, aucpr))
}

get_proj_error <- function(full_matr, genes_interest, proj_vector, test_idx, target_values){

  # run test
  test_matr = hif_full_data[test_idx,]
  test_matr$Y = target_values[test_idx]

  colnames_avail = genes_interest[which(genes_interest %in% colnames(test_matr))]
  #dca_data = as.matrix(scale(log10(test_matr[,colnames_avail]+1)))
  dca_data = as.matrix((log10(test_matr[,colnames_avail]+1)))

  dca_proj = dca_data %*% proj_vector
  dca_res_match = data.frame(dca_res=dca_proj)
  dca_res_match$row_id = row.names(dca_res_match)
  dca_res_match$Y = test_matr$Y

  dca_res_match = dca_res_match[order(dca_res_match$dca_res),]

  prediction_res = get_pred_perf(dca_res_match, "dca_res")
  auc_dca = prediction_res[[3]]
  auc_roc = unlist(auc_dca@y.values)
  if(auc_roc < 0.5){
    auc_roc = 1-auc_roc
  }
  auc_pr = prediction_res[[4]]

  return(list(auc_roc, auc_pr, perf_roc=prediction_res[[1]], perf_pr=prediction_res[[2]]))

}

do_dca_proj <- function(full_matr, genes_interest, train_idx, test_idx, target_values, useScale=T){

  full_matr$Y = target_values

  train_matr = full_matr[train_idx,]
  test_matr = full_matr[test_idx,]

  colnames_avail = unique(genes_interest[which(genes_interest %in% colnames(full_matr))])

  dca_train_matr = train_matr[,colnames_avail]
  chunks = rep(1, nrow(train_matr))
  chunks[train_matr$Y!=0] = 2
  neglinks = matrix(c(0, 1, 1, 0), 2, 2)
  #dca_res = dca(data=scale(log10(dca_train_matr+1)), chunks=chunks, neglinks=neglinks)

  if(useScale){
    dca_res = dca(data=t(scale(log10(t(dca_train_matr+1)))), chunks=chunks, neglinks=neglinks)
  }else{
    dca_res = dca(data=(log10(dca_train_matr+1)), chunks=chunks, neglinks=neglinks)
  }
  proj_vector = t(as.matrix(dca_res$DCA))

  # run test
  dca_test_matr = test_matr[,colnames_avail]+1
  #dca_data = as.matrix(scale(log10(dca_test_matr)))
  if(useScale){
    dca_data = as.matrix(t(scale(log10(t(dca_test_matr+1)))))
  }else{
    dca_data = as.matrix(log10(dca_test_matr+1))
  }
  dca_data[is.nan(dca_data)] <- 0
  dca_proj = dca_data %*% proj_vector
  dca_res_match = data.frame(dca_res=dca_proj)
  dca_res_match$row_id = row.names(dca_res_match)
  dca_res_match$Y = test_matr$Y

  dca_res_match = dca_res_match[order(dca_res_match$dca_res),]

  #prediction_res = get_pred_perf(dca_res_match, "dca_res")
  #auc_dca = prediction_res[[3]]
  #auc_roc = unlist(auc_dca@y.values)
  #if(auc_roc < 0.5){
    #auc_roc = 1-auc_roc
  #}
  #auc_pr = prediction_res[[4]]


  return(list(proj_vector, dca_res_match))

}


get_cv_auc <- function(full_matr, genes_interest, idx_to_use, target_values){

  colnames_avail = genes_interest[which(genes_interest %in% colnames(full_matr))]

  full_matr$Y = target_values
  full_matr = full_matr[idx_to_use,]
  #full_matr = full_matr[full_matr$is_normal == FALSE,]

  idx_neg = which(full_matr$Y == 0)
  idx_pos = which(full_matr$Y == 1)

  # do 80/20 split
  idx_length_neg = floor(length(idx_neg)/10)
  idx_length_pos = floor(length(idx_pos)/10)
  avg_cv = c()
  avg_pr = c()
  for(curr_split in 1:10){

    start_idx_neg = (curr_split-1)*idx_length_neg+1
    end_idx_neg = curr_split*idx_length_neg
    test_idx_neg = idx_neg[start_idx_neg:end_idx_neg]

    if(length(idx_pos) < 100){
      start_idx_pos = 1
      end_idx_pos = floor(length(idx_pos)/2)
      test_idx_pos = idx_pos[start_idx_pos:end_idx_pos]

    }else{
      start_idx_pos = (curr_split-1)*idx_length_pos+1
      end_idx_pos = curr_split*idx_length_pos
      test_idx_pos = idx_pos[start_idx_pos:end_idx_pos]
    }

    test_idx = c(test_idx_neg, test_idx_pos)
    train_idx = which(!1:nrow(full_matr) %in% test_idx)

    curr_res = do_dca_proj(full_matr, colnames_avail, train_idx, test_idx, full_matr$Y)
    avg_cv = c(avg_cv, curr_res[[3]])
    avg_pr = c(avg_pr, curr_res[[4]])

  }

  print(avg_cv)
  perf_roc = curr_res[[5]]
  perf_pr = curr_res[[6]]
  return(list(mean(avg_cv), mean(avg_pr), perf_roc, perf_pr))

}
