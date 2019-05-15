#source("https://bioconductor.org/biocLite.R")

#library("biomaRt")
library(VennDiagram)

make_venn_diagram <-function(list_of_elements, filename, title, col_vec){
  venn.diagram(list_of_elements, filename, height = 3000, width = 6000, resolution = 500, units = "px", compression = "lzw", na = "stop", main = title, fill=col_vec)

}


getGeneName <- function(ens_names, merge_on_col, EG2EZ_name){
  geneNames = merge(EG2EZ_name, ens_names, by=merge_on_col)
  return(geneNames)
}

remove_period <- function(GeneID_vector){
  myInterestingGenes = as.character(GeneID_vector)
  myInterestingGenes=gsub("[.][0-9]*", "", myInterestingGenes)
  return(myInterestingGenes)
}

translate_NAMETYPE_TO_ENS <-function(genes_to_trans, org="hsapiens_gene_ensembl", name_type='hgnc_symbol'){
  # select mart and data set
  bm <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  bm <- useDataset(org, mart=bm)

  # Get ensembl gene ids and hgnc_symbol
  # this is our translation table
  EG2HGNC_name <- getBM(mart=bm, attributes=c('ensembl_gene_id', name_type))
  EG2HGNC_name <- na.omit(EG2HGNC_name)

  colnames(genes_to_trans)=c(name_type)

  genes_to_trans = merge(EG2HGNC_name, genes_to_trans, by=name_type, all.y=TRUE)
  genes_to_trans=unique(genes_to_trans)

  genes_to_trans$hgnc_symbol[is.na(genes_to_trans$hgnc_symbol)] <- as.character(genes_to_trans$ensembl_gene_id[is.na(genes_to_trans$hgnc_symbol)])
  genes_to_trans$hgnc_symbol[genes_to_trans$hgnc_symbol==""] <- as.character(genes_to_trans$ensembl_gene_id[genes_to_trans$hgnc_symbol==""])

  return(genes_to_trans)
}

translate_using_tables <-function(genes_to_trans, table_to_use, col_interest, col_merge){

  #read in reference file
  genename_table = read.table(table_to_use, header=T)
  genename_table = genename_table[, c(col_interest, col_merge)]

  # format the input file
  genes_to_trans_id = remove_period(genes_to_trans$id)
  temp=colnames(genes_to_trans)
  genes_to_trans=cbind(genes_to_trans, genes_to_trans_id)
  colnames(genes_to_trans)=c(temp, col_merge)

  # merge
  genes_to_trans = merge(genename_table, genes_to_trans, by=col_merge, all.y=TRUE)
  genes_to_trans = unique(genes_to_trans)

  #genes_to_trans[is.na(genes_to_trans$mgi_symbol), col_interest] <- as.character(genes_to_trans[is.na(genes_to_trans$mgi_symbol), "id"])
  #genes_to_trans[genes_to_trans$mgi_symbol=="", col_interest] <- as.character(genes_to_trans[genes_to_trans$mgi_symbol=="", "id"])
  #genes_to_trans[genes_to_trans$mgi_symbol=="NA", col_interest] <- as.character(genes_to_trans[genes_to_trans$mgi_symbol=="NA", "id"])

  return(genes_to_trans)


}


translate_ENS_to_MGI <-function(genes_to_trans, org="mmusculus_gene_ensembl", name_type='mgi_symbol'){
  # select mart and data set
  bm <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  bm <- useDataset(org, mart=bm)

  # Get ensembl gene ids and hgnc_symbol
  # this is our translation table
  EG2mgi_name <- getBM(mart=bm, attributes=c('ensembl_gene_id', name_type))
  EG2mgi_name <- na.omit(EG2mgi_name)

  ##DE GENES
  #get the DE genes
  #remove the .xxx at the end of the gene names
  genes_to_trans_id = remove_period(genes_to_trans$id)

  temp=colnames(genes_to_trans)
  genes_to_trans=cbind(genes_to_trans, genes_to_trans_id)
  colnames(genes_to_trans)=c(temp, "ensembl_gene_id")

  genes_to_trans = merge(EG2mgi_name, genes_to_trans, by="ensembl_gene_id", all.y=TRUE)
  genes_to_trans=unique(genes_to_trans)

  genes_to_trans$mgi_symbol[is.na(genes_to_trans$mgi_symbol)] <- as.character(genes_to_trans$ensembl_gene_id[is.na(genes_to_trans$mgi_symbol)])
  genes_to_trans$mgi_symbol[genes_to_trans$mgi_symbol==""] <- as.character(genes_to_trans$ensembl_gene_id[genes_to_trans$mgi_symbol==""])

  return(genes_to_trans)
}


translate_ENS_to_HGNC <-function(genes_to_trans, org="hsapiens_gene_ensembl", name_type='hgnc_symbol'){

  # remove NA from result
  genes_to_trans = na.omit(genes_to_trans)

  # select mart and data set
  bm <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  bm <- useDataset(org, mart=bm)

  # Get ensembl gene ids and hgnc_symbol
  # this is our translation table
  EG2HGNC_name <- getBM(mart=bm, attributes=c('ensembl_gene_id', name_type))
  EG2HGNC_name <- na.omit(EG2HGNC_name)

  ##DE GENES
  #get the DE genes
  #remove the .xxx at the end of the gene names
  genes_to_trans_id = remove_period(genes_to_trans$id)

  temp=colnames(genes_to_trans)
  genes_to_trans=cbind(genes_to_trans, genes_to_trans_id)
  colnames(genes_to_trans)=c(temp, "ensembl_gene_id")

  genes_to_trans = merge(data.table(EG2HGNC_name), data.table(genes_to_trans), by="ensembl_gene_id", all.y=TRUE)
  genes_to_trans = data.frame(genes_to_trans)
  genes_to_trans=unique(genes_to_trans)

  genes_to_trans$hgnc_symbol[is.na(genes_to_trans$hgnc_symbol)] <- as.character(genes_to_trans$ensembl_gene_id[is.na(genes_to_trans$hgnc_symbol)])
  genes_to_trans$hgnc_symbol[genes_to_trans$hgnc_symbol==""] <- as.character(genes_to_trans$ensembl_gene_id[genes_to_trans$hgnc_symbol==""])

  return(genes_to_trans)
}


translate_REFSEQ_to_ENS <-function(genes_to_trans){
  # select mart and data set
  bm <- useMart("ensembl")
  bm <- useDataset("hsapiens_gene_ensembl", mart=bm)

  # Get ensembl gene ids and hgnc_symbol
  # this is our translation table
  EG2HGNC_name <- getBM(mart=bm, attributes=c('ensembl_gene_id', 'hgnc_symbol', 'refseq_mrna'))
  EG2HGNC_name <- na.omit(EG2HGNC_name)

  ##DE GENES
  #get the DE genes
  #remove the .xxx at the end of the gene names
  #genes_to_trans_id = remove_period(genes_to_trans$id)
  genes_to_trans = other_countsTable
  colnames(genes_to_trans)[colnames(genes_to_trans) == "RefSeqID"] = "refseq_mrna"

  genes_to_trans = merge(EG2HGNC_name, genes_to_trans, by="refseq_mrna", all.y=TRUE)
  genes_to_trans=unique(genes_to_trans)

  genes_to_trans$hgnc_symbol[is.na(genes_to_trans$hgnc_symbol)] <- as.character(genes_to_trans$ensembl_gene_id[is.na(genes_to_trans$hgnc_symbol)])

  return(genes_to_trans)
}

write_out_genes_mouse <-function(filename_prefix, de_table, pval, first, second){
  ### filename_prefix where the files will be written to, and a prefix for the filename
  ### de_table is the table of DE genes
  ### pval is the adjusted pval cutoff

  filename=paste(filename_prefix, first, second, "MA.pdf", sep="_")
  pdf(filename)
  plotMA(de_table)
  dev.off()

  de_table = de_table[,c(1, 2, 3, 4, 6, 8)]

  tot_res_all = de_table
  tot_res = de_table[de_table$padj < pval,]

  filename=paste(filename_prefix, first, second, "-.txt", sep="_")
  neg_res = tot_res[tot_res[,c("log2FoldChange")]<=0,]
  neg_res = translate_ENS_to_MGI(neg_res)
  write.table( neg_res[ order(neg_res[,c("padj")]), ], file=filename, quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\r\n")

  filename=paste(filename_prefix, first, second, "+.txt", sep="_")
  pos_res = tot_res[tot_res[,c("log2FoldChange")]>0,]
  pos_res = translate_ENS_to_MGI(pos_res)
  if(! is.null(dim(pos_res))){
    write.table( pos_res[ order(pos_res[,c("padj")]), ], file=filename, quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n")
  }

  filename=paste(filename_prefix, first, second, "all.txt", sep="_")
  tot_res = translate_ENS_to_MGI(tot_res)
  write.table( tot_res[ order(tot_res[,c("padj")]), ], file=filename, quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\r\n")

  filename=paste(filename_prefix, first, second, "ranked.txt", sep="_")
  translated_genes = translate_ENS_to_MGI(tot_res_all)
  tot_res_all = translated_genes
  write.table( tot_res_all[ order(tot_res_all[,c("padj")]), c("mgi_symbol")], file=filename, quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\r\n")


}



#a function to write out DE genes nicely.
write_out_genes <-function(filename_prefix, de_table, pval, first, second, translateHumanGenes=F){
  ### filename_prefix where the files will be written to, and a prefix for the filename
  ### de_table is the table of DE genes
  ### pval is the adjusted pval cutoff

  filename=paste(filename_prefix, first, second, "MA.pdf", sep="_")
  pdf(filename)
  plotMA(de_table)
  dev.off()

  de_table = de_table[,c(1, 2, 3, 4, 6, 8)]

  tot_res_all = de_table
  tot_res = de_table[de_table$padj < pval,]

  filename=paste(filename_prefix, first, second, "-.txt", sep="_")
  neg_res = tot_res[tot_res[,c("log2FoldChange")]<=0,]
  neg_res = translate_ENS_to_HGNC(neg_res)
  write.table( neg_res[ order(neg_res[,c("padj")]), ], file=filename, quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\r\n")

  filename=paste(filename_prefix, first, second, "+.txt", sep="_")
  pos_res = tot_res[tot_res[,c("log2FoldChange")]>0,]
  pos_res = translate_ENS_to_HGNC(pos_res)
  if(! is.null(dim(pos_res))){
    write.table( pos_res[ order(pos_res[,c("padj")]), ], file=filename, quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n")
  }

  filename=paste(filename_prefix, first, second, "all.txt", sep="_")
  tot_res = translate_ENS_to_HGNC(tot_res)
  write.table( tot_res[ order(tot_res[,c("padj")]), ], file=filename, quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\r\n")

  filename=paste(filename_prefix, first, second, "ranked.txt", sep="_")
  if(translateHumanGenes){
    translated_genes = translate_ENS_to_HGNC(tot_res_all)
    tot_res_all = translated_genes
    write.table( tot_res_all[ order(tot_res_all[,c("padj")]), c("hgnc_symbol")], file=filename, quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\r\n")

  }else{
    write.table( tot_res_all[ order(tot_res_all[,c("padj")]), c("id")], file=filename, quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\r\n")
  }

}

#method for finding the DE genes between 2 conditions
call_DE_single_cond <-function(cds, first, second, pval){
  ### the DeSeq object
  ### first is the first condition
  ### second is the second condition
  ### pval is the adjusted pval cutoff

  res = nbinomTest( cds, first, second)
  resSig = na.omit(res[ res$padj < pval, ])
  cat("number of genes DE ", dim(resSig))
  return(res)

}

call_DE_many_cond <-function(cds, all_pairs_conds, outfile, pval, translateHumanGenes=F, translateMouseGenes=F, testInteraction=F){
  ### cds is the DeSeq object
  ### all_pairs_conds is a matrix of condition comparisons. the first row will be compared to the second row
      ### ex = cond1 cond2 cond3
      ###      cond2 cond3 cond4
      ### this will make 3 comparisons, cond1 vs cond2, cond2 vs cond3, cond3 vs cond4

  ### outfile is where you want the gene lists written to
  ### pval is the adjusted pval cutoff

  ### RETURN TYPE: a vector of named lists.
  ###              Each named list has first, second, res.
  ###              First and second are the comparisons, res is the DeSeq result

  total_results=list()

  for(i in 1:dim(all_pairs_conds)[2]){
    first = all_pairs_conds[1,i]
    second = all_pairs_conds[2,i]

    print(i)

    if(!testInteraction){
    	res=call_DE_single_cond(cds, first, second, pval)
    }else{
    	res=call_interactionDE_single_cond(cds, first, second, pval)

    }
    output = list(first=c(as.matrix(first)), second=c(as.matrix(second)), res=res)
    total_results = c(total_results, list(output))
    print(first)
    print(second)
    cat(i, " out of ", dim(all_pairs_conds)[2], "\n")
    if(translateMouseGenes){
      write_out_genes_mouse(outfile, res, pval, first, second)
    }else{
      write_out_genes(outfile, res, pval, first, second, translateHumanGenes)

    }

  }

  return(total_results)
}
