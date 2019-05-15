# Paths for HIF project

# The following lines are added to script to make sure the R script is executable as a CommandLineTool
# Therefore it can be used for CWL tool development

args = commandArgs(trailingOnly=TRUE)
RESULTS_FOLDER     =  paste(as.character(args[1]), "/", sep='')
cancertype         =  as.character(args[2])
sample_id          =  as.character(args[3])
count_expression   =  as.character(args[4])
lib_size           =  as.character(args[5])
labkey_prefix     =  as.character(args[6])

HDF5_LIBRARY_SIZE_GTEX_FILE = "/cluster/work/grlab/projects/tumor_profiler/data/TCGA/hypoxia_score_resources/counts_gtex_2015-03-29.libsize.tsv"
HDF5_LIBRARY_SIZE_FILE      = "/cluster/work/grlab/projects/tumor_profiler/data/TCGA/hypoxia_score_resources/expression_counts.whitelisted.libsize.tsv"
HIF_GENE_SIG_FILE           = "/cluster/work/grlab/projects/tumor_profiler/data/TCGA/hypoxia_score_resources/hif_sig_expr.tsv"
HIF_GENE_SIG_GTEX_FILE      = "/cluster/work/grlab/projects/tumor_profiler/data/TCGA/hypoxia_score_resources/hif_sig_gtex_expr.tsv"
GTEX_META_DATA_FILE         = "/cluster/work/grlab/projects/tumor_profiler/data/TCGA/hypoxia_score_resources/SraRunTable.txt"
PSI_TCGA_FILE               = "/cluster/work/grlab/projects/tumor_profiler/data/TCGA/hypoxia_score_resources/event_positions_v.2_tcga_PSI_andre_new.txt"
EXON_SKIP_926_PATH          = "/cluster/work/grlab/projects/tumor_profiler/data/TCGA/hypoxia_score_resources/exon_skip_926_count.txt"
GENDER_FILE                 = "/cluster/work/grlab/projects/tumor_profiler/data/TCGA/hypoxia_score_resources/genderData.tsv"
PSI_TCGA_FILE               = "/cluster/work/grlab/projects/tumor_profiler/data/TCGA/hypoxia_score_resources/event_positions_v.2_tcga_PSI_andre_new.txt"
PSI_GTEX_FILE               = "/cluster/work/grlab/projects/tumor_profiler/data/TCGA/hypoxia_score_resources/event_positions_v.2_gtex_PSI.txt"
survival_file               = "/cluster/work/grlab/projects/tumor_profiler/data/TCGA/hypoxia_score_resources/clinical_var.tsv"


tcga_meta_data_file  = "/cluster/work/grlab/projects/TCGA/PanCancer/metadata/2015-08-07/full_metadata.overlap_RNA_WXS.RNA.tsv"
PSI_ICGC_NORMAL_FILE = "/cluster/work/grlab/projects/Krek/results/spladder_normal/sample_events_normal_tcga_PSI.tsv"
PSI_GTEX_NORMAL_FILE = "/cluster/work/grlab/projects/Krek/results/spladder_normal/sample_events_normal_gtex_PSI.tsv"
normal_sample_dir = "/cluster/work/grlab/projects/Krek/results/spladder_normal/"

ENSEMBL_TO_TRANS_FILE = "/cluster/work/grlab/projects/Krek/annotation/GRCh38_p10_ens_to_trans.txt"
ENSEMBL_TO_ENTREZ_FILE = "/cluster/work/grlab/projects/Krek/annotation/GRCh38_p12_ens_trans_entrez.tsv"

HDF5_EXPRESSION_FILE = "/cluster/work/grlab/projects/TCGA/PanCancer/rerun2018_hdf5/expression_counts.whitelisted.hdf5"
HIF_GENE_RAND_FILE = "/cluster/work/grlab/projects/Krek/results/tcga_events/hif_rand_expr.tsv"
TCGA_MUTSIG_FILE = "/cluster/work/grlab/projects/TCGA/PanCanAtlas/variants_mc3/mutload.tsv.gz"

HIF_GENE_RAND_GTEX_FILE = "/cluster/work/grlab/projects/Krek/results/tcga_events/hif_rand_gtex_expr.tsv"

HDF5_EXPRESSION_GTEX_FILE = "/cluster/work/grlab/projects/ICGC/counts_gtex_2015-03-29.hdf5"
TCGA_CNV_FILE = "/cluster/work/grlab/projects/TCGA/PanCancer/hdf5_10Kset/cnvdata.hdf5"


DCARS_FILE = "/cluster/work/grlab/projects/Krek/results/dcars/full_dcars_res.tsv"

EXPRESSION_DE_FOLDER = "/cluster/work/grlab/projects/Krek/data/expression_counts/"

KEGG_PATHWAYS_DIR = "/cluster/home/nrd44/checkouts/kegg_pathways/"

KEGG_PATHWAY_SPIA =  "/cluster/work/grlab/projects/Krek/spia/"

RAS_ACTIVATION_FILE = "/cluster/work/grlab/projects/Krek/annotation/ras_activation_status.csv"
RESISTANCE_FILE = "/cluster/work/grlab/projects/Krek/annotation/E-MTAB-3645-A-AFFY-141-normalized-expressions.tsv"

OV_CISPLATIN_FILE = "/cluster/work/grlab/projects/Krek/annotation/GSE98230_counts.tsv"

GENE_WEIGHTS_FILE =  "/cluster/work/grlab/projects/Krek/results/tcga_events_strat/patient_proj/gene_weights_hypoxia.tsv"

TCGA_GTEX_HIF_SCORES = "/cluster/work/grlab/projects/Krek/results/tcga_events_strat/paper_figures_tables/frozen_files/sample_scores.tsv"

HIF_SCORE_RESULTS_DIR = "/cluster/work/grlab/projects/Krek/results/tcga_events_strat/paper_figures_tables/"
GENE_WEIGHTS_FILE = paste(HIF_SCORE_RESULTS_DIR, "gene_weights.tsv", sep="")
GTEX_HIF_METADATA_FILE = paste(HIF_SCORE_RESULTS_DIR, "gtex_metadata.tsv", sep="")
TCGA_HIF_METADATA_FILE = paste(HIF_SCORE_RESULTS_DIR, "tcga_metadata.tsv", sep="")

HIF_EXPR_COUNTS_DIR = "/cluster/work/grlab/projects/Krek/data/expression_counts/"

VHL_MUTANT_FILE = "/cluster/work/grlab/projects/Krek/annotation/GSE74958_rcc_read_count_mat.txt"
TCGA_VHL_MUTANT_FILE = "/cluster/work/grlab/projects/Krek/annotation/tcga_ccRCC_vlh_mut.tsv"

HRD_SCORE_ID_FILE = "/cluster/work/grlab/projects/Krek/annotation/DDR_score_ids.tsv"
HRD_SAMPLE_ID_FILE = "/cluster/work/grlab/projects/Krek/annotation/DDR_sample_ids.tsv"
HRD_SCORE_FILE = "/cluster/work/grlab/projects/Krek/annotation/DDR_scores.tsv"

CLONAL_INFO_FILE = "/cluster/work/grlab/projects/Krek/annotation/cirello_tumor_clones_tcga.csv"

OV_SIMPLE_FILE = "/cluster/work/grlab/projects/Krek/annotation/GSE52695_Ovarian_RefSeq_RPKM_Values.csv"
BRCA_RADIO_FILE = "/cluster/work/grlab/projects/Krek/annotation/GSE120798_log2_quantile_normalized_counts.txt"
COLOR_TISSUE_FILE = "/cluster/work/grlab/projects/Krek/results/tcga_events_strat/paper_figures_tables/tissue_colors.tsv"

BRCA_RECUR_FILE = "/cluster/work/grlab/projects/Krek/annotation/GSE119937_concat.tsv"
BRCA_RECUR_METADATA_FILE = "/cluster/work/grlab/projects/Krek/annotation/sample_id_recurr.txt"
#tumor_total_necrosis_percent
#specific_tumor_total_necrosis_percent

GBM_RESISTANT_FILE = "/cluster/work/grlab/projects/Krek/annotation/GSE79671_CountMatrix.txt"
GBM_RESISTANT_METADATA_FILE = "/cluster/work/grlab/projects/Krek/annotation/GSE79671_metadata.tsv"
GBM_DRUG_FILE = "/cluster/work/grlab/projects/Krek/annotation/gbm_drug.csv"

is.nan.data.frame <- function(x){ do.call(cbind, lapply(x, is.nan)) }

HIF_GENES = c("ENSG00000148926", "ENSG00000109107", "ENSG00000176171",
                  "ENSG00000104765", "ENSG00000074410", "ENSG00000107159",
                  "ENSG00000130635", "ENSG00000047457", "ENSG00000168209",
                  "ENSG00000129521", "ENSG00000111674", "ENSG00000104812",
                  "ENSG00000159399", "ENSG00000100292", "ENSG00000134333",
                  "ENSG00000113083", "ENSG00000123384", "ENSG00000185499",
                  "ENSG00000119950", "ENSG00000104419", "ENSG00000185633",
                  "ENSG00000124785", "ENSG00000152256", "ENSG00000114268",
                  "ENSG00000204531", "ENSG00000119938", "ENSG00000139832",
                  "ENSG00000141526", "ENSG00000117394", "ENSG00000103257",
                  "ENSG00000113739", "ENSG00000265972", "ENSG00000112715",
                  "ENSG00000186918", "ENSG00000117289")

house_keeping_genes = c("ENSG00000111640", "ENSG00000075624", "ENSG00000165704",
                        "ENSG00000166710", "ENSG00000137818", "ENSG00000112592")
