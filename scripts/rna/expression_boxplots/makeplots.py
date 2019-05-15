import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import seaborn as sbn
import scipy as sp
import h5py
import sys
import os
import pdb

relgenes = ["ENSG00000167601","ENSG00000157764","ENSG00000012048","ENSG00000139618",
            "ENSG00000110092","ENSG00000118971","ENSG00000112576",
            "ENSG00000135446","ENSG00000105810","ENSG00000147889",
            "ENSG00000138798","ENSG00000146648","ENSG00000078098",
            "ENSG00000138685","ENSG00000077782","ENSG00000066468","ENSG00000068078",
            "ENSG00000088256","ENSG00000156052","ENSG00000019991","ENSG00000157404",
            "ENSG00000049130","ENSG00000169032","ENSG00000100030",
            "ENSG00000105976","ENSG00000120215","ENSG00000198793",
            "ENSG00000196712","ENSG00000213281","ENSG00000134853",
            "ENSG00000121879","ENSG00000154229","ENSG00000163932",
            "ENSG00000171862","ENSG00000160307","ENSG00000185664",
            "ENSG00000100146","ENSG00000125398","ENSG00000197122","ENSG00000107165"]
relgenes_plain = ["AXL","BRAF","BRCA1","BRCA2",
                  "CCND1","CCND2","CCND3",
                  "CDK4","CDK6","CDKN2A",
                  "EGF","EGFR","FAP",
                  "FGF2","FGFR1","FGFR2","FGFR3",
                  "GNA11","GNAQ","HGF","KIT",
                  "KITLG","MAP2K1","MAPK1",
                  "MET","MLANA","MTOR",
                  "NF1","NRAS","PDGFRA",
                  "PIK3CA","PRKCA","PRKCD",
                  "PTEN","S100B","PMEL",
                  "SOX10","SOX9","SRC",
                  "TYRP1"]
                  
dict_gene = dict(zip(relgenes, relgenes_plain))

fn_tcga = "/cluster/work/grlab/projects/TCGA/PanCancer/rerun2018_hdf5/expression_counts.hdf5"
fn_skcmids = "/cluster/work/grlab/projects/TCGA/PanCanAtlas/variants_mc3/pancan.merged.v0.2.6.PUBLIC.maf.donors_SKCM"
if __name__ == "__main__":
    fn_input = sys.argv[1]
    skcmids = sp.loadtxt(fn_skcmids, delimiter = '\t', dtype = 'string')[:,0]


    samplename = fn_input.split('/')[-1].split('__')[0]
    print samplename
    IN = h5py.File(fn_tcga, 'r')
    gids = sp.array([x.split('.')[0] for x in IN['gids'][:]])
    midx = sp.in1d(gids,relgenes)

    counts = IN['counts'][:]
    libsize =  sp.percentile(counts,75,axis=0)
    libsize = libsize.astype('float')
    libsize = counts/libsize

    gids = gids[midx]
    sids = IN['sids'][:]
    sids = sp.array([x.split('.')[0].rsplit('-',4)[0] for x in sids])
    smidx = sp.in1d(sids, skcmids)
    counts = counts[:,smidx]
    IN.close()

    data = sp.loadtxt(fn_input, delimiter = '\t', dtype = 'string')
    gids_ref = sp.array([x.split('.')[0] for x in data[:,0]])
    counts_ref = data[:,1].astype('float')
    libsize = sp.percentile(counts_ref, 75)
    counts_ref = counts_ref / libsize

    midx = sp.in1d(gids_ref, relgenes)
    gids_ref = gids_ref[midx]
    counts_ref = counts_ref[midx]
    assert (sp.all(gids == gids_ref)),"not matching genes"
    
    for i,g in enumerate(gids):
        fig = plt.figure(figsize=(4,10))
        ax = fig.add_subplot(111)
        sbn.boxplot(counts[i,:],orient='v',showfliers=False, ax = ax)
        ax.plot(0,counts_ref[i],marker='*',markersize=20,color="red")
        # print g, 
        # print dict_gene[g]
        # if relgenes_plain[i] == "CCND1":
        #     pdb.set_trace()
        
        plt.savefig('%s__boxplot_%s_former_%s.png' % (samplename,dict_gene[g],relgenes_plain[i]))
        plt.clf()
