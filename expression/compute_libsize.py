import sys
import scipy as sp
import h5py
import os
import re
import pdb

if len(sys.argv) < 3:
    print >> sys.stderr, 'Usage: %s <expression.tsv> <outfile.tsv>' % sys.argv[0]
    sys.exit(1)
fname = sys.argv[1]
fn_out = sys.argv[2]
### get list of protein coding genes
coding = sp.loadtxt('/cluster/work/grlab/projects/TCGA/PanCancer/annotation/gencode.v19.annotation.coding_gene_list.tab', dtype='str', delimiter='\t')
coding = coding[:, [0, 4]]
### filter for autosomes
k_idx = sp.where(~sp.in1d(coding[:, 0], sp.array(['chrMT', 'chrX', 'chrY'])))[0]
coding = coding[k_idx, :]
coding = coding[:, 1]

print 'loading expression data from ' + fname

data = sp.loadtxt(fname, delimiter = '\t', dtype = 'string')
expression = data[:,1].astype('float')
genes = data[:,0]
strains = [fname]


k_idx = sp.where(sp.in1d(genes, coding))[0]
genes = genes[k_idx]
print genes.shape[0]
expression = expression[k_idx]
libsize_uq = sp.percentile(expression, 75)
libsize_tc = expression.sum(axis=0)

s_idx = sp.argsort(libsize_uq)[::-1]
out = open(fn_out,'w')#open(re.sub('.tsv$', '', fname) + '.libsize.tsv', 'w')
print >> out, 'sample\tlibsize_75percent\tlibsize_total_count'
print >> out, '\t'.join([strains[0], str(libsize_uq), str(libsize_tc)])
out.close()
