import scanpy as sc
import numpy as np
import pandas as pd
from collections import Counter
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, facecolor='white')

cd16 = sc.read_h5ad('cd16.h5ad')

##infer the differentiation trajectory of all tumor-infiltrating CD56dimCD16hi NK cells following the workflow in Scanpy
cd16t = cd16[cd16.obs['meta_tissue']=='Tumor']
##built diffusion map  
sc.tl.diffmap(cd16t,n_comps=10)

#selected the CD56brightCD16hi NK cells as the root because they located at the start of the directed streamline inferred by RNA velocity
cd16t.uns['iroot'] = np.flatnonzero(cd16t.obs['celltype']  == 'CD56brightCD16hi')[0]
sc.tl.dpt(cd16t)


##fit the expression profile of each gene to the pseudotime
def make_gene_df(gene_list,adata):
    df = pd.DataFrame()
    for each_gene in gene_list:
        df[each_gene]=adata.raw[:,each_gene].X.todense().A1
    df['Dpt_pseudotime']=adata.obs['dpt_pseudotime'].tolist()
    return df

geneset1 = ["GZMA","GZMB","PRF1"]
geneset2 = ['HSP90AA1','HSPA1A','DNAJB1']
geneset3 = ["KIR3DL1","KIR3DL2"]
geneset4 = ['NR4A1','NR4A2','NR4A3']

sns.set_style("white")


df = pd.concat([make_gene_df(i,cd16t) for i in geneset1])
sns.lmplot(scatter=False,data=df,x='Dpt_pseudotime',y='Gene Expression',
           hue="Gene",
           order=2,        
          )
plt.savefig('cd16_cyto.pdf')


df = pd.concat([make_gene_df(i,cd16t) for i in geneset2])
sns.lmplot(scatter=False,data=df,x='Dpt_pseudotime',y='Gene Expression',
           hue="Gene",
           order=2,        
          )
plt.savefig('cd16_dys.pdf')


df = pd.concat([make_gene_df(i,cd16t) for i in geneset3])
sns.lmplot(scatter=False,data=df,x='Dpt_pseudotime',y='Gene Expression',
           hue="Gene",
           order=2,        
          )
plt.savefig('cd16_inhi.pdf')


df = pd.concat([make_gene_df(i,cd16t) for i in geneset4])
sns.lmplot(scatter=False,data=df,x='Dpt_pseudotime',y='Gene Expression',
           hue="Gene",
           order=2,        
          )
plt.savefig('cd16_NR4A.pdf')


