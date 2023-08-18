import pandas as pd
import numpy as np
import scanpy as sc
from collections import Counter
import os
import anndata as ann
import bbknn

def check_data(test):
    sc.pp.normalize_total(test, target_sum=1e4)
    sc.pp.log1p(test)
    sc.pp.highly_variable_genes(test, min_mean=0.0125, max_mean=3, min_disp=0.5)
    test.raw = test                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
    test = test[:, test.var.highly_variable]

    sc.pp.neighbors(test, n_neighbors=10, n_pcs=40)
    sc.tl.umap(test)
    sc.tl.leiden(test)
    sc.pl.umap(test, color=['leiden', 'NCAM1','NKG7','FCGR3A',"KLRF1"], color_map="Reds",ncols=5 )
    return test

#batch correction
adata = sc.read_h5ad("combine.h5ad")
adata.obs['batch'] = adata.obs['patients'] + "_" + adata.obs['datasets']

patient_name =list(adata.obs.groupby("batch").count()[list(adata.obs.groupby("batch").count().iloc[:,0]>5)].index)
adata =adata[[i in patient_name for i in adata.obs['batch']]]

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5,batch_key="cancer_types")
adata.raw = adata
adata = adata[:,adata.var.highly_variable]
print(adata)
sc.tl.pca(adata,n_comps=40)

bbknn.bbknn(adata,batch_key='batch')

sc.tl.umap(adata)
sc.tl.leiden(adata)

sc.tl.rank_genes_groups(adata, 'leiden', method='logreg',use_raw=False)
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names']}).head(5)

# Excluding other cell types
sc.pl.umap(adata, color=[ 'NCAM1','NKG7','FCGR3A',"KLRF1","LYZ",'CST3','C1QA',
                        "EPCAM","KRT18","KRT19","CD68","CD163","CD79A", "CD19", "MS4A1",'leiden'], color_map="Reds",legend_loc="on data" )
adata = adata[(adata.obs['leiden'] !='10') &
              (adata.obs['leiden']!='9') & (adata.obs['leiden']!='5')]

# Seperation of CD56brightCD16lo and CD56dimCD16hi NK cells based on the expression of NCAM1 and FCGR3A 
sc.pl.umap(adata, color=[ 'NCAM1','FCGR3A'], color_map="Reds",legend_loc="on data" )

adata.obs['majorType'] = "CD56lowCD16high"
adata.obs['majorType'] [(adata.obs['leiden']=="8") | (adata.obs['leiden']=="2")] = "CD56brightCD16lo"

adata = adata.raw.to_adata()

CD16low = adata.obs[adata.obs['majorType']== "CD56brightCD16lo"]
CD16high = adata.obs[adata.obs['majorType']!= "CD56brightCD16lo"]


##Clustering for CD56brightCD16lo NK cells
datalist =list(CD16low.obs.groupby("batch").count()[(CD16low.obs.groupby("batch").count()>5)['sampleID']].index)

CD16low = CD16low[[i in datalist for i in CD16low.obs['batch'] ]]
sc.pp.highly_variable_genes(CD16low, min_mean=0.0125, max_mean=3, min_disp=0.5,batch_key="meta_histology")
CD16low.raw = CD16low
CD16low = CD16low[:,CD16low.var.highly_variable]
print(CD16low )

sc.tl.pca(CD16low,n_comps=40)
bbknn.bbknn(CD16low,batch_key='batch')
sc.tl.umap(CD16low)
sc.tl.leiden(CD16low)
sc.pl.umap(CD16low, color=["meta_tissue",'datasets',"meta_histology","leiden"], color_map="Reds",ncols=2 )

sc.tl.rank_genes_groups(CD16low, 'leiden', method='logreg',use_raw=False)
result = CD16low.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names']}).head(5)

##Clustering for CD56dimCD16hi NK cells
datalist = list(CD16high.obs.groupby("batch").count()[(CD16high.obs.groupby("batch").count()>5)['sampleID']].index)
CD16high = CD16high[[i in datalist for i in CD16high.obs['batch'] ]]

sc.pp.highly_variable_genes(CD16high, min_mean=0.0125, max_mean=3, min_disp=0.5,batch_key="meta_histology")
CD16high.raw = CD16high
CD16high = CD16high[:,CD16high.var.highly_variable]

sc.tl.pca(CD16high,n_comps=40)
bbknn.bbknn(CD16high,batch_key='batch')
sc.tl.umap(CD16high)
sc.tl.leiden(CD16high)
sc.pl.umap(CD16high, color=["meta_tissue",'datasets',"meta_histology"], color_map="Reds",ncols=2 )

sc.tl.rank_genes_groups(CD16high, 'leiden', method='logreg',use_raw=False)
result = CD16high.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names']}).head(5)


#Save
CD16low.write("CD56brightCD16lo.h5ad")
CD16high.write("CD56dimCD16hi.h5ad")
