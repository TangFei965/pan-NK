import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import bbknn
import os
import re
import scvelo as scv

scv.logging.print_version()

sc.settings.verbosity = 3  
sc.logging.print_versions()
sns.despine()
sns.set_style("whitegrid")
sc.settings.set_figure_params(dpi=80, dpi_save=300, fontsize = 14, frameon = True,figsize=(3,3))

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo')

cd16 = sc.read_h5ad("cd16.h5ad")
adata_toplot = sc.read_h5ad("NK_expr_withVelo.h5ad")

adata_toplot = adata_toplot[adata_toplot.obs['cellID'].isin(list(cd16.obs.index))]
adata_toplot.obs =  cd16[list(adata_toplot.obs['cellID'])].obs
adata_toplot.obsm['X_umap'] =  cd16[list(adata_toplot.obs['cellID'])].obsm['X_umap']

scv.pp.filter_genes(adata_toplot, min_shared_counts=20)
sc.pp.normalize_total(adata_toplot, target_sum=1e4) 
sc.pp.log1p(adata_toplot)
sc.pp.highly_variable_genes(adata_toplot, min_mean=0.0125, max_mean=3, min_disp=0.5)



adata_toplot1 = adata_toplot[:,adata_toplot.var.highly_variable]

scv.pp.moments(adata_toplot1, n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(adata_toplot1, n_jobs = 36)
scv.tl.velocity(adata_toplot1, group_by = 'celltype', mode = 'dynamical')

scv.tl.velocity_graph(adata_toplot1) 
scv.tl.terminal_states(adata_toplot1)

adata_toplot1.uns["celltype_colors"] =['#EBD57C',
'#F47A79',
'#ED4437',
'#FE9E37',
'#B24646',
'#96873B',
'#B49D99',
'#B37557',
]

scv.pl.velocity_embedding_stream(adata_toplot1, basis='umap', color = ['celltype'],density = 1.5, alpha=0.4,linewidth=0.8)


