#!/usr/bin/env python
# coding: utf-8

# In[1]:


import import_ipynb
import Scanpy_functions_v03262021 as sc_pipe
import scvelo as scv
scv.logging.print_version()
import warnings
import scirpy as ir
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
import bbknn
import logging
from sklearn.mixture import GaussianMixture
from scipy.stats     import norm
import glob
import os
import hvplot.pandas
import docx
from docx import Document
from docx.shared import Inches
from docx.shared import Pt
from scipy import sparse
import scanpy.external as sce
import holoviews as hv
import panel as pn
import bokeh
from bokeh.resources import INLINE
import scanorama
import gseapy


# In[2]:


# define sample metadata. Usually read from a file.
exclude_genes = []

samples = {
'GA23428':{},
'GA23429':{},
'GA18593':{},
'GA18594':{},
'GA11556':{},
'GA11557':{}
}

In_path ='/user/ifrec/liuyuchen/scRNASeq_DATA/scRNASeq_BD_Shimizu'
out_path ='/user/ifrec/liuyuchen/scRNASeq_DATA/scRNASeq_BD_Shimizu'
multi_sample = True
sc.settings.verbosity = 3
sc.settings.figdir = out_path+'/figures'


# In[3]:


sc.set_figure_params(dpi=300, dpi_save=300)


# In[4]:


ref=sc.read_mtx('/user/ifrec/liuyuchen/GSE185711/SCP738expression616423b1771a5b086db80551umis.mtx')


# In[5]:


ref


# In[6]:


meta= pd.read_csv('/user/ifrec/liuyuchen/GSE185711/SCP738metadatametaData.txt',sep='\t')


# In[7]:


meta=meta.drop(0, axis=0)


# In[8]:


meta.index= meta['NAME']


# In[9]:


meta


# In[10]:


ref=ref.T


# In[11]:


ref.obs=meta


# In[12]:


ref.obs


# In[13]:


genes =pd.read_csv('/user/ifrec/liuyuchen/GSE185711/SCP738expression616423b1771a5b086db80551umis_genes.tsv',header=None)


# In[14]:


ref.var_names=genes[0].tolist()


# In[15]:


fibro= sc.read_h5ad(out_path+'/Fibroblasts_sub_0602.h5ad')


# In[16]:


adata=fibro.copy()


# In[17]:


adata_ref =ref


# In[18]:


var_names = adata_ref.var_names.intersection(adata.var_names)
adata_ref = adata_ref[:, var_names]
adata = adata[:, var_names]


# In[19]:


sc.pp.pca(adata_ref)
sc.pp.neighbors(adata_ref)
sc.tl.umap(adata_ref)


# In[20]:


sc.tl.leiden(adata_ref)


# In[21]:


sc.pl.umap(adata_ref, color='leiden')


# In[22]:


sc.pl.umap(adata_ref, color='cell_type_integrated')


# In[23]:


sc.pl.umap(adata_ref, color='disease__ontology_label')


# In[24]:


sc.pl.umap(adata_ref, color='organ__ontology_label')


# In[25]:


adata_ref


# In[26]:


sc.pp.neighbors(adata)


# In[27]:


sc.tl.ingest(adata, adata_ref, obs=['disease__ontology_label','organ__ontology_label'])


# In[28]:


sc.pl.umap(adata,color=['disease__ontology_label','organ__ontology_label'])


# In[29]:


adata_concat = adata_ref.concatenate(adata, batch_categories=['ref', 'manuscript'])
sc.pl.umap(adata_concat, color='batch')
sc.pl.umap(adata_concat,color= 'disease__ontology_label')
sc.pl.umap(adata_concat,color='organ__ontology_label')


# In[30]:


fibro.obs=adata.obs


# In[31]:


sc.pl.umap(fibro,color='disease__ontology_label')


# In[32]:


sc.pl.umap(fibro,color='organ__ontology_label')


# In[33]:


#adata_ref.write('/user/ifrec/liuyuchen/GSE185711/GSE185711.h5ad')


# In[34]:


sc.pl.umap(adata, color ='leiden')


# In[35]:


sc.pl.umap(adata, color ='leiden',legend_loc='on data')


# In[36]:


adata


# In[10]:


sc.tl.rank_genes_groups(adata, groupby='leiden',use_raw=False)
rank_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
rank_pavlues = pd.DataFrame(adata.uns['rank_genes_groups']['pvals'])
rank_pavlues.rename(columns=lambda x: x+'_pvals', inplace=True)
rank_fold = pd.DataFrame(adata.uns['rank_genes_groups']['logfoldchanges'])
rank_fold.rename(columns=lambda x: x+'_logfoldchanges', inplace=True)
pd.concat([rank_genes, rank_pavlues,rank_fold], axis=1, join='inner').to_csv(out_path+'/Fibroblasts_By_subclusters_0602.csv')



# In[11]:


x = pd.concat([rank_genes, rank_pavlues,rank_fold], axis=1, join='inner')


# In[18]:


sc.get.rank_genes_groups_df(adata,group=None,log2fc_min=5).to_csv('sub_log2fold_over_5.csv')


# In[9]:


sc.pl.rank_genes_groups_heatmap(adata,n_genes=100,standard_scale='var')


# In[ ]:





# In[10]:


clusters=adata.obs.leiden.unique().tolist()
for c in clusters:
    q= rank_genes[c].head(200).tolist()
    sc.queries.enrich(q, org="hsapiens",gprofiler_kwargs={'no_evidences':False}).to_csv(out_path+'/cluster_'+c+'_fibroblast_enrichment.csv')


# In[11]:


gos=pd.read_csv('GO_term_summary_20230213_192237.csv')


# In[12]:


markers = gos['Symbol'].str.upper().tolist()


# In[13]:


set(rank_genes['6'].head(100).tolist()).intersection(set(markers))


# In[14]:


for c in ['1','7','6','10','11']:
    q= list(set(rank_genes[c].head(100).tolist()).intersection(set(markers)))
    sc.pl.stacked_violin(adata,q,groupby='leiden',title='cluster '+str(c))


# In[15]:


adata.obs['Leiden']=adata.obs['leiden']


# In[16]:


adata.write(out_path+'/Fibroblasts_sub_0602.h5ad')


# In[17]:


from pyBCS import scanpy2bcs
scanpy2bcs.format_data(out_path+'/Fibroblasts_sub_0602.h5ad', out_path+'/Fibroblasts_sub_0602.bcs',
                        input_format="h5ad", graph_based="leiden")


# In[15]:


top_15_Keratinocytes = list(set(sum(rank_genes.head(15).values.tolist(),[])))


# In[ ]:


sc.pl.stacked_violin(sub,top_15_Keratinocytes,groupby='Conditions')


# In[ ]:


for g in top_15_Keratinocytes:
    sc.pl.violin(sub,g,groupby='Conditions')


# In[ ]:


sub = adata[adata.obs['Cell_typing_Panglao']=='Endothelial cells'].copy()


# In[ ]:


sc.tl.rank_genes_groups(sub, groupby='Conditions',reference='Skin',use_raw=False)


# In[ ]:


rank_genes = pd.DataFrame(sub.uns['rank_genes_groups']['names'])
rank_pavlues = pd.DataFrame(sub.uns['rank_genes_groups']['pvals'])
rank_pavlues.rename(columns=lambda x: x+'_pvals', inplace=True)
rank_fold = pd.DataFrame(sub.uns['rank_genes_groups']['logfoldchanges'])
rank_fold.rename(columns=lambda x: x+'_logfoldchanges', inplace=True)
pd.concat([rank_genes, rank_pavlues,rank_fold], axis=1, join='inner').to_csv(out_path+'/Endothelial_Bycondition_1220.csv')


# In[ ]:


sc.pl.rank_genes_groups(sub, groupby='Conditions',reference='Skin',use_raw=False)


# In[ ]:


x =  pd.DataFrame(sub.uns['rank_genes_groups']['names'])
y =  pd.DataFrame(sub.uns['rank_genes_groups']['scores'])
z =  pd.DataFrame(sub.uns['rank_genes_groups']['logfoldchanges'])
gene_table = pd.concat([x,y,z], axis=1, join='inner').head(10)
gene_table.columns =['genes','scores','logfoldchanges']
gene_table['scores'] = gene_table['scores'].astype(int)
gene_table.index =gene_table['genes']
gene_table.plot(kind="bar")
sc.pl.rank_genes_groups(sub, groupby='Conditions',reference='Skin',use_raw=False)


# In[ ]:


top_15_Endothelial = list(set(sum(rank_genes.head(15).values.tolist(),[])))


# In[ ]:


sc.pl.stacked_violin(sub,top_15_Endothelial,groupby='Conditions')


# In[ ]:


for g in top_15_Endothelial:
    sc.pl.violin(sub,g,groupby='Conditions')


# In[ ]:


sub = adata[adata.obs['Cell_typing_Panglao']=='Fibroblasts'].copy()
sc.tl.rank_genes_groups(sub, groupby='Conditions',reference='Skin',use_raw=False)
rank_genes = pd.DataFrame(sub.uns['rank_genes_groups']['names'])
rank_pavlues = pd.DataFrame(sub.uns['rank_genes_groups']['pvals'])
rank_pavlues.rename(columns=lambda x: x+'_pvals', inplace=True)
rank_fold = pd.DataFrame(sub.uns['rank_genes_groups']['logfoldchanges'])
rank_fold.rename(columns=lambda x: x+'_logfoldchanges', inplace=True)
pd.concat([rank_genes, rank_pavlues,rank_fold], axis=1, join='inner').to_csv(out_path+'/Fibroblasts_Bycondition_1220.csv')
x =  pd.DataFrame(sub.uns['rank_genes_groups']['names'])
y =  pd.DataFrame(sub.uns['rank_genes_groups']['scores'])
z =  pd.DataFrame(sub.uns['rank_genes_groups']['logfoldchanges'])
gene_table = pd.concat([x,y,z], axis=1, join='inner')

gene_table.columns =['genes','scores','logfoldchanges']
gene_table['scores'] = gene_table['scores'].astype(int)
gene_table.index =gene_table['genes']
gene_table=gene_table[gene_table['genes'].isin(['IGF2','SFRP2','INHBA','PRSS23','EGFL6','INHBB','INHA'])]
gene_table.plot(kind="bar")
sc.pl.rank_genes_groups(sub, groupby='Conditions',reference='Skin',use_raw=False)


# In[ ]:


gene_table[gene_table['genes'].isin(['IGF2','SFRP2','INHBA','PRSS23','EGFL6'])]


# In[ ]:


top_15_Fibroblasts = list(set(sum(rank_genes.head(15).values.tolist(),[])))


# In[ ]:


sc.pl.stacked_violin(sub,top_15_Fibroblasts,groupby='Conditions')


# In[ ]:


for g in top_15_Fibroblasts:
    sc.pl.violin(sub,g,groupby='Conditions')


# In[ ]:


sc.pl.umap(adata, color ='INHBA', use_raw=False,cmap='OrRd')


# In[ ]:


sc.pl.umap(adata, color =['IGF2','SFRP2','INHBA','PRSS23','EGFL6','INHBB'], use_raw=False,cmap='OrRd')


# In[ ]:


sc.pl.umap(adata, color ='INHA',cmap='OrRd')


# In[ ]:




