
import scanpy as sc
import pandas as pd

##
sc.settings.verbosity = 1  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor="white")
##
# this is an earlier version of the dataset from the pbmc3k tutorial
adata_ref = sc.datasets.pbmc3k_processed()
adata = sc.datasets.pbmc68k_reduced()
## define the same var
var_names = adata_ref.var_names.intersection(adata.var_names)
adata_ref = adata_ref[:, var_names]
adata = adata[:, var_names]
##
sc.pp.pca(adata_ref)
sc.pp.neighbors(adata_ref)
sc.tl.umap(adata_ref)
# sc.pl.umap(adata_ref, color="louvain")
## Mapping
sc.tl.ingest(adata, adata_ref, obs="louvain")
adata.uns["louvain_colors"] = adata_ref.uns["louvain_colors"]  # fix colors
sc.pl.umap(adata, color=["louvain", "bulk_labels"], wspace=0.5)
##
adata_concat = adata_ref.concatenate(adata, batch_categories=["ref", "new"])  # concat both
adata_concat.obs.louvain = adata_concat.obs.louvain.astype("category")
# fix category ordering
adata_concat.obs.louvain.cat.reorder_categories(
    adata_ref.obs.louvain.cat.categories
)
# fix category colors
adata_concat.uns["louvain_colors"] = adata_ref.uns["louvain_colors"]
sc.pl.umap(adata_concat, color=["batch", "louvain"])

## working on a smaple dataset to see batch effect
# note that this collection of batches is already intersected on the genes
adata_all = sc.read(
    "data/pancreas.h5ad",
    backup_url="https://www.dropbox.com/s/qj1jlm9w10wmt0u/pancreas.h5ad?dl=1",
)
counts = adata_all.obs.celltype.value_counts()
# filter data
minority_classes = counts.index[-5:].tolist()  # get the minority classes
adata_all = adata_all[~adata_all.obs.celltype.isin(minority_classes)]  # actually subset
adata_all.obs.celltype.cat.reorder_categories(  # reorder according to abundance
    counts.index[:-5].tolist()
)
# plot embedding
sc.pp.pca(adata_all)
sc.pp.neighbors(adata_all)
sc.tl.umap(adata_all)
sc.pl.umap(
    adata_all, color=["batch", "celltype"], palette=sc.pl.palettes.vega_20_scanpy
)
## Removing batch via BBKNN
sc.external.pp.bbknn(adata_all, batch_key="batch")
sc.tl.umap(adata_all)
sc.pl.umap(adata_all, color=["batch", "celltype"])

## Back to scanpy
adata_ref = adata_all[adata_all.obs.batch == "0"]
sc.pp.pca(adata_ref)
sc.pp.neighbors(adata_ref)
sc.tl.umap(adata_ref)
sc.pl.umap(adata_ref, color="celltype")

adatas = [adata_all[adata_all.obs.batch == i].copy() for i in ["1", "2", "3"]]
sc.settings.verbosity = 2  # a bit more logging
for iadata, adata in enumerate(adatas):
    print(f"... integrating batch {iadata+1}")
    adata.obs["celltype_orig"] = adata.obs.celltype  # save the original cell type
    sc.tl.ingest(adata, adata_ref, obs="celltype")


adata_concat = adata_ref.concatenate(adatas)
adata_concat.obs.celltype = adata_concat.obs.celltype.astype("category")
# fix category ordering
adata_concat.obs.celltype.cat.reorder_categories(
    adata_ref.obs.celltype.cat.categories)
# fix category coloring
adata_concat.uns["celltype_colors"] = adata_ref.uns["celltype_colors"]

sc.pl.umap(adata_concat, color=["batch", "celltype"])
## check consitancy
adata_query = adata_concat[adata_concat.obs.batch.isin(["1", "2", "3"])]
sc.pl.umap(adata_query, color=["batch", "celltype", "celltype_orig"], wspace=0.4)