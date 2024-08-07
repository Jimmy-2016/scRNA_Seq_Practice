import numpy as np
import matplotlib.pyplot as pl
import scanpy as sc

##
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = "./write/paul15.h5ad"
# low dpi (dots per inch) yields small inline figures
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor="white")
adata = sc.datasets.paul15()
# this is not required and results will be comparable without it
adata.X = adata.X.astype("float64")
## pre_processing
sc.pp.recipe_zheng17(adata)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata, color="paul15_clusters", legend_loc="on data")
## diffuion map for prunning the graph
sc.tl.diffmap(adata)
sc.pp.neighbors(adata, n_neighbors=10, use_rep="X_diffmap")
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata, color="paul15_clusters", legend_loc="on data")
##  clustring and paga
sc.tl.louvain(adata, resolution=1.0)
sc.tl.paga(adata, groups="louvain")
sc.pl.paga(adata, color=["louvain", "Hba-a2", "Elane", "Irf8"])
sc.pl.paga(adata, color=["louvain", "Itga2b", "Prss34", "Cma1"])
##
adata.obs["louvain"].cat.categories
adata.obs["louvain_anno"] = adata.obs["louvain"]
adata.obs["louvain_anno"].cat.categories = [
    *["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"],
    *["10/Ery", "11", "12", "13", "14", "15"],
    *["16/Stem", "17", "18"],
    *[
        "19/Neu",
        "20/Mk",
        "21",
    ],
    *["22/Baso", "23", "24/Mo"],
]


zeileis_colors = np.array(sc.pl.palettes.zeileis_28)
new_colors = np.array(adata.uns["louvain_anno_colors"])
new_colors[[16]] = zeileis_colors[[12]]  # Stem colors / green
new_colors[[10, 17, 5, 3, 15, 6, 18, 13, 7, 12]] = zeileis_colors[  # Ery colors / red
    [5, 5, 5, 5, 11, 11, 10, 9, 21, 21]
]
new_colors[[20, 8]] = zeileis_colors[[17, 16]]  # Mk early Ery colors / yellow
new_colors[[4, 0]] = zeileis_colors[[2, 8]]  # lymph progenitors / grey
new_colors[[22]] = zeileis_colors[[18]]  # Baso / turquoise
new_colors[[19, 14, 2]] = zeileis_colors[[6, 6, 6]]  # Neu / light blue
new_colors[[24, 9, 1, 11]] = zeileis_colors[[0, 0, 0, 0]]  # Mo / dark blue
new_colors[[21, 23]] = zeileis_colors[[25, 25]]  # outliers / grey
adata.uns["louvain_anno_colors"] = new_colors

adata.uns["iroot"] = np.flatnonzero(adata.obs["louvain_anno"] == "16/Stem")[0]
sc.tl.dpt(adata)
gene_names = [
    *["Gata2", "Gata1", "Klf1", "Epor", "Hba-a2"],  # erythroid
    *["Elane", "Cebpe", "Gfi1"],  # neutrophil
    *["Irf8", "Csf1r", "Ctsg"],  # monocyte
]
adata_raw = sc.datasets.paul15()
sc.pp.log1p(adata_raw)
sc.pp.scale(adata_raw)
adata.raw = adata_raw
sc.pl.draw_graph(adata, color=["louvain_anno", "dpt_pseudotime"], legend_loc="on data")