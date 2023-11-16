import scanpy as sc
import anndata
import warnings

#Define wrappers to run scANVI and scVI (which is required to pre-train scANVI)

def scvi(adata, batch, dims, hvg=None, return_model=False, max_epochs=None):
    """scVI wrapper function

    Adapted from scib package version 1.1.3, based on scvi-tools version >=0.16.0 (available through `conda <https://docs.scvi-tools.org/en/stable/installation.html>`_)

    .. note::
        scVI expects only non-normalized (count) data on highly variable genes!

    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param dims: number of dimensions to use for integration
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix as well as an embedding representation of the
        corrected data
    """
    import numpy as np
    try:
        from scvi.model import SCVI
    except ModuleNotFoundError as e:
        raise OptionalDependencyNotInstalled(e)

  #  scib.utils.check_sanity(adata, batch, hvg)

    # Check for counts data layer
    if "counts" not in adata.layers:
        raise TypeError(
            "Adata does not contain a `counts` layer in `adata.layers[`counts`]`"
        )
        
    n_latent = dims #Defaults from SCVI github tutorials scanpy_pbmc3k and harmonization is 30
    
    # Defaults from SCVI github tutorials scanpy_pbmc3k and harmonization
    n_hidden = 128
    n_layers = 2

    # copying to not return values added to adata during setup_anndata
    net_adata = adata.copy()
    if hvg is not None:
        net_adata = adata[:, hvg].copy()
    SCVI.setup_anndata(net_adata, layer="counts", batch_key=batch)

    vae = SCVI(
        net_adata,
        gene_likelihood="nb",
        n_layers=n_layers,
        n_latent=n_latent,
        n_hidden=n_hidden,
    )
    train_kwargs = {"train_size": 1.0}
    if max_epochs is not None:
        train_kwargs["max_epochs"] = max_epochs
    vae.train(**train_kwargs)
    adata.obsm["X_emb"] = vae.get_latent_representation()

    if not return_model:
        return adata
    else:
        return vae
      

def scanvi(adata, batch, dims, labels, hvg=None, max_epochs=None):
    """scANVI wrapper function

    Adapted from scib package version 1.1.3, based on scvi-tools version >=0.16.0 (available through `conda <https://docs.scvi-tools.org/en/stable/installation.html>`_)

    .. note::
        Use non-normalized (count) data for scANVI!

    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param dims: number of dimensions to use for integration
    :param labels: label key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix as well as an embedding representation of the
        corrected data
    """
    import numpy as np
    try:
        from scvi.model import SCANVI
    except ModuleNotFoundError as e:
        raise OptionalDependencyNotInstalled(e)


    if max_epochs is None:
        n_epochs_scVI = int(np.min([round((20000 / adata.n_obs) * 400), 400]))
        n_epochs_scANVI = int(np.min([np.max([2, round(n_epochs_scVI / 3.0)]), 2]))
    else:
        n_epochs_scVI = max_epochs
        n_epochs_scANVI = max_epochs

    vae = scvi(adata, batch, dims, hvg, return_model=True,max_epochs=n_epochs_scVI)

    # STEP 1: RUN scVI to initialize scANVI
    scanvae = SCANVI.from_scvi_model(
        scvi_model=vae,
        labels_key=labels,
        unlabeled_category="unknown"
    )
    # STEP 2: RUN scANVI
    scanvae.train(max_epochs=n_epochs_scANVI, train_size=1.0)
    adata.obsm["X_emb"] = scanvae.get_latent_representation()

    return adata

def calc_hvg(adata, nhvg, batch):

	# remove HVG if already precomputed
	if 'highly_variable' in adata.var:
		del adata.var['highly_variable']

	h = sc.pp.highly_variable_genes(
		adata,
		flavor="seurat_v3",
		n_top_genes=nhvg,
		layer="counts",
		batch_key=batch,
		subset=False,
		inplace=False
	)
	h = h[h.highly_variable==True]
	return list(h.index)


if __name__ == '__main__':
	file = "cache/merged.h5ad"
	batch = "SampleLabel"
	labels = "scGate_multi"
	dims = 50
	nhvg = 800
	outPath = "out/scANVI_integrated_CD8.h5ad"

	adata = anndata.read_h5ad(file)
	adata.layers["counts"] = adata.raw.X.copy()

#	hvg = list(adata.var.index)
	hvg = calc_hvg(adata, nhvg, batch)
	print(hvg[0:10])

	integrated = scanvi(adata, batch, dims, labels, hvg)

	sc.pp.neighbors(integrated, use_rep="X_emb")
	sc.tl.leiden(integrated)
	sc.tl.umap(integrated)
	del integrated.layers

	sc.write(outPath, integrated)



