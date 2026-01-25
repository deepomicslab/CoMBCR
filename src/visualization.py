import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import pandas as pd
import numpy as np
import torch
import torch.nn.functional as F
import anndata
import scanpy as sc
sc.settings.verbosity = 0
import os
from typing import Optional, Union, Tuple, List, Literal
from matplotlib.figure import Figure
import warnings
warnings.filterwarnings('ignore')


def create_joint_embedding_adata(
    bcr_emb_path: str,
    gex_emb_path: str,
    metadata: Optional[Union[str, pd.DataFrame]] = None,
    normalize: bool = True,
    n_neighbors: int = 10,
    n_pcs: int = 30,
    compute_umap: bool = True,
    return_separate_layers: bool = True
) -> anndata.AnnData:
    """
    Create an AnnData object containing joint embeddings
    
    Parameters
    ----------
    bcr_emb_path : str
        Path to BCR embedding CSV file
    gex_emb_path : str
        Path to GEX embedding CSV file
    metadata : Optional[Union[str, pd.DataFrame]], default=None
        Metadata, can be a file path or DataFrame
        If provided, index should be barcode
    normalize : bool, default=True
        Whether to apply L2 normalization to joint embeddings
    n_neighbors : int, default=10
        Number of neighbors for computing neighbor graph
    n_pcs : int, default=30
        Number of principal components for computing neighbor graph
    compute_umap : bool, default=True
        Whether to compute UMAP
    return_separate_layers : bool, default=True
        Whether to save separate BCR and GEX embeddings in layers
    
    Returns
    -------
    anndata.AnnData
        AnnData object containing:
        - X: joint embedding
        - layers['bcr']: BCR embedding (if return_separate_layers=True)
        - layers['gex']: GEX embedding (if return_separate_layers=True)
        - obs: metadata (if provided)
        - obsm['X_umap']: UMAP coordinates (if compute_umap=True)
    
    Examples
    --------
    >>> # Basic usage
    >>> adata = create_joint_embedding_adata(
    ...     bcr_emb_path="path/to/bcrembeddings.csv",
    ...     gex_emb_path="path/to/gexembedding.csv"
    ... )
    
    >>> # With metadata
    >>> adata = create_joint_embedding_adata(
    ...     bcr_emb_path="path/to/embedding.csv",
    ...     gex_emb_path="path/to/gexembedding.csv",
    ...     metadata="annotation.csv"
    ... )
    
    >>> # Without UMAP and separate layers
    >>> adata = create_joint_embedding_adata(
    ...     bcr_emb_path="path/to/embedding.csv",
    ...     gex_emb_path="path/to/gexembedding.csv",
    ...     compute_umap=False,
    ...     return_separate_layers=False
    ... )
    """
    
    # Load embeddings
    print("Loading embeddings...")
    bcremb = pd.read_csv(bcr_emb_path, index_col='barcode')
    gexemb = pd.read_csv(gex_emb_path, index_col='barcode')
    
    # Ensure both embeddings have the same barcodes
    assert bcremb.index.equals(gexemb.index), "BCR and GEX embeddings must have the same barcodes"
    
    # Concatenate embeddings
    print("Concatenating embeddings...")
    embeddings = np.concatenate([bcremb.values, gexemb.values], axis=1)
    
    # Normalize
    if normalize:
        print("Normalizing joint embeddings...")
        embeddings = torch.from_numpy(embeddings)
        embeddings = F.normalize(embeddings, dim=-1).numpy()
    
    # Create AnnData object with named index
    print("Creating AnnData object...")
    obs_df = pd.DataFrame(index=bcremb.index.tolist())
    obs_df.index.name = 'barcode' 
    
    adata = anndata.AnnData(
        X=embeddings,
        obs=obs_df
    )
    
    # Store separate embedding layers
    if return_separate_layers:
        print("Storing separate embeddings in layers...")
        adata.obsm['CoMBCR_bcr'] = bcremb.values
        adata.obsm['CoMBCR_gex'] = gexemb.values
    
    # Add metadata
    if metadata is not None:
        print("Adding metadata...")
        if isinstance(metadata, str):
            metadata_df = pd.read_csv(metadata, index_col="barcode")
        else:
            metadata_df = metadata.copy()
        
        # Ensure index matching
        common_barcodes = adata.obs.index.intersection(metadata_df.index)
        if len(common_barcodes) < len(adata.obs):
            print(f"Warning: Only {len(common_barcodes)}/{len(adata.obs)} barcodes found in metadata")
        
        adata.obs = adata.obs.join(metadata_df, how='left')
    
    # Compute neighbor graph
    print(f"Computing neighbor graph for Joint embedding ...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    # Compute UMAP
    if compute_umap:
        print("Computing UMAP for Joint embedding ...")
        sc.tl.umap(adata)
    
    return adata

def create_sub_embedding_adatas(
    adata: anndata.AnnData,
    compute_umap: bool = True,
    copy_obs: bool = True,
    copy_uns: bool = False
) -> Tuple[anndata.AnnData, anndata.AnnData]:
    """
    Create separate AnnData objects for BCR and GEX embeddings from joint embedding AnnData
    
    This function extracts the BCR and GEX embeddings stored in adata.obsm and creates
    independent AnnData objects with the same neighbor graph parameters as the original.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Joint embedding AnnData object created by create_joint_embedding_adata()
        Must contain 'CoMBCR_bcr' and 'CoMBCR_gex' in obsm
    compute_umap : bool, default=True
        Whether to compute UMAP for each sub-embedding
    copy_obs : bool, default=True
        Whether to copy observation annotations (metadata) to sub-embeddings
    copy_uns : bool, default=False
        Whether to copy unstructured annotations to sub-embeddings
    
    Returns
    -------
    Tuple[anndata.AnnData, anndata.AnnData]
        bcr_adata, gex_adata - Two AnnData objects containing BCR and GEX embeddings
        Each will have:
        - X: the respective embedding
        - obs: metadata (if copy_obs=True)
        - uns: unstructured data (if copy_uns=True)
        - obsp['connectivities']: neighbor graph
        - obsp['distances']: distance matrix
        - obsm['X_umap']: UMAP coordinates (if compute_umap=True)
    
    Raises
    ------
    ValueError
        If 'CoMBCR_bcr' or 'CoMBCR_gex' not found in adata.obsm
    KeyError
        If neighbor graph parameters not found in adata.uns
    
    Examples
    --------
    >>> # Create joint embedding first
    >>> adata = create_joint_embedding_adata(
    ...     bcr_emb_path="bcr.csv",
    ...     gex_emb_path="gex.csv",
    ...     metadata="metadata.csv",
    ...     n_neighbors=15
    ... )
    
    >>> # Extract sub-embeddings with same parameters
    >>> bcr_adata, gex_adata = create_sub_embedding_adatas(adata)
    
    >>> # Visualize BCR embedding
    >>> sc.pl.umap(bcr_adata, color='celltypes')
    
    >>> # Visualize GEX embedding
    >>> sc.pl.umap(gex_adata, color='tissue')
    
    >>> # Without UMAP computation
    >>> bcr_adata, gex_adata = create_sub_embedding_adatas(
    ...     adata, 
    ...     compute_umap=False
    ... )
    """
    
    # Check if required embeddings exist
    if 'CoMBCR_bcr' not in adata.obsm:
        raise ValueError(
            "'CoMBCR_bcr' not found in adata.obsm. "
            "Make sure to create the joint embedding with return_separate_layers=True"
        )
    if 'CoMBCR_gex' not in adata.obsm:
        raise ValueError(
            "'CoMBCR_gex' not found in adata.obsm. "
            "Make sure to create the joint embedding with return_separate_layers=True"
        )
    
    # Extract neighbor graph parameters from original adata
    try:
        n_neighbors = adata.uns['neighbors']['params']['n_neighbors']
    except KeyError:
        print("Warning: Could not find neighbor parameters in adata.uns, using default n_neighbors=10")
        n_neighbors = 10
    
    # Create BCR AnnData
    print("Creating BCR AnnData...")
    bcr_adata = anndata.AnnData(
        X=adata.obsm['CoMBCR_bcr'].copy()
    )
    bcr_adata.obs_names = adata.obs_names
    bcr_adata.obs.index.name = 'barcode'
    
    if copy_obs:
        bcr_adata.obs = adata.obs.copy()
    
    if copy_uns:
        bcr_adata.uns = adata.uns.copy()
    
    # Store embedding info
    bcr_adata.uns['embedding_type'] = 'BCR'
    if 'CoMBCR' in adata.uns:
        bcr_adata.uns['embedding_dim'] = adata.uns['CoMBCR']['bcr_dim']
    
    # Create GEX AnnData
    print("Creating GEX AnnData...")
    gex_adata = anndata.AnnData(
        X=adata.obsm['CoMBCR_gex'].copy()
    )
    gex_adata.obs_names = adata.obs_names
    gex_adata.obs.index.name = 'barcode'
    
    if copy_obs:
        gex_adata.obs = adata.obs.copy()
    
    if copy_uns:
        gex_adata.uns = adata.uns.copy()
    
    # Store embedding info
    gex_adata.uns['embedding_type'] = 'GEX'
    if 'CoMBCR' in adata.uns:
        gex_adata.uns['embedding_dim'] = adata.uns['CoMBCR']['gex_dim']
    
    # Compute neighbor graphs with same parameters
    print(f"Computing neighbor graph for BCR embedding ...")
    sc.pp.neighbors(bcr_adata, n_neighbors=n_neighbors, use_rep='X')
    
    print(f"Computing neighbor graph for GEX embedding ...")
    sc.pp.neighbors(gex_adata, n_neighbors=n_neighbors, use_rep='X')
    
    # Compute UMAP if requested
    if compute_umap:
        print("Computing UMAP for BCR embedding...")
        sc.tl.umap(bcr_adata)
        
        print("Computing UMAP for GEX embedding...")
        sc.tl.umap(gex_adata)
    
    print("Done!")
    return bcr_adata, gex_adata


def plot_training_loss(
    log_path: str,
    mode: Literal['earlystopping', 'save_epoch'] = 'earlystopping',
    save_path: Optional[str] = None
) -> Figure:
    """
    Plot training loss curves with analysis
    
    Parameters
    ----------
    log_path : str
        Path to the training log file (e.g., 'example_outdir/CoMBCR.pth.log')
        The log file should be tab-separated with columns: Epoch, loss, loss_cmc, 
        loss_p2p, loss_b2b
    mode : {'earlystopping', 'save_epoch'}, default='earlystopping'
        Training mode to determine which epoch's checkpoint was saved:
        - 'earlystopping': Model saved at epoch with minimum total loss
        - 'save_epoch': Model saved at the last epoch
    save_path : Optional[str], default=None
        Path to save the figure. If None, figure will be displayed but not saved.
        Recommended format: .png, .pdf, or .svg
        
    Returns:
    --------
    fig : matplotlib.figure.Figure
        The generated figure object
    """
    
    # Read the log file
    if not os.path.exists(log_path):
        raise FileNotFoundError(f"Log file not found: {log_path}")
    
    logfile = pd.read_csv(log_path, sep='\t')
    
    # Convert to numpy arrays
    epochs = logfile.Epoch.values
    loss = np.array(logfile.loss.values)
    loss_cmc = np.array(logfile.loss_cmc.values)
    loss_p2p = np.array(logfile.loss_p2p.values)
    loss_b2b = np.array(logfile.loss_b2b.values)
    
    # Determine the saved checkpoint epoch based on mode
    if mode == 'earlystopping':
        # Find epoch with minimum total loss
        saved_epoch_idx = np.argmin(loss)
        saved_epoch = epochs[saved_epoch_idx]
        mode_label = "Early Stop"
    elif mode == 'save_epoch':
        # Last epoch
        saved_epoch_idx = len(epochs) - 1
        saved_epoch = epochs[saved_epoch_idx]
        mode_label = "Specified Epoch"
    else:
        raise ValueError(f"Invalid mode: {mode}. Must be 'earlystopping' or 'save_epoch'")
    
    # Get loss values at saved epoch
    saved_loss = loss[saved_epoch_idx]
    saved_cmc = loss_cmc[saved_epoch_idx]
    saved_p2p = loss_p2p[saved_epoch_idx]
    saved_b2b = loss_b2b[saved_epoch_idx]
    
    # Create the plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle(f'Training Loss Curves Analysis ({mode_label})', 
                 fontsize=16, fontweight='bold')
    
    # Shared legend label
    legend_label = f'Saved checkpoint (Epoch {saved_epoch})'
    
    # Total Loss
    axes[0, 0].plot(epochs, loss, 'b-', linewidth=2)
    axes[0, 0].scatter(saved_epoch, saved_loss, color='red', s=40, 
                       zorder=5, label=legend_label, marker='o', edgecolors='darkred', linewidths=2)
    axes[0, 0].axhline(y=saved_loss, color='red', linestyle='--', alpha=0.3, linewidth=1.5)
    axes[0, 0].set_title('Cross + Profile + BCR Loss', fontsize=12, fontweight='bold')
    axes[0, 0].set_xlabel('Epoch')
    axes[0, 0].set_ylabel('Loss')
    axes[0, 0].grid(True, alpha=0.3)
    
    # Cross Loss (CMC Loss)
    axes[0, 1].plot(epochs, loss_cmc, 'g-', linewidth=2)
    axes[0, 1].scatter(saved_epoch, saved_cmc, color='red', s=40, 
                       zorder=5, marker='o', edgecolors='darkred', linewidths=2)
    axes[0, 1].axhline(y=saved_cmc, color='red', linestyle='--', alpha=0.3, linewidth=1.5)
    axes[0, 1].set_title('Cross Loss', fontsize=12, fontweight='bold')
    axes[0, 1].set_xlabel('Epoch')
    axes[0, 1].set_ylabel('Loss')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Profile Loss (P2P Loss)
    axes[1, 0].plot(epochs, loss_p2p, 'orange', linewidth=2)
    axes[1, 0].scatter(saved_epoch, saved_p2p, color='red', s=40, 
                       zorder=5, marker='o', edgecolors='darkred', linewidths=2)
    axes[1, 0].axhline(y=saved_p2p, color='red', linestyle='--', alpha=0.3, linewidth=1.5)
    axes[1, 0].set_title('Profile Loss', fontsize=12, fontweight='bold')
    axes[1, 0].set_xlabel('Epoch')
    axes[1, 0].set_ylabel('Loss')
    axes[1, 0].grid(True, alpha=0.3)
    
    # BCR Loss (B2B Loss)
    axes[1, 1].plot(epochs, loss_b2b, 'purple', linewidth=2)
    axes[1, 1].scatter(saved_epoch, saved_b2b, color='red', s=40, 
                       zorder=5, marker='o', edgecolors='darkred', linewidths=2)
    axes[1, 1].axhline(y=saved_b2b, color='red', linestyle='--', alpha=0.3, linewidth=1.5)
    axes[1, 1].set_title('BCR Loss', fontsize=12, fontweight='bold')
    axes[1, 1].set_xlabel('Epoch')
    axes[1, 1].set_ylabel('Loss')
    axes[1, 1].grid(True, alpha=0.3)
    
    # Add shared legend at the top
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.96), 
               ncol=1, fontsize=11, frameon=True, shadow=True)
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Save figure if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"\n✅ Figure saved to: {save_path}")
    
    return fig

