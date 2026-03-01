# -*- coding: utf-8 -*-
"""
Functionality for loading single-cell RNA sequencing data, surface-
based spatial brain phenotypes, parcellation and data frame.
"""

import scanpy as sc
import anndata as ad
import scipy.sparse as sp
import nibabel as nib
import numpy as np
import pandas as pd
import os
from typing import Union

def load_sc(
    adata: Union[os.PathLike,ad.AnnData],
    flag_sparse: bool
) -> ad.AnnData:

    """
    Load and check scRNA-seq file (.h5ad).

    Parameters
    ----------
    adata : str or AnnData
        Path to .h5ad file or AnnData object     
    flag_sparse : bool
        If True, enforce adata.X to be a sparse matrix

    Returns
    -------
    adata : ad.AnnData
        AnnData object containing scRNA-seq data.
    """
    if isinstance(adata, str):
        adata = sc.read_h5ad(adata)
    elif isinstance(adata, ad.AnnData):
        pass
    else:
        raise TypeError(
                f"AnnData must be a filepath (str) or ad.AnnData. Got {type(df)}."
            )
    
    # Check input:
    # 1. If there is any NaN value;
    # 2. If there is any negative value;
    if np.isnan(adata.X.sum()):
        raise ValueError(
            "Single-cell gene expression matrix should not contain any NaN value. Please check your data."
        )    
    if (adata.X < 0).sum() > 0:
        raise ValueError(
            "Single-cell gene expression matrix should not contain any negative value. Please check your data."
        )      
    
    # Convert gene expression matrix into sparse matrix
    if flag_sparse:
        if not sp.issparse(adata.X):
            adata.X = sp.csr_matrix(adata.X)    
    return adata

# A quick function to grep vertex number for each atlas
def get_vertex_num(atlas, density, hemi):
    atlas_vertex_dict = {'fsLR_32k': 32492, 
                         'fsLR_164k': 163842,
                         'fsaverage_10k': 10242, 
                         'fsaverage_41k': 40962, 
                         'fsaverage_164k': 163842,
                         'civet_41k': 40962}
    
    if hemi in ['L','lh','left']:
        return atlas_vertex_dict[f'{atlas}_{density}']
    elif hemi == 'both':
        return atlas_vertex_dict[f'{atlas}_{density}']*2

def load_gii(
    gii: Union[os.PathLike, tuple, nib.gifti.gifti.GiftiImage],
    atlas: str, density: str, hemi: str
) -> np.ndarray:
    
    """
    Load GIFTI image.

    Parameters
    ----------
    gii : path_like str, tuple or nib.gifti.gifti.GiftiImage
        Recognized GIFTI image
        - If hemi in ['L','lh','left']:
            - path_like str: Filepaths to .gii
            - nib.gifti.gifti.GiftiImage
        - Otherwise, use both hemisphere:
            - tuple: must be organized as L/R pair
                - (str, str): Filepaths to .gii
                - (GiftiImage, GiftiImage)
    atlas : {'fsaverage', 'fSLR'} optional,
        Name of surface atlas on which `gii` is defined.
    density : str, optional
        Density of surface mesh on which `gii` is defined. Must becompatible 
        with specified `atlas`.
    hemi : {'L','lh','left','both'} optional,
        Hemisphere used to perform downstream analyses.

    Returns
    -------
    gii_data : np.ndarray
        Data extracted from given GIFTI image.
    """
    
    # Fetch expected vertex number for given space
    vertex_num = get_vertex_num(atlas, density, hemi)
    
    # Load GIFTI image for left hemisphere analysis
    if hemi in ['L','lh','left']:
        if isinstance(gii, str):
            gii_img = nib.load(gii)
        elif isinstance(gii, nib.gifti.gifti.GiftiImage):
            gii_img = gii
        else:
            raise TypeError(
                f"For left hemisphere data, GIFTI image must be a filepath (str) or GiftiImage. Got {type(gii)}."
            )
        gii_data = gii_img.darrays[0].data
    
    # Load GIFTI image for both hemisphere analysis
    elif hemi == 'both':
        gii_data_lr = []
        if isinstance(gii, tuple):
            if not gii or len(gii) != 2: # Assume tuples are for L/R pair
                raise ValueError("Tuple must contain exactly two elements (for L and R).")            
            # Determine type of elements in tuple
            str_tuple = all(isinstance(p, str) and p.endswith('.gii') for p in gii)
            gii_tuple = all(isinstance(p, nib.gifti.gifti.GiftiImage) for p in gii)
            if str_tuple:
                lh_gii_img = nib.load(gii[0])
                rh_gii_img = nib.load(gii[1])
                gii_data_lr.append(lh_gii_img.darrays[0].data)
                gii_data_lr.append(rh_gii_img.darrays[1].data)
            elif gii_tuple:
                gii_data_lr.append(gii[0].darrays[0].data)
                gii_data_lr.append(gii[1].darrays[1].data)
            else:
                raise TypeError(
                    f"For both hemisphere data, GIFTI image must be a tuple containg filepaths or GIFTI images. Got {type(gii[0])} and {type(gii[1])}."
                )
        else:
            raise TypeError(
                f"For both hemisphere data, GIFTI image must be a tuple containg filepaths or GIFTI images. Got {type(gii)}."
            )
        gii_data = np.concatenate(gii_data_lr)
        
    if gii_data.shape[0] != vertex_num:
        raise ValueError(f"For the {atlas} {density} ({hemi} hemisphere) standard space, the expected number of vertices is {vertex_num}, whereas the actual number obtained is {gii_data.shape[0]}.")        
    return gii_data
    
def load_parc(
    parcellation: Union[os.PathLike, tuple, nib.gifti.gifti.GiftiImage],
    atlas: str, density: str, hemi: str
) -> np.ndarray:
    
    """
    Load parcellation image.

    Parameters
    ----------
    parcellation : path_like str, tuple, nib.gifti.gifti.GiftiImage
        Recognized parcellation image
        - If hemi in ['L','lh','left']:
            - path_like str: Filepaths to .gii or .annot
            - nib.gifti.gifti.GiftiImage
        - Otherwise, use both hemisphere:
            - tuple: must be organized as L/R pair
                - (str, str): Filepaths to .gii or .annot
                - (nib.gifti.gifti.GiftiImage, nib.gifti.gifti.GiftiImage)
    atlas: {'fsaverage', 'fSLR', 'civet'} optional,
        Name of surface atlas on which `parcellation` is defined.
    density : str, optional
        Density of surface mesh on which `parcellation` is defined. Must becompatible 
        with specified `atlas`.
    hemi: {'L','lh','left','both'} optional,
        Hemisphere used to perform downstream analyses.

    Returns
    -------
    gii_data : np.ndarray
        Data extracted from given GIFTI image.
    """
    
    # Fetch expected vertex number for given space
    vertex_num = get_vertex_num(atlas, density, hemi)
    
    # Load parcellation image for left hemisphere analysis
    if hemi in ['L','lh','left']:
        if isinstance(parcellation, str):
            if parcellation.endswith('.gii'):
                parc_data = load_gii(parcellation,atlas,density,hemi)
            elif parcellation.endswith('.annot'):
                parc_data, _, _ = nib.freesurfer.io.read_annot(parcellation)
            else:
                raise ValueError(
                    f"Filepath must end with .gii or .annot. Got: {parcellation}"
                )
        elif isinstance(parcellation, nib.gifti.gifti.GiftiImage):
            parc_data = load_gii(parcellation,atlas,density,hemi)
        else:
            raise TypeError(
                f"For left hemisphere, parcellation image must be a filepath (str) or GiftiImage. Got {type(parcellation)}."
            )
            
    # Load parcellation image for both hemisphere analysis
    elif hemi == 'both':
        if isinstance(parcellation, tuple):
            if not parcellation or len(parcellation) != 2: # Assume tuples are for L/R pair
                raise ValueError("Tuple must contain exactly two elements (for L and R).")            
            # Determine type of elements in tuple
            str_gii_tuple = all(isinstance(p, str) and p.endswith('.gii') for p in parcellation)
            str_annot_tuple = all(isinstance(p, str) and p.endswith('.annot') for p in parcellation)
            gii_tuple = all(isinstance(p, nib.gifti.gifti.GiftiImage) for p in parcellation)
            
            if str_gii_tuple or gii_tuple:
                parc_data = load_gii(parcellation,atlas,density,hemi)
            elif str_annot_tuple:
                lh_parc_data, _, _ = nib.freesurfer.io.read_annot(parcellation[0])
                rh_parc_data, _, _ = nib.freesurfer.io.read_annot(parcellation[1])
                parc_data = np.concatenate([lh_parc_data,rh_parc_data])
            else:
                raise TypeError(
                    f"For both hemisphere data, parcellation image must be a tuple containg filepaths ending with '.gii' or '.annot', or GIFTI images. Got {type(gii[0])} and {type(gii[1])}."
                )
        else:
            raise TypeError(
                f"For both hemisphere, parcellation image must be a tuple containg filepaths or GIFTI images. Got {type(parcellation)}."
            )
        
    if parc_data.shape[0] != vertex_num:
        raise ValueError(f"For the {atlas} {density} ({hemi} hemisphere) standard space, the expected number of vertices is {vertex_num}, whereas the actual number obtained is {parc_data.shape[0]}.")        
    return parc_data

def load_df(df: Union[os.PathLike, pd.DataFrame]) -> pd.DataFrame:
    
    """
    Load parcellation image.

    Parameters
    ----------
    file : path_like str, pd.DataFrame
        Recognized data frame file
        - path_like str: Filepath to .csv, .tsv
        - pd.DataFrame

    Returns
    -------
    df : pd.DataFrame
        Data frame for analysis.
    """
    
    if isinstance(df, str):
        if df.endswith('.csv'):
            df = pd.read_csv(df, index_col=0)
        elif df.endswith('.tsv'):    
            df = pd.read_csv(df, sep='\t', index_col=0)
        else:
            raise ValueError(
                    f"Filepath must end with .csv or .tsv. Got: {df}"
                )
    elif isinstance(df, pd.DataFrame):
        pass
    else:
        raise TypeError(
                f"DataFrame must be a filepath (str) or pd.DataFrame. Got {type(df)}."
            )    
    return df