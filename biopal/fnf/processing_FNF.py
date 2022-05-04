# SPDX-FileCopyrightText: BioPAL <biopal@esa.int>
# SPDX-License-Identifier: MIT
"""
Processing (core) functions for the BIOMASS Forest / Non-Forest map generation

@author: Alberto Alonso-Gonzalez
"""

import numpy as np



def KPSVD_singular_values(Tm, npol=3):
    """
    Performs the Kronecker Product Singular Value Decomposition (KPSVD) of the
    provided MPMB covariance matrix.
    The full set of singular values (Si) is returned.
    
    Note that no whitening is applied here as it performs better over surfaces
    when the scattering level at some polarizations is close to the NESZ.
    
    Parameters
    ----------
    Tm : ndarray, shape (Naz, Nrg, npol*N, npol*N)
        Multi-polarization multi-baseline covariance matrix. The matrix for
        each pixel is located on the last 2 dimensions (npol*N, npol*N) where
        npol is the number of polarizations and N the number of acquisitions
    npol : int, optional, default=3
        Number of polarizations available
    
    Returns
    -------
    Si : ndarray, shape (Naz, Nrg, Ns)
        Vector of singular values for each pixels, sorted in desending order
    """
    Np = Tm.shape[-1]
    Nbaselines = Np // npol
    Ns = np.minimum(npol*npol, Nbaselines*Nbaselines)
    
    Si = np.empty(Tm.shape[:2] + (Ns, ), dtype=Tm.real.dtype)
    # Process line by line instead of full image at once to avoid
    # too much memory consumption on very large images
    # --> Is this really required?
    for i, Tline in enumerate(Tm):
        # Using reshape + transpose + reshape trick to get tilde matrices fast
        Tmt = np.transpose(Tline.reshape(-1, Nbaselines, npol, Nbaselines, npol),
                           (0,1,3,2,4)).reshape(-1, Nbaselines*Nbaselines, npol*npol)
        # Do not compute singular vectors to improve speed
        S = np.linalg.svd(Tmt, compute_uv=False)
        Si[i,...] = S
    return Si

def structure_matrix_algebraic_measures(Rm):
    """
    Get Coherence standard deviation and Probability of the first eigenvalue
    of the structure matrix measures from the given multibaseline covariance
    matrix Rm.

    Parameters
    ----------
    Rm : ndarray, shape (Naz, Nrg, N, N)
        Multi-baseline covariance matrix. The matrix for each pixel is located
        on the last 2 dimensions (N, N) where N is the number of acquisitions.

    Returns
    -------
    coh_astd : ndarray, shape (Naz, Nrg)
        Coherence standard deviation measure for each pixel.
    pi : ndarray, shape (Naz, Nrg)
        Probability of the first eigenvalue of the structure matrix measure
        for each pixel.

    """
    coh_astd = np.zeros(Rm.shape[:-2])
    pi = np.zeros(Rm.shape[:-2])
    # Mask to get all the upper diagonal coherences
    msku = np.triu(np.ones(np.asarray(Rm.shape[-1]), dtype=bool), 1)
    # Process line by line instead of full image at once to avoid
    # too much memory consumption on very large images
    # --> Is this really required?
    for i, Rline in enumerate(Rm):
        # Get the diagonal ^ -0.5
        di = 1 / np.sqrt(np.einsum('...ii->...i', Rline.real))
        # Generate whitening matrix
        Ti = np.einsum("...i,...j->...ij", di, di)
        # Build the matrix of coherences - Structure matrix
        coh = Ti * Rline
        # Set invalid values to 0 to avoid linalg problems (required?)
        coh[np.isnan(coh)] = 0
        Wl = np.linalg.eigvalsh(coh)
        # Coherence standard deviation
        coh_astd[i,...] = np.std(np.abs(coh[...,msku]), axis=-1)
        # Probability of the first eigenvalue of the structure matrix
        pi[i,...] = Wl[...,-1] / np.sum(Wl,axis=-1)
    return coh_astd, pi

def get_centered_logreg_measures(Tm):
    """
    Get all the centered features employed by logistic regression from the
    given Multi-polarization multi-baseline covariance matrix.

    Parameters
    ----------
    Tm : ndarray, shape (Naz, Nrg, npol*N, npol*N)
        Multi-polarization multi-baseline covariance matrix. The matrix for
        each pixel is located on the last 2 dimensions (npol*N, npol*N) where
        npol is the number of polarizations and N the number of acquisitions.

    Returns
    -------
    Msv1 : ndarray, shape (Naz, Nrg)
        Centered Relative power of the first KPSVD component.
    Mcohstd : ndarray, shape (Naz, Nrg)
        Centered Coherence standard deviation.
    Mp1 : ndarray, shape (Naz, Nrg)
        Centered Probability of the first eigenvalue of the structure matrix.
    Mp1p2coh : ndarray, shape (Naz, Nrg)
        Centered Pauli1 - Pauli2 coherence.
    Mhhvvcoh : ndarray, shape (Naz, Nrg)
        Centered HH-VV coherence.
    Mspan : ndarray, shape (Naz, Nrg)
        Centered Trace (span) mean.

    """
    # Matrix to pass from Pauli to Lexicographic (HH,HV,VV) basis
    N = np.asarray([[ 1,  0,  1],
                    [ 1,  0, -1],
                    [ 0,  2,  0]]) / np.sqrt(2)
    Si = KPSVD_singular_values(Tm)
    # Relative power of the first KPSVD component
    sv1 = Si[...,0] / np.sum(Si, axis=-1)
    # Get measures from the Pauli1 structure matrix
    coh_astd, pi = structure_matrix_algebraic_measures(Tm[..., ::3, ::3])
    # Get mean coherency matrix between all acquisitions
    mcov = np.mean([Tm[..., i:i+3, i:i+3] for i in range(0, Tm.shape[-1], 3)], axis=0)
    # Pauli1 - Pauli2 coherence
    p12coh = np.abs(mcov[...,0,1])/np.sqrt(mcov[...,0,0].real*mcov[...,1,1].real)
    # Get mean covariance matrix between all acquisitions
    mcov = np.mean([N.T @ Tm[..., i:i+3, i:i+3] @ N for i in range(0, Tm.shape[-1], 3)], axis=0)
    # HH-VV coherence
    hhvvcoh = np.abs(mcov[...,0,2])/np.sqrt(mcov[...,0,0].real*mcov[...,2,2].real)
    # Trace mean
    trimg = 10*np.log10(np.einsum('...ii->...', Tm.real) / (Tm.shape[-1]//3))
    # Return centered measures
    return (2*sv1 - 1,       4*coh_astd - 1,   2*pi - 1,
            2*p12coh - 1,    2*hhvvcoh - 1,    trimg/20 + 0.5)

def logistic_regression(x, intercept, weights):
    """
    Perform logistic regression to estimate probability.
    This corresponds to the logistic function applied to the linear
    combination of the given features multiplied by weights plus intercept.

    Parameters
    ----------
    x : ndarray, shape (..., nfeatures)
        Input data with nfeatures for each pixel.
    intercept : ndarray, shape (...)
        Intercept value for logistic regression.
    weights : ndarray, shape (..., nfeatures)
        Logistic regression coefficients.

    Returns
    -------
    p : ndarray, shape (...)
        Estimated probability for each pixel.

    """
    # This exponential is equivalent to this one:
    # slm = np.exp(intercept + np.sum(x*weights, axis=-1))
    # return (slm/(1+slm))
    return 1 / (1 + np.exp(-intercept - np.sum(x*weights, axis=-1)))

def estimate_pixel_coefficients(kzs, log_reg_coeffs, log_reg_maxkz):
    """
    Estimate the logistic regression coefficients for each pixel based on its
    vertical wavenumber according to the coefficients table.

    Parameters
    ----------
    kzs : ndarray, shape (Naz, Nrg, N)
        Vertical wavenumber for each pixel and acquisition.
    log_reg_coeffs : ndarray, shape (K, 7)
        Logistic regression coefficients table.
        It provides the optimum coefficients at K different vertical
        wavenumbers.
    log_reg_maxkz : ndarray, shape (K)
        Maximum vertical wavenumber corresponding to each ot the K lines of
        the logistic regression table.

    Returns
    -------
    pixel_coefs : ndarray, shape (Naz, Nrg, 7)
        The logistic regression weights for each pixel.

    """
    # Prepare input data for np.interp:
    # it requires input x coordinates to be in increasing order
    idx = np.argsort(log_reg_maxkz)
    log_reg_maxkz = log_reg_maxkz[idx]
    log_reg_coeffs = log_reg_coeffs[idx,:]
    # List to contain an image for each of the coefficients
    pixel_coefs = []
    max_kz = kzs.max(axis=-1) - kzs.min(axis=-1)
    # Linear interpolate for each of the coefficients
    for i in range(log_reg_coeffs.shape[-1]):
        pixel_coefs.append(np.interp(max_kz, log_reg_maxkz, log_reg_coeffs[:,i]))
    # Generate a numpy array with all the coefficients
    pixel_coefs = np.stack(pixel_coefs, axis=-1)
    return pixel_coefs

