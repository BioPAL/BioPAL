# SPDX-FileCopyrightText: BioPAL <biopal@esa.int>
# SPDX-License-Identifier: MIT

import logging
import numpy as np
from scipy.stats import chi2


def get_rho_omega(j, p, n):
    # j is the number of covariance matrices being considered
    # p is the dimension of the covariance matrix (.ie 2 for dual pol, 3 for full pol)
    # n is the number of looks
    rho = 1 - (2 * p ** 2 - 1) * (1 + 1 / (j * (j - 1))) / (6 * p * n)
    omega2 = -(p ** 2) * (1 - 1 / rho) ** 2 / 4 + p ** 2 * (p ** 2 - 1) * (
        1 + (2 * j - 1) / (j ** 2 * (j - 1) ** 2)
    ) / (24 * n ** 2 * rho ** 2)
    return rho, omega2


def funR_wishart_SU(Y, X, j, n, p):

    # Here X is the next time point and Y is the sum of the previous X's after
    # the last change was detected. j is the total number of data being
    # considered ie. those which have contributed to Y plus one.

    lnR = n * (
        p * (j * np.log(j) - (j - 1) * np.log(j - 1))
        + (j - 1) * np.log(np.linalg.det(Y))
        + np.log(np.linalg.det(X))
        - j * np.log(np.linalg.det(Y + X))
    )
    lnR = np.real(lnR)
    return lnR


def alg_wishart_SU(X, Y, j, p, n, alpha, verbose):

    change = 0
    rho, omega2 = get_rho_omega(j + 1, p, n)
    testStat = -2 * rho * funR_wishart_SU(Y, X, j + 1, n, p)

    if testStat > 0:
        probR = 1 - (
            chi2.cdf(testStat, p ** 2) + omega2 * (chi2.cdf(testStat, p ** 2 + 4) - chi2.cdf(testStat, p ** 2))
        )
        if verbose:
            logging.info("Marginal Hypothesis %0.4f {} \n".format(probR))

        prob = 1 - probR
        if probR < alpha:
            # Change detected
            Y = X
            j = 1
            change = 1
        else:
            Y = Y + X
            j = j + 1
    else:
        prob = 0

    return Y, j, change, prob, testStat


def _apply_index_(a, idx, axis=-1):
    '''
    Funtion to apply argsort indexes to an array over an axis.
    Modified from:
        https://stackoverflow.com/questions/11253495/numpy-applying-argsort-to-an-array
    '''
    i = list(np.ogrid[[slice(x) for x in a.shape]])
    i[axis] = idx
    return a[tuple(i)]
    

def generalized_eigendecomp(A, B, eps=1e-6):
    """
    Computes the generalized eigendecomposition of the two covariance
    matrices A and B.
    
    :math:`\mathbf{A} \mathbf{v} = \lambda \mathbf{B} \mathbf{v}`
    
    This functions accepts numpy-like broadcasting.

    Parameters
    ----------
    A : ndarray, shape (..., N, N)
        Covariance matrix or matrices A to perform the generalized
        eigendecomposition
    B : ndarray, shape (..., N, N)
        Covariance matrix or matrices B to perform the generalized
        eigendecomposition
    eps : float, optional
        Preloading to apply to the main diagonal to avoid numerical errors
        when covariance matrices are semi-definite or zeros.
        The default value is 1e-6.
        It can be set to 0 to disable preloading.

    Returns
    -------
    W : ndarray, shape (..., N)
        Generalized eigenvalues (real) between A and B matrices
    V : ndarray, shape (..., N, N)
        Generalized eigenvectors between A and B matrices
    
    References
    ----------
    Alonso González, A., López Martínez, C., Papathanassiou, K., & Hajnsek, I.
    (2020). Polarimetric SAR time series change analysis over agricultural
    areas. IEEE transactions on geoscience and remote sensing, 58(10),
    7317-7330.

    """
    # Apply main diagonal preloading
    if eps > 0:
         A = A + eps * np.eye(A.shape[-1])
         B = B + eps * np.eye(A.shape[-1])
    # Perform generalized eigendecompsition
    W,V = np.linalg.eig(np.linalg.solve(A, B))
    # Note: eigenvalues/vectors are not sorted, we need to sort them
    ind = np.argsort(W.real, axis=-1)
    # Apply argsort to eigenvalues & eigenvectors
    W = _apply_index_(W.real, ind, axis=-1)
    V = _apply_index_(V, ind[...,np.newaxis,:], axis=-1)
    return W, V
    

def generate_pinc_pdec_RGBimages(W, V, p_min=1, p_max=10):
    """
    Generate polarimetric change increase and decrease Pauli RGB images.

    Parameters
    ----------
    W : ndarray, shape (..., N)
        Generalized eigenvalues (should be real).
    V : ndarray, shape (..., N, N)
        Generalized eigenvectors.
    p_min : float, optional
        Minimum value (in dB) for the scaling of the increasing/decreasing 
        Pauli RGB images. The default is 1dB.
    p_max : float, optional
        Maximum value (in dB) for the scaling of the increasing/decreasing 
        Pauli RGB images. The default is 10dB.

    Returns
    -------
    pinc : ndarray, shape (..., 3)
        Change representation corresponding to the Pauli RGB of the increasing
        polarization states within the range (p_min, p_max) in dB.
    pdec : ndarray, shape (..., 3)
        Change representation corresponding to the Pauli RGB of the decreasing
        polarization states within the range (p_min, p_max) in dB.
    
    References
    ----------
    Alonso González, A., López Martínez, C., Papathanassiou, K., & Hajnsek, I.
    (2020). Polarimetric SAR time series change analysis over agricultural
    areas. IEEE transactions on geoscience and remote sensing, 58(10),
    7317-7330.

    """
    # Get generalized eigenvalues in dB
    Wl = 10*np.log10(W)
    # Compute pinc and pdec vectors
    pinc = np.linalg.norm(np.abs(V) * ((Wl * (Wl > 0))[..., np.newaxis,: ]),
                          axis=-1)
    pdec = np.linalg.norm(np.abs(V) * ((-Wl * (Wl < 0))[..., np.newaxis, :]),
                          axis=-1)
    # NOTE: The channel ordering has to be changed to match the classical
    # Pauli RGB representation
    pinc = pinc[...,(1,2,0)]
    pdec = pdec[...,(1,2,0)]
    # Clip images to the given change range
    pinc = (np.clip(pinc, p_min,p_max)-p_min)/(p_max-p_min)
    pdec = (np.clip(pdec, p_min,p_max)-p_min)/(p_max-p_min)
    return pinc, pdec
    