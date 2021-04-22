import numpy as np
from biopal.fd.processing_FD import generalized_eigendecomp
from biopal.fd.processing_FD import _apply_index_
from biopal.fd.processing_FD import generate_pinc_pdec_RGBimages
from scipy.linalg import eig as sp_eig

def generate_cov_matrices_change(contrast_vector):
    '''
    Generates a simple change by scaling the eigenvalues of a randomly
    generated covariance matrix.
    This means that both generated A snd B matrices will share eigenvectors,
    i.e. they will commute.

    Parameters
    ----------
    contrast_vector : sequence of length p
        The scaling (factor) to the individual eigenvalues
    
    Returns
    -------
    A : ndarray, shape (p, p)
        Random covariance matrix
    B : ndarray, shape (p, p)
        A matrix with eigenvalues scaled according to contrast_vector

    '''
    p = len(contrast_vector)
    # Generate random A matrix
    A = np.random.randn(p, p) + 1j * np.random.randn(p, p)
    # Make it Hermitian & positive definite
    A = A.T.conj().dot(A)
    
    # Generate a change on matrix A by scaling its eigenvalues
    W, V = np.linalg.eigh(A)
    # Scale eigenvalues by given factor
    W *= contrast_vector
    # Construct matrix B with given changes
    B = V.dot(np.diag(W)).dot(V.T.conj())
    
    return A, B

def generate_random_cov_matrices(p=3):
    '''
    Generates random covariance matrices with the given size.

    Parameters
    ----------
    p : integer, optional
        The size of the covariance matrix, the default is 3.

    Returns
    -------
    A : ndarray, shape (p, p)
        Random covariance matrix
    B : ndarray, shape (p, p)
        Random covariance matrix

    '''
    # Generate random A matrix
    A = np.random.randn(p, p) + 1j * np.random.randn(p, p)
    # Make it Hermitian & positive definite
    A = A.T.conj().dot(A)
    # Generate random B matrix
    B = np.random.randn(p, p) + 1j * np.random.randn(p, p)
    # Make it Hermitian & positive definite
    B = B.T.conj().dot(B)
    return A, B

def test_generalized_eigenvalue_commutingAB():
    '''
    Test the generalized eigendecomposition against a generated simple
    change on dummy covariance matrices.
    The increase and decrease images p_inc and p_dec are also tested for the
    generated change.

    Returns
    -------
    None.

    '''
    # Polarimetric contrasts to generate
    contrast_vector = [1, 2, 3]
    # Generate random covariance matrices with the given change
    A, B = generate_cov_matrices_change(contrast_vector)
    # Call generalized eigendecomp
    W, V = generalized_eigendecomp(A, B, eps=0)
    # Check if the ontained contrast matches the generated one
    assert np.allclose(np.sort(W), np.sort(contrast_vector)), "The obtained generalized eigenvalues do not match with generated ones"
    
    
    # Now check the generation of pinc and pdec images
    pmin = 0
    pmax = np.max(10*np.log10(contrast_vector))
    # Generate pinc and pdec using BioPAL function
    pinc, pdec = generate_pinc_pdec_RGBimages(W, V, p_min=pmin, p_max=pmax)
    
    # Generate reference pinc and pdec RGB images manually
    pinc_ref = np.zeros_like(pinc)
    pdec_ref = np.zeros_like(pdec)
    # Since both matrices share eigenvectors the generalized eigenvectors
    # are the same than each covariance matrix eigenvectors in this case
    Wa, Va = np.linalg.eigh(A)
    for i in range(3):
        lw = 10*np.log10(contrast_vector[i])
        if lw > 0:
            # increase
            lw = np.clip((lw - pmin) / (pmax - pmin), 0, None)
            pinc_ref += (lw*np.abs(Va[(1, 2, 0), i]))**2
        else:
            # decrease
            lw = np.clip((-lw - pmin) / (pmax - pmin), 0, None)
            pdec_ref += (lw*np.abs(Va[(1, 2, 0), i]))**2
    pinc_ref = np.clip(np.sqrt(pinc_ref), 0, 1)
    pdec_ref = np.clip(np.sqrt(pdec_ref), 0, 1)
    
    assert np.allclose(pinc, pinc_ref), "The obtained increase RGB image p_inc does not match with the reference one"
    assert np.allclose(pdec, pdec_ref), "The obtained decrease RGB image p_dec does not match with the reference one"
    
def test_generalized_eigenvalue_scipy():
    '''
    Test the generalized eigendecomposition with respect to the scipy
    implementation with randomly generated matrices.
    This function also test if broadcasting works properly to compute 2D image
    blocks at once.

    Returns
    -------
    None.

    '''
    # Size of covariance matrices
    p = 3
    # Generate random matrices
    A, B = generate_random_cov_matrices(p)
    # Compute generalized eigenvalues
    W, V = generalized_eigendecomp(A, B, eps=0)
    # Compute generalized eigenvalues with scipy (NOTE the different order of 
    # the parameters)
    Ws, Vs = sp_eig(B, A)
    # Note: eigenvalues/vectors are not sorted, we need to sort them
    ind = np.argsort(Ws.real, axis=-1)
    # Apply argsort to eigenvalues & eigenvectors
    Ws = _apply_index_(Ws.real, ind, axis=-1)
    Vs = _apply_index_(Vs, ind[...,np.newaxis,:], axis=-1)
    
    # Check if generalized eigenvalues match
    assert np.allclose(W, Ws), "The generalized eigenvalues obtained do not match with scipy implementation"
    # Check if generalized eigenvectors match, up to a phase change
    assert np.allclose(np.abs(V), np.abs(Vs)), "The generalized eigenvectors (in absolute value) obtained do not match with scipy implementation"
    
    # Now make the test with a 2D image, in order to test if the broadcasting
    # rules work properly to compute image blocks at once
    img_shape = (20, 10)
    A = np.zeros(img_shape + (p,p), dtype=np.complex)
    B = np.zeros_like(A)
    Ws = np.zeros(img_shape + (p,), dtype=np.complex)
    Vs = np.zeros_like(A)
    
    # Generate a 2D image of random covariance matrices
    for i in range(img_shape[0]):
        for j in range(img_shape[1]):
            A[i, j], B[i, j] = generate_random_cov_matrices(p)
            Wsi, Vsi = sp_eig(B[i, j], A[i, j])
            # Note: eigenvalues/vectors are not sorted, we need to sort them
            ind = np.argsort(Wsi.real, axis=-1)
            # Apply argsort to eigenvalues & eigenvectors
            Wsi = _apply_index_(Wsi.real, ind, axis=-1)
            Vsi = _apply_index_(Vsi, ind[...,np.newaxis,:], axis=-1)
            Ws[i, j], Vs[i, j] = Wsi, Vsi
    
    # Compute generalized eigendecomp for the whole image at once
    W, V = generalized_eigendecomp(A, B, eps=0)
    # Check if generalized eigenvalues match
    assert np.allclose(W, Ws), "The generalized eigenvalues obtained for a 2D image block do not match with scipy implementation"
    # Check if generalized eigenvectors match, up to a phase change
    assert np.allclose(np.abs(V), np.abs(Vs)), "The generalized eigenvectors (in absolute value) obtained for a 2D image block do not match with scipy implementation"
    
    