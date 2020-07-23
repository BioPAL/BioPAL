import numpy as np
from   scipy.sparse        import csr_matrix
from   scipy.signal        import firwin, oaconvolve
from   scipy.constants     import c as LIGHTSPEED
from   biomassL2.constants import RANGE_BANDWIDTH_HZ, CARRIER_FREQUENCY_HZ
import logging


def main_covariance_estimation_SR(data_stack, cov_est_window_size, pixel_spacing_slant_rg, pixel_spacing_az, incidence_angle_rad):
    """ Covariance estimation:  inputs are in Slant Range radar coordinates 
        see also main_covariance_estimation_GR """
    
    # compute the pixels spacing in ground range, from slant range:
    pixel_spacing_grd_x = pixel_spacing_slant_rg/np.sin(incidence_angle_rad)
    pixel_spacing_grd_y = pixel_spacing_az
    
    slant_range_spacing_s = pixel_spacing_slant_rg/LIGHTSPEED
    
    MPMB_covariance, rg_vec_subs, az_vec_subs, subs_F_r, subs_F_a = main_covariance_estimation_GR(data_stack, cov_est_window_size, pixel_spacing_grd_x, pixel_spacing_grd_y, slant_range_spacing_s)

    return MPMB_covariance,  rg_vec_subs,  az_vec_subs, subs_F_r, subs_F_a


def main_covariance_estimation_SSF_SR(data_stack, cov_est_window_size, pixel_spacing_slant_rg, pixel_spacing_az, incidence_angle_rad, R, look_angles, ground_slope):
    """ Covariance estimation:  inputs are in Slant Range radar coordinates 
        see also main_covariance_estimation_GR """
    
    # compute the pixels spacing in ground range, from slant range:
    pixel_spacing_grd_x = pixel_spacing_slant_rg/np.sin(incidence_angle_rad)
    pixel_spacing_grd_y = pixel_spacing_az
    
    slant_range_spacing_s = pixel_spacing_slant_rg/LIGHTSPEED
    
    MPMB_covariance,  rg_vec_subs,  az_vec_subs, subs_F_r, subs_F_a = main_covariance_estimation_SSF_GR(data_stack, cov_est_window_size, pixel_spacing_grd_x, pixel_spacing_grd_y, R, look_angles, ground_slope, slant_range_spacing_s)

    return MPMB_covariance, rg_vec_subs, az_vec_subs, subs_F_r, subs_F_a


def main_covariance_estimation_GR(data_stack, cov_est_window_size, pixel_spacing_grd_x, pixel_spacing_grd_y, slant_range_spacing_s=None):
    """ Covariance estimation:  inputs are in Ground Range coordinates 
        see also main_covariance_estimation_SR """
    
    cov_est_opt_str = covariance_options_struct_filler(data_stack, cov_est_window_size, pixel_spacing_grd_x, pixel_spacing_grd_y, slant_range_spacing_s)
    
    MPMB_covariance, rg_vec_subs, az_vec_subs, subs_F_r, subs_F_a = MPMBCovarianceEstimation(data_stack, cov_est_opt_str)
    
    return MPMB_covariance, rg_vec_subs, az_vec_subs, subs_F_r, subs_F_a 

 
def main_covariance_estimation_SSF_GR(data_stack, cov_est_window_size, pixel_spacing_grd_x, pixel_spacing_grd_y, R, look_angles, ground_slope, slant_range_spacing_s=None):
    """ Covariance estimation:  inputs are in Ground Range coordinates 
        see also main_covariance_estimation_SR """
        
    cov_est_opt_str = covariance_options_struct_filler(data_stack, cov_est_window_size, pixel_spacing_grd_x, pixel_spacing_grd_y, slant_range_spacing_s)
    
    MPMB_covariance,  rg_vec_subs,  az_vec_subs, subs_F_r, subs_F_a = MPMBCovarianceEstimationSSF(data_stack, cov_est_opt_str, R, look_angles, ground_slope )
    
    return MPMB_covariance, rg_vec_subs, az_vec_subs, subs_F_r, subs_F_a  


def covariance_options_struct_filler(data_stack, cov_est_window_size, pixel_spacing_grd_x, pixel_spacing_grd_y, slant_range_spacing_s=None):
    
    first_level_key = next(iter(data_stack.keys()))
    
    if first_level_key == 'hh' or first_level_key == 'hv' or first_level_key == 'vh' or first_level_key == 'vv':
        # notched stack, has not acquisitions dict
        pol_names      = list(  data_stack.keys() )
        Num_of_pols    = len( pol_names )
        
    else:
        acq_names      = list( data_stack.keys() )
        first_acq_dict = data_stack[ acq_names[0] ]
        pol_names      = list(  first_acq_dict.keys() )
        Num_of_pols    = len( pol_names )
    
    class cov_est_opt_str:
        pass
    
    cov_est_opt_str.NW_r              = np.int64( cov_est_window_size / pixel_spacing_grd_x / 2 ) * 2 + 1
    cov_est_opt_str.NW_a              = np.int64( cov_est_window_size / pixel_spacing_grd_y / 2 ) * 2 + 1 
    
    cov_est_opt_str.NWmono_r          = np.floor( (cov_est_opt_str.NW_r - 1)/2 ).astype('int64')
    cov_est_opt_str.NWmono_a          = np.floor( (cov_est_opt_str.NW_a - 1)/2 ).astype('int64')
    
    cov_est_opt_str.subs_F_r          = cov_est_opt_str.NWmono_r
    cov_est_opt_str.subs_F_a          = cov_est_opt_str.NWmono_a
    
    cov_est_opt_str.polarimetric_mask = np.ones((Num_of_pols, Num_of_pols))
    
    
    cov_est_opt_str.wavelenght = LIGHTSPEED/CARRIER_FREQUENCY_HZ
    cov_est_opt_str.B  = RANGE_BANDWIDTH_HZ
    cov_est_opt_str.f0 = CARRIER_FREQUENCY_HZ
    cov_est_opt_str.dt = slant_range_spacing_s
    
    return cov_est_opt_str


def main_correlation_estimation_SR(data_stack, cov_est_window_size, pixel_spacing_slant_rg, pixel_spacing_az, incidence_angle_rad):

    MPMB_covariance,  rg_vec_subs,  az_vec_subs, subs_F_r, subs_F_a  = main_covariance_estimation_SR(data_stack, cov_est_window_size, pixel_spacing_slant_rg, pixel_spacing_az, incidence_angle_rad)
    
    # Normalizing covariance matrix
    MPMB_correlation, norm_diag = Covariance4D2Correlation4D(MPMB_covariance)
    
    return MPMB_correlation,  rg_vec_subs,  az_vec_subs, subs_F_r, subs_F_a 


def main_correlation_estimation_SSF_SR(data_stack, cov_est_window_size, pixel_spacing_slant_rg, pixel_spacing_az, incidence_angle_rad, R, look_angles, ground_slope):
 
    MPMB_covariance,  rg_vec_subs,  az_vec_subs, subs_F_r, subs_F_a  = main_covariance_estimation_SSF_SR(data_stack, cov_est_window_size, pixel_spacing_slant_rg, pixel_spacing_az, incidence_angle_rad, R, look_angles, ground_slope)

    # Normalizing covariance matrix
    MPMB_correlation, norm_diag = Covariance4D2Correlation4D(MPMB_covariance)
    
    return MPMB_correlation,  rg_vec_subs,  az_vec_subs, subs_F_r, subs_F_a 


def main_correlation_estimation_GR(data_stack, cov_est_window_size, pixel_spacing_grd_x, pixel_spacing_grd_y):
    
    MPMB_covariance,  rg_vec_subs,  az_vec_subs, subs_F_r, subs_F_a  = main_covariance_estimation_GR(data_stack, cov_est_window_size, pixel_spacing_grd_x, pixel_spacing_grd_y)
        
    # Normalizing covariance matrix
    MPMB_correlation, norm_diag = Covariance4D2Correlation4D(MPMB_covariance)
    
    return MPMB_correlation,  rg_vec_subs,  az_vec_subs, subs_F_r, subs_F_a 


def main_correlation_estimation_SSF_GR(data_stack, cov_est_window_size, pixel_spacing_grd_x, pixel_spacing_grd_y, R, look_angles, ground_slope):
    
    MPMB_covariance,  rg_vec_subs,  az_vec_subs, subs_F_r, subs_F_a  = main_covariance_estimation_SSF_GR(data_stack, cov_est_window_size, pixel_spacing_grd_x, pixel_spacing_grd_y, R, look_angles, ground_slope)
        
    # Normalizing covariance matrix
    MPMB_correlation, norm_diag = Covariance4D2Correlation4D(MPMB_covariance)
    
    return MPMB_correlation,  rg_vec_subs,  az_vec_subs, subs_F_r, subs_F_a 


def MPMBCovarianceEstimation(D_in, opt_str):
    """ MPMBCovarianceEstimation
        
        INPUTS:
            D_in: can be a:
                      dictionary with acquisition keys (Num_Baselines) each containing a dictionary with hh hv vh and vv 
                      polarization keys containing arrays with dimensions [Nrg x Naz]
                or (in case of ground notched input) :
                      dictionary with hh hv vh and vv 
                      polarization keys containing arrays with dimensions [Nrg x Naz] 
            opt_str:
                NWmono_r: One-sided range average window [pixel]
                subs_F_r: Range subsampling factor
                NWmono_a: One-sided azimuth average window [pixel]
                subs_F_a: Azimuth subsampling factor
                polarimetric_mask: (optional) [Npol x Npol] logical matrix
                                              stating the polarization combination to
                                              be computed
     
        OUTPUTS:
            Cov_MPMB: Multi-polarimetric multi-baseline covariance matrix
            rg_out: range subsampled axis
            az_out: azimuth subsampled axis 
    """
        
    
    first_level_key = next(iter(D_in.keys()))
    
    if first_level_key == 'hh' or first_level_key == 'hv' or first_level_key == 'vh' or first_level_key == 'vv':
        # notched stack, has not acquisitions dict, so we need to encapsulate the input in an external dummy dictionary
        D = {'single_acq': D_in}
    else:
        D = D_in
    del D_in
        
    acq_names      = list( D.keys() )
    N_imm          = len(acq_names)
    first_acq_dict = D[ acq_names[0] ]
    pol_names      = list(  first_acq_dict.keys() )
    num_pols       = len( pol_names )
    N_rg, N_az     = first_acq_dict[ pol_names[0] ].shape
    cell_data = True
    double_data = first_acq_dict[ pol_names[0] ].dtype == 'complex128'

    try:
        Mask = opt_str.polarimetric_mask
    except AttributeError:
        Mask = np.ones((num_pols, num_pols))
    
    # Range filter matrix
    class filtering_matrix_opt_str:
        pass
    filtering_matrix_opt_str.Nin = N_rg
    filtering_matrix_opt_str.NWmono = opt_str.NWmono_r
    filtering_matrix_opt_str.subs_F = opt_str.subs_F_r
    Fr, rg_out, Rnorm = build_filtering_matrix(filtering_matrix_opt_str)
    N_rg_out = rg_out.size
    Fr_normalized = Rnorm @ Fr
    
    # Azimuth filter matrix
    filtering_matrix_opt_str.Nin = N_az
    filtering_matrix_opt_str.NWmono = opt_str.NWmono_a
    filtering_matrix_opt_str.subs_F = opt_str.subs_F_a
    Fa, az_out, Anorm = build_filtering_matrix(filtering_matrix_opt_str)
    N_az_out = az_out.size
    Fa_normalized_transposed = (Anorm @ Fa).T
    
    # Init
    Cov_MPMB = np.zeros( (num_pols*N_imm, num_pols*N_imm, N_rg_out, N_az_out), dtype=np.complex64)
    
    logging.info('    MPMB covariance estimation...')
    #nan_mask = np.zeros( D[0].shape , dtype=bool) N_rg, N_az
    nan_mask = np.zeros( (N_rg, N_az, N_imm) , dtype=bool)
    
    for pol_name_curr in pol_names:
        all_acq_data = np.zeros( (N_rg, N_az, N_imm), dtype=complex )
        for acq_idx_curr, acq_curr in enumerate(acq_names):
            all_acq_data[:,:,acq_idx_curr] = D[acq_curr][pol_name_curr]
        nan_mask = nan_mask + np.isnan( all_acq_data )
      
    # this has been moved from the innner cycle (referring to PM9)
    for acq_idx_curr, acq_curr in enumerate( acq_names):
            for pol_name_curr in pol_names:
                D[acq_curr][pol_name_curr][nan_mask[:,:,acq_idx_curr]] = 0
                
    for ch_i, pol_id_i in enumerate( pol_names ):
        
        ind_i = np.arange(N_imm) + ch_i*N_imm    
        for ch_j in np.arange(ch_i, num_pols):
            pol_id_j = pol_names[ch_j]
            if Mask[ch_i,ch_j] == 1:
                print('    .')
                ind_j = np.arange(N_imm) + ch_j*N_imm
                I = np.zeros( (N_imm, N_imm, N_rg_out, N_az_out), dtype=np.complex64)
                for n_idx, acq_id_n in enumerate( acq_names ):
                    if ch_i == ch_j:
                        m_min = n_idx#?
                    else:
                        m_min = 0
                    
                    for m_idx in np.arange(m_min, N_imm):
                        acq_id_m = acq_names[m_idx]
                        if cell_data:
                            if double_data:
                                temp = D[ acq_id_n ][ pol_id_i ]*np.conjugate(D[ acq_id_m ][ pol_id_j ])
                            else:
                                temp = D[ acq_id_n ][ pol_id_i ].astype('complex128')*np.conjugate(D[ acq_id_m ][ pol_id_j ]).astype('complex128')
                        else:
                            if double_data:
                                temp = D[ acq_id_n ][ pol_id_i ]*np.conjugate(D[ acq_id_m ][ pol_id_j ])
                            else:
                                temp = D[ acq_id_n ][ pol_id_i ].astype('complex128')*np.conjugate(D[ acq_id_m ][ pol_id_j ]).astype('complex128')
                        temp = Fr_normalized @ temp
                        temp = temp @ Fa_normalized_transposed
                        I[n_idx, m_idx, :, :] = temp.astype('complex64')
                Cov_MPMB[ind_i[:, np.newaxis], ind_j[np.newaxis, :], :, :] = I

    
    # Symmetric part generation
    diag_mask = np.tile((np.eye(num_pols*N_imm).reshape(num_pols*N_imm, num_pols*N_imm, 1, 1))>0, [1, 1, N_rg_out, N_az_out])
    Cov_MPMB = Cov_MPMB + np.conjugate(np.moveaxis(Cov_MPMB, [1, 0, 2, 3], [0, 1, 2, 3]))
    Cov_MPMB[diag_mask] = Cov_MPMB[diag_mask]/2
    
    
    logging.info('    ...done.\n')
    
    return Cov_MPMB, rg_out, az_out, opt_str.subs_F_r, opt_str.subs_F_a 


def MPMBCovarianceEstimationSSF(D, opt_str, R, look_angles, ground_slope):
    """ MPMBCovarianceEstimation
        
        INPUTS:
            D: list of lists with dimensions [Num_pol x 1][Nrg x Naz x Num_Baselines] 
            R: [Nrg x Naz x Nimm] distances (to bring to baseband)
            look_angles: [Nrg x Naz x Nimm] offnadir angles
            ground_slope: [Nrg x Naz] ground slope
            opt_str:
                NWmono_r: One-sided range average window [pixel]
                subs_F_r: Range subsampling factor
                NWmono_a: One-sided azimuth average window [pixel]
                subs_F_a: Azimuth subsampling factor
                polarimetric_mask: (optional) [Npol x Npol] logical matrix
                                              stating the polarization combination to
                                              be computed
                B: signal bandwidth
                dt: sampling interval
                f0: central frequency
                lambda: wavelength
         
        OUTPUTS:
            Cov_MPMB: Multi-polarimetric multi-baseline covariance matrix
            rg_out: range subsampled axis
            az_out: azimuth subsampled axis 
    """
        
    
    acq_names      = list( D.keys() )
    N_imm          = len(acq_names)
    first_acq_dict = D[ acq_names[0] ]
    pol_names      = list(  first_acq_dict.keys() )
    num_pols       = len( pol_names )
    N_rg, N_az     = first_acq_dict[ pol_names[0] ].shape
    cell_data = True
    double_data = first_acq_dict[ pol_names[0] ].dtype == 'complex128'

    try:
        Mask = opt_str.polarimetric_mask
    except AttributeError:
        Mask = np.ones((num_pols, num_pols))
    
    # Range filter matrix
    class filtering_matrix_opt_str:
        pass
    filtering_matrix_opt_str.Nin = N_rg
    filtering_matrix_opt_str.NWmono = opt_str.NWmono_r
    filtering_matrix_opt_str.subs_F = opt_str.subs_F_r
    Fr, rg_out, Rnorm = build_filtering_matrix(filtering_matrix_opt_str)
    N_rg_out = rg_out.size
    Fr_normalized = Rnorm @ Fr
    
    # Azimuth filter matrix
    filtering_matrix_opt_str.Nin = N_az
    filtering_matrix_opt_str.NWmono = opt_str.NWmono_a
    filtering_matrix_opt_str.subs_F = opt_str.subs_F_a
    Fa, az_out, Anorm = build_filtering_matrix(filtering_matrix_opt_str)
    N_az_out = az_out.size
    Fa_normalized_transposed = (Anorm @ Fa).T
    
    # Init
    Cov_MPMB = np.zeros( (num_pols*N_imm, num_pols*N_imm, N_rg_out, N_az_out), dtype=np.complex64)
	
    class ssf_opt_str:
        pass
    ssf_opt_str.B = opt_str.B
    ssf_opt_str.f0 = opt_str.f0
    ssf_opt_str.dt = opt_str.dt

    logging.info('    MPMB covariance estimation (with spectral shift enabled)...')
    nan_mask = np.zeros( (N_rg, N_az, N_imm) , dtype=bool)
    
    for pol_name_curr in pol_names:
        all_acq_data = np.zeros( (N_rg, N_az, N_imm), dtype=complex )
        for acq_idx_curr, acq_curr in enumerate(acq_names):
            all_acq_data[:,:,acq_idx_curr] = D[acq_curr][pol_name_curr]
        nan_mask = nan_mask + np.isnan( all_acq_data )
    
    # this has been moved from the innner cycle (referring to PM9)
    for acq_idx_curr, acq_curr in enumerate( acq_names):
            for pol_name_curr in pol_names:
                D[acq_curr][pol_name_curr][nan_mask[:,:,acq_idx_curr]] = 0
            
    for ch_i, pol_id_i in enumerate( pol_names ):
    
        ind_i = np.arange(N_imm) + ch_i*N_imm
        
        for ch_j in np.arange(ch_i, num_pols):
            pol_id_j = pol_names[ch_j]
            if Mask[ch_i,ch_j] == 1:
                print('    .')
                ind_j = np.arange(N_imm) + ch_j*N_imm
                I = np.zeros((N_imm, N_imm, N_rg_out, N_az_out), dtype=np.complex64)
                for n_idx, acq_id_n in enumerate( acq_names ):
                    if ch_i == ch_j:
                        m_min = n_idx#?
                    else:
                        m_min = 0
                    
                    for m_idx in np.arange(m_min, N_imm):
                        acq_id_m = acq_names[m_idx]
					
                        if cell_data:
                            if double_data:
                                I1 = D[ acq_id_n ][ pol_id_i ]
                                I2 = D[ acq_id_m ][ pol_id_j ]
                            else:
                                I1 =  D[ acq_id_n ][ pol_id_i ]. astype('complex128')
                                I2 = (D[ acq_id_m ][ pol_id_j ]).astype('complex128')
                        else:
                            if double_data:
                                I1 = D[ acq_id_n ][ pol_id_i ]
                                I2 = D[ acq_id_m ][ pol_id_j ]
                            else:
                                I1 =  D[ acq_id_n ][ pol_id_i ]. astype('complex128')
                                I2 = (D[ acq_id_m ][ pol_id_j ]).astype('complex128')
                    
                    # Demodulation to baseband
                    I1 = I1*np.exp(-1j*4*np.pi/opt_str.wavelenght*R[acq_id_n])
                    I2 = I2*np.exp(-1j*4*np.pi/opt_str.wavelenght*R[acq_id_m])

                    # Spectral shift filtering
                    offnadir = np.zeros( ( look_angles[acq_id_m].shape[0],look_angles[acq_id_m].shape[1], 2  ) )
                    offnadir[:,:,0] = look_angles[acq_id_n]
                    offnadir[:,:,1] = look_angles[acq_id_m]
                    [I1, I2] = SpectralShiftFiltering(I1, I2,offnadir, ground_slope, ssf_opt_str)

                    # Modulation
                    I1 = I1*np.exp(1j*4*np.pi/opt_str.wavelenght*R[acq_id_n])
                    I2 = I2*np.exp(1j*4*np.pi/opt_str.wavelenght*R[acq_id_m])
                    
                    temp = Fr_normalized @ (I1*np.conjugate(I2))
                    temp = temp @ Fa_normalized_transposed
                    I[n_idx, m_idx, :, :] = temp.astype('complex64')
                    
                Cov_MPMB[ind_i[:, np.newaxis], ind_j[np.newaxis, :], :, :] = I

    
    # Symmetric part generation
    diag_mask = np.tile((np.eye(num_pols*N_imm).reshape(num_pols*N_imm, num_pols*N_imm, 1, 1))>0, [1, 1, N_rg_out, N_az_out])
    Cov_MPMB = Cov_MPMB + np.conjugate(np.moveaxis(Cov_MPMB, [1, 0, 2, 3], [0, 1, 2, 3]))
    Cov_MPMB[diag_mask] = Cov_MPMB[diag_mask]/2
    
    
    print('    done.\n')
    
    return Cov_MPMB, rg_out, az_out, opt_str.subs_F_r, opt_str.subs_F_a 

    
def SpectralShiftFiltering(I1, I2, look_angles, ground_slope, opt_str):
    """
     
    INPUT
        I1: [Nr x Nc] baseband SLC
        I2: [Nr x Nc] baseband SLC
        look_angles: [Nr x Nc x 2] look angles associated with each SLC
        ground_slope: [Nr x Nc] ground slope in range direction
        opt_str.
              dt: sampling interval
              f0: central frequency
              B: bandwidth
     
    OUTPUT
        I1_filt: [Nr x Nc] baseband spectral shift filtered SLC
        I2_filt: [Nr x Nc] baseband spectral shift filtered SLC"""

    # Spectral shift (Hz)
    Df = opt_str.f0* np.squeeze( np.diff(look_angles, axis=2) )/np.tan(np.mean(look_angles, axis=2) - ground_slope)
    
    # Space-varying demodulation phases (half)
    phi_half = np.pi*np.cumsum(Df*opt_str.dt, axis=0)

    # Demodulation
    I1_filt = I1*np.exp(-1j*phi_half)
    I2_filt = I2*np.exp(1j*phi_half)

    # Building filter
    if not(hasattr(opt_str, 'Nfilter')):
        opt_str.Nfilter = int( np.round(I1.shape[0]/4) )
    
    curr_filter = firwin(numtaps=opt_str.Nfilter, cutoff=opt_str.B*opt_str.dt)
    curr_filter = curr_filter/np.sqrt(np.sum(np.abs(curr_filter)**2))
    curr_filter = curr_filter.reshape(opt_str.Nfilter, 1)
    
    # Filtering
    I1_filt = oaconvolve(I1_filt, curr_filter, mode='same', axes=0)
    I2_filt = oaconvolve(I2_filt, curr_filter, mode='same', axes=0)

    # Back to baseband
    I1_filt = I1_filt*np.exp(1j*phi_half)
    I2_filt = I2_filt*np.exp(-1j*phi_half)

    return I1_filt, I2_filt



def MPMBshuffle(MPMB_correlation,  rg_vec_subs,  az_vec_subs, Npol, N):
    
    Nrg_subs =  rg_vec_subs.size
    Naz_subs =  az_vec_subs.size
    MBMP_correlation = np.zeros((Npol*N, Npol*N, Nrg_subs, Naz_subs), dtype=np.complex64)
    r = MPMB_correlation.shape[0]
    I = np.eye(r)    #Identity matrix
    S = I * 0
    for i in np.arange(N):
        S[i*Npol:Npol+i*Npol, :] = I[i:r:N, :]
        
    for r in np.arange(Nrg_subs):
        for a in np.arange(Naz_subs):
            MBMP_correlation[:, :, r, a] = S @ MPMB_correlation[:, :, r, a] @ S.T
            
    return MBMP_correlation  


def Covariance4D2Correlation4D(Cov_MPMB):
    """It normalizes each element of the 4-D Multi-Polarimetric Multi-Baseline
    % with respect to the corresponding diagonal terms.
    % 
    % INPUT
    %      Cov_MPMB: [Nimm*Npol x Nimm*Npol x Nr x Na] Covariance matrices
    % 
    % OUTPUT
    %       Corr_MPMB: [Nimm*Npol x Nimm*Npol x Nr x Na] normalized covariance
    %                  matrices
    %       varargout{1}: [Nr x Na x Nimm*Npol] diagonal entries of Cov_MPMB"""
    temp, N, Nr, Na = Cov_MPMB.shape
    diag_mask = np.tile((np.eye(N).reshape(N, N, 1, 1))>0, [1, 1, Nr, Na])
    norm_util = Cov_MPMB[diag_mask]
    norm_diag = np.moveaxis(norm_util.reshape((N, Nr, Na)), [1, 2, 0], [0, 1, 2])
    norm_util = 1/np.sqrt(norm_util + np.spacing(1))
    
    norm_util = norm_util.reshape((N, 1, Nr, Na))*norm_util.reshape((1, N, Nr, Na))
    Corr_MPMB = norm_util*Cov_MPMB
    
    return Corr_MPMB, norm_diag


def Covariance2D2Correlation2D(Cov2D):
    """It normalizes each element of the covariance matrix with respect to the
    % corresponding elements on the main diagonal.
    % 
    % INPUT
    %      Cov2D: [N x N] covariance matrix
    % 
    % OUTPUT
    %       Corr2D: [N x N] normalized covariance matrix"""
    temp = np.diag(1/np.sqrt(np.diag(Cov2D)))
    Corr2D = temp@Cov2D@temp
    
    return Corr2D




def build_filtering_matrix(opt_str):
    """ It builds a sparse matrix to carry out one-dimensional moving average.
        It deals with regularly sampled data.
     
        INPUTS:
           opt_str:
                Nin:    number of rows of the matrix to be filtered
                NWmono: one-sided length of the average window
                subs_F: subsampling factor of the filtered signal
     
         OUTPUTS:
            F:     sparse filtering matrix
            xout:  axis of the filtered signal
            Rnorm: normalizing matrix
    """
    
    xin = np.arange(opt_str.Nin)
    Nxin = xin.size
    xout = np.arange(0, opt_str.Nin, opt_str.subs_F)
    Nxout = xout.size
    
    col = np.kron(xout, np.ones((1, 2*opt_str.NWmono + 1))) + np.kron(np.ones((1, Nxout)), np.arange(-opt_str.NWmono, opt_str.NWmono+1))
    row = np.kron(np.arange(Nxout), np.ones((1, 2*opt_str.NWmono + 1)))
    
    ok_mask = (col >= 0) * (col < Nxin)
    Nok_mask = np.sum(ok_mask)
    
    F = csr_matrix((np.ones((1, Nok_mask)).flatten(), (row[ok_mask], col[ok_mask])), shape=(Nxout, Nxin))
    
    Rnorm = csr_matrix((1/np.sum(F.toarray(), axis = 1), (np.arange(Nxout), np.arange(Nxout))))
    
    return F, xout, Rnorm


def covariance_matrix_mat2vec(matrix_in):
    """ covariance matrix has herimtian symmetry, so the elements under the 
        diagonal are not saved in the matrix_out
        matrix_in:  [3x3xNrgxNaz] or [3x3]
        
        matrix_out: [6xNrgxNaz] or [6,]
        
        Elements in matrix in are stored in the optput one with following indexing convention:
        
        Those are the 3x3 matrix_in elements (for each Nrg, Naz) indexes of extraction:
        [ 0  1  2
          /  4  5
          /  /  8 ]
         
        Those are the 6 matrix_out elements (for each Nrg, Naz)  indexes of insertion
        [ 0 1 2 4 5 8 ]
    """
    extraction_indexes = [0,1,2,4,5,8]
    
    shapes = matrix_in.shape
    if len(shapes) == 4:
        Nrg = shapes[2]
        Naz = shapes[3]
        
        matrix_out = np.zeros( ( 6, Nrg, Naz ), dtype=matrix_in.dtype )
        
        for rg_idx in np.arange(Nrg):
            for az_idx in np.arange(Naz):
    
                matrix_out[ :, rg_idx, az_idx] = matrix_in[:,:,rg_idx, az_idx].reshape( (9) )[extraction_indexes]
    else:

        matrix_out = matrix_in.reshape( (9) )[extraction_indexes]
    
    return matrix_out


def covariance_matrix_vec2mat( vector_in ):
    # vector_in has dimension 6
    # matrix_out has dimensions 3x3
    L = len(vector_in.shape)

#        [ 0         1     2
#          cj(1)     3     4
#          cj(2)  cj(4)    5 ]
    if L==3:
        
        _, Nrg, Naz = vector_in.shape
        matrix_out   = np.zeros( ( 3, 3, Nrg, Naz ) , dtype=vector_in.dtype  )
        
        for rg_idx in np.arange(Nrg):
            for az_idx in np.arange(Naz):
                matrix_out[:,:,rg_idx, az_idx] = [ 
                [ vector_in[0, rg_idx, az_idx],          vector_in[1, rg_idx, az_idx],          vector_in[2, rg_idx, az_idx] ], 
                [ np.conj(vector_in[1, rg_idx, az_idx]), vector_in[3, rg_idx, az_idx],          vector_in[4, rg_idx, az_idx] ],
                [ np.conj(vector_in[2, rg_idx, az_idx]), np.conj(vector_in[4, rg_idx, az_idx]), vector_in[5, rg_idx, az_idx] ]
                            ]
    elif L==1:
        
        matrix_out = np.zeros( ( 3, 3 ) , dtype=vector_in.dtype )
        
        matrix_out[:,:] = [ 
        [ vector_in[0],          vector_in[1],          vector_in[2] ], 
        [ np.conj(vector_in[1]), vector_in[3],          vector_in[4] ],
        [ np.conj(vector_in[2]), np.conj(vector_in[4]), vector_in[5] ]
                    ]
    return matrix_out


def SKPD_processing(Cov_MPMB, opt_str):
    """SKPD_processing
    % 
    % Kernel of the Sum of Kronecker Product Decomposition.
    % 
    % INPUT
    %      Cov_MPMB: [N*Npol x N*Npol] covariance matrix
    %       opt_str.
    %                        N: number of images
    %                     Npol: number of polarizations
    %               Nsubspaces: number of desired subspaces
    %                   Nparam: number of parameters sampling the "a" and "b"
    %                           range of physical admissibility
    % 
    % OUTPUT
    %       out_str.
    %               error: [1 x 4] logical vector.
    %                      1. Cov_MPMB has nan of inf values
    %                      2. Singular matrices have nan or inf values
    %                      3. Singular matrices are not full rank
    %                      4. Could not determine "a" or "b" ranges
    %               lambda_SVD: singular values of the SKPD
    %               Ro: tomographic part of the decomposition (normalized wrt
    %                   Ro{k}(1, 1))
    %               Co: polarimetric part of the decomposition
    %               ax: vector gathering the admissible values for parameter a
    %               bx: vector gathering the admissible values for parameter b
    %               Ra: corresponding tomographic covariance matrices (at the
    %                   edges of the possible range of values of "a")
    %               Rb: corresponding tomographic covariance matrices (at the
    %                   edges of the possible range of values of "b")
    %               Ca: corresponding polarimetric covariance matrices (at the
    %                   edges of the possible range of values of "a")
    %               Cb: corresponding polarimetric covariance matrices (at the
    %                   edges of the possible range of values of "b")
    %               full_rank_Ra: true if Ra is full rank
    %               full_rank_Rb: true if Rb is full rank
    %               full_rank_Ca: true if Ca is full rank
    %               full_rank_Cb: true if Cb is full rank
    %               Ra_coherence: tomographic coherence of Ra
    %               Rb_coherence: tomographic coherence of Rb
    %               Ca_coherence: tomographic coherence of Ca
    %               Cb_coherence: tomographic coherence of Cb
    %               Ra_more_coherent_than_Rb: [1 x 1] logical
    %               ax_R_max_coh: [1 x 2] logical. True if the parameter edge
    %                             provides greater coherence
    %               bx_R_max_coh: [1 x 2] logical. True if the parameter edge
    %                             provides greater coherence
    %               R_most_coeherent: most coherent R
    %               Rcoh_thin: Most coherent scattering mechanism, most
    %                          coherent among the range of admissibility
    %               Rcoh_fat: Most coherent scattering mechanism, least
    %                          coherent among the range of admissibility
    %               Rincoh_thin: Least coherent scattering mechanism, most
    %                          coherent among the range of admissibility
    %               Rincoh_fat: Least coherent scattering mechanism, least
    %                          coherent among the range of admissibility
    % 
    % DEPENDENCIES
    %             LowRankKronApproximation
    %             SemiDefinitePositivenessJointRange
    %                             JointDiagonalization"""
    class out_str:
        pass
    out_str.error = np.zeros((4, 1))
    
    if np.any(np.isnan(Cov_MPMB) | np.isinf(Cov_MPMB)):
        out_str.error[0] = 1
        return out_str
    class LRKA_str:
        pass
    LRKA_str.A    = Cov_MPMB
    LRKA_str.rank = opt_str.Nsubspaces
    LRKA_str.Nr1  = opt_str.num_pols
    LRKA_str.Nc1  = opt_str.num_pols
    LRKA_str.Nr2  = opt_str.num_acq
    LRKA_str.Nc2  = opt_str.num_acq
    LRKA_out_str  = LowRankKronApproximation(LRKA_str)
    Co = LRKA_out_str.B
    Ro = LRKA_out_str.C
    lambda_SVD = LRKA_out_str.lambda_SVD
    
    # Demanding power to polarimetric signature
    for k in np.arange(opt_str.Nsubspaces):
        temp = Ro[k][0, 0]
        Ro[k] = Ro[k]/temp
        Co[k] = Co[k]*temp*lambda_SVD[k]
    
    out_str.lambda_SVD = lambda_SVD[np.arange(opt_str.Nsubspaces)]
    out_str.Ro = Ro
    out_str.Co = Co
    
    if np.any(np.isnan(Ro)) or np.any(np.isnan(Co)) | np.any(np.isinf(Ro)) or np.any(np.isinf(Co)):
        out_str.error[1] = 1
        return out_str
    
    if not(np.linalg.matrix_rank(Ro[0]) == opt_str.num_acq) or not(np.linalg.matrix_rank(Ro[1]) == opt_str.num_acq) or\
            not(np.linalg.matrix_rank(Co[0]) == opt_str.num_pols) or not(np.linalg.matrix_rank(Co[1]) == opt_str.num_pols):
        out_str.error[2] = 1
        return out_str
    
    class SDP_opt_str:
        pass
    SDP_opt_str.Ro = Ro
    SDP_opt_str.Co = Co
    SDP_opt_str.N = opt_str.Nparam
    SDP_out_str = SemiDefinitePositivenessJointRange(SDP_opt_str)
    if SDP_out_str.error:
        out_str.error[3] = 1
        return out_str
    ax = SDP_out_str.a
    bx = SDP_out_str.b
    
    out_str.ax = ax
    out_str.bx = bx
    
    out_str.Ra = list([ax[0]*Ro[0] + (1 - ax[0])*Ro[1],\
        ax[-1]*Ro[0] + (1 - ax[-1])*Ro[1]])
    out_str.Rb = list([bx[0]*Ro[0] + (1 - bx[0])*Ro[1],\
        bx[-1]*Ro[0] + (1 - bx[-1])*Ro[1]])
    
    out_str.Ca = list([-(1 - ax[0])*Co[0] + ax[0]*Co[1],\
        -(1 - ax[-1])*Co[0] + ax[-1]*Co[1]])
    out_str.Cb = list([(1 - bx[0])*Co[0] - bx[0]*Co[1],\
        (1 - bx[-1])*Co[0] - bx[-1]*Co[1]])
    
    out_str.full_rank_Ra = [np.linalg.cond(out_str.Ra[0]) < np.linalg.cond(out_str.Ra[1]), False]
    out_str.full_rank_Ra[1] = not(out_str.full_rank_Ra[0])
    out_str.full_rank_Rb = [np.linalg.cond(out_str.Rb[0]) < np.linalg.cond(out_str.Rb[1]), False]
    out_str.full_rank_Rb[1] = not(out_str.full_rank_Rb[0])
    
    out_str.full_rank_Ca = [np.linalg.cond(out_str.Ca[0]) < np.linalg.cond(out_str.Ca[1]), False]
    out_str.full_rank_Ca[1] = not(out_str.full_rank_Ca[0])
    out_str.full_rank_Cb = [np.linalg.cond(out_str.Cb[0]) < np.linalg.cond(out_str.Cb[1]), False]
    out_str.full_rank_Cb[1] = not(out_str.full_rank_Cb[0])
    
    triu_mask_R = np.triu(np.ones((opt_str.num_acq, opt_str.num_acq)), 1)>0
    Ntriu_mask_R = (opt_str.num_acq**2 - opt_str.num_acq)/2
    triu_mask_C = np.triu(np.ones((opt_str.num_pols, opt_str.num_pols)), 1)>0
    Ntriu_mask_C = (opt_str.num_pols**2 - opt_str.num_pols)/2
    out_str.Ra_coherence = np.zeros((2, 1))
    out_str.Rb_coherence = np.zeros((2, 1))
    out_str.Ca_coherence = np.zeros((2, 1))
    out_str.Cb_coherence = np.zeros((2, 1))
    for k in np.arange(2):
        O = np.diag(1/np.sqrt(np.diag(out_str.Ra[k])))
        temp = O@out_str.Ra[k]@O
        out_str.Ra_coherence[k] = np.sum(np.abs(temp[triu_mask_R]))/Ntriu_mask_R
        O = np.diag(1/np.sqrt(np.diag(out_str.Rb[k])))
        temp = O@out_str.Rb[k]@O
        out_str.Rb_coherence[k] = np.sum(np.abs(temp[triu_mask_R]))/Ntriu_mask_R
        O = np.diag(1/np.sqrt(np.diag(out_str.Ca[k])))
        temp = O@out_str.Ca[k]@O
        out_str.Ca_coherence[k] = np.sum(np.abs(temp[triu_mask_C]))/Ntriu_mask_C
        O = np.diag(1/np.sqrt(np.diag(out_str.Cb[k])))
        temp = O@out_str.Cb[k]@O
        out_str.Cb_coherence[k] = np.sum(np.abs(temp[triu_mask_C]))/Ntriu_mask_C
    
    out_str.Ra_more_coherent_than_Rb = np.max(out_str.Ra_coherence) > np.max(out_str.Rb_coherence)
    out_str.ax_R_max_coh = [out_str.Ra_coherence[0] > out_str.Ra_coherence[1], False]
    out_str.ax_R_max_coh[1] = not(out_str.ax_R_max_coh[0])
    out_str.bx_R_max_coh = [out_str.Rb_coherence[0] > out_str.Rb_coherence[1], False]
    out_str.bx_R_max_coh[1] = not(out_str.bx_R_max_coh[0])
    
    out_str.R_most_coherent = (out_str.Ra_more_coherent_than_Rb)*(\
        out_str.Ra[0]*out_str.ax_R_max_coh[0] + out_str.Ra[1]*out_str.ax_R_max_coh[1]) +\
        (1-out_str.Ra_more_coherent_than_Rb)*(\
        out_str.Rb[0]*out_str.bx_R_max_coh[0] + out_str.Rb[1]*out_str.bx_R_max_coh[1])
    
    out_str.Rcoh_thin = (out_str.Ra_more_coherent_than_Rb)*(\
        out_str.Ra[0]*out_str.ax_R_max_coh[0] + out_str.Ra[1]*out_str.ax_R_max_coh[1]) +\
        (1-out_str.Ra_more_coherent_than_Rb)*(\
        out_str.Rb[0]*out_str.bx_R_max_coh[0] + out_str.Rb[1]*out_str.bx_R_max_coh[1])
    out_str.Rcoh_fat = (out_str.Ra_more_coherent_than_Rb)*(\
        out_str.Ra[0]*out_str.ax_R_max_coh[1] + out_str.Ra[1]*out_str.ax_R_max_coh[0]) +\
        (1-out_str.Ra_more_coherent_than_Rb)*(\
        out_str.Rb[0]*out_str.bx_R_max_coh[1] + out_str.Rb[1]*out_str.bx_R_max_coh[0])
    out_str.Rincoh_thin = (1-out_str.Ra_more_coherent_than_Rb)*(\
        out_str.Ra[0]*out_str.ax_R_max_coh[0] + out_str.Ra[1]*out_str.ax_R_max_coh[1]) +\
        (out_str.Ra_more_coherent_than_Rb)*(\
        out_str.Rb[0]*out_str.bx_R_max_coh[0] + out_str.Rb[1]*out_str.bx_R_max_coh[1])
    out_str.Rincoh_fat = (1-out_str.Ra_more_coherent_than_Rb)*(\
        out_str.Ra[0]*out_str.ax_R_max_coh[1] + out_str.Ra[1]*out_str.ax_R_max_coh[0]) +\
        (out_str.Ra_more_coherent_than_Rb)*(\
        out_str.Rb[0]*out_str.bx_R_max_coh[1] + out_str.Rb[1]*out_str.bx_R_max_coh[0])
    
    return out_str


def LowRankKronApproximation(opt_str):
    """It returns the sum of Kronecker product decomposition of the matrix A
    % just as the singular value decomposition. Matrix A is decomposed as:
    % for k = 1:rank
    %     A = A + kron(B{k}, C{k});
    % end
    % 
    % INPUT
    %      opt_str.
    %              A: [Nr x Nc] matrix to be decomposed and approximated
    %              rank: rank of the approximation
    %              Nr1: number of rows of the matrix B
    %              Nc1: number of columns of the matrix B
    %              Nr2: number of rows of the matrix C
    %              Nr2: number of columns of the matrix C
    % 
    % OUTPUT
    %       out_str.
    %               lambda: singular values of the decomposition
    %               B: [rank x 1] cell of [Nr1 x Nc1] matrices
    %               C: [rank x 1] cell of [Nr2 x Nc2] matrices"""
    RA = np.zeros((opt_str.Nr1*opt_str.Nc1, opt_str.Nr2*opt_str.Nc2), dtype=np.complex128)
    for p in np.arange(opt_str.Nc1):
        Ap = np.zeros((opt_str.Nr1, opt_str.Nr2*opt_str.Nc2), dtype=np.complex128)
        ind_p = np.arange(p*opt_str.Nc2, (p+1)*opt_str.Nc2)
        for q in np.arange(opt_str.Nr1):
            ind_q = np.arange(q*opt_str.Nr2, (q+1)*opt_str.Nr2)
            Aqp = opt_str.A[ind_q[:, np.newaxis], ind_p[np.newaxis, :]]
            Ap[q, :] = Aqp.transpose().flatten()
        ind = np.arange(p*opt_str.Nr1, (p+1)*opt_str.Nr1)
        RA[ind, :] = Ap
    
    U, S, VH = np.linalg.svd(RA)
    class out_str:
        pass
    out_str.lambda_SVD = S
    
    out_str.B = list()
    out_str.C = list()
    for k in np.arange(opt_str.rank):
        out_str.B.append(np.reshape(U[:, k], (opt_str.Nr1, opt_str.Nc1), order='F'))
        out_str.C.append(np.reshape(np.transpose(VH)[:, k], (opt_str.Nr2, opt_str.Nc2), order='F'))# warning !
    
    return out_str


def SemiDefinitePositivenessJointRange(opt_str):
    """Given the relationships:
    %               Rg = a*Ro{1} + (1 - a)*Ro{2}
    %               Rv = b*Ro{2} + (1 - b)*Ro{2}
    %               Cg = (1 - b)*C{1} - b*C{2}
    %               Cv = -(1 - a)*C{1} + a*C{2}
    % this routine returns the range of values for parameters "a" and "b" such
    % that matrices Rg, Rv, Cg, Cv are semi definite positive.
    % 
    % INPUT
    %      opt_str.
    %              Ro: [2 x 1] cell of [N x N] Hermitian matrices
    %              Co: [2 x 1] cell of [Npol x Npol] Hermitian matrices
    %               N: number of points sampling the admissible range of values
    % 
    % OUTPUT
    %       out_str.
    %               a: [1 x opt_str.N] vector
    %               b: [1 x opt_str.N] vector
    %               error: true if an error has occurred
    % 
    % DEPENDENCIES
    %             JointDiagonalization"""
    class out_str:
        pass
    out_str.error = False
    
    # Semi positive definiteness of the interferometric covariance matrices
    U = JointDiagonalization(opt_str.Ro[0:2])
    d1 = np.real(np.diag(U.conj().transpose() @ opt_str.Ro[0] @ U))
    d2 = np.real(np.diag(U.conj().transpose() @ opt_str.Ro[1] @ U))
    
    dd = d1 - d2
    indm = dd < 0
    indp = dd > 0
    temp = -d2/dd
    
    try:
        ab_lo = np.max(temp[indp])
    except ValueError:
        ab_lo = -np.inf
    try:
        ab_up = np.min(temp[indm])
    except ValueError:
        ab_up = np.inf
    # Semi positive definiteness of the polarimetric covariance matrices
    U = JointDiagonalization(opt_str.Co[0:2])
    d1 = np.real(np.diag(U.conj().transpose() @ opt_str.Co[0] @ U))
    d2 = np.real(np.diag(U.conj().transpose() @ opt_str.Co[1] @ U))
    
    dd = d1 + d2
    indm = dd < 0
    indp = dd > 0
    temp = d1/dd
    
    try:
        a_lo = np.max(temp[indp])
    except ValueError:
        a_lo = -np.inf
    try:
        a_up = np.min(temp[indm])
    except ValueError:
        a_up = np.inf
    
    dd = -d1-d2
    indm = dd < 0
    indp = dd > 0
    temp = -d1/dd
    
    try:
        b_lo = np.max(temp[indp])
    except ValueError:
        b_lo = -np.inf
    try:
        b_up = np.min(temp[indm])
    except ValueError:
        b_up = np.inf
    
    if a_lo < b_up:
        out_str.error = True
    # Semi positive definiteness of both covariance matrices
    a_lo = np.max((a_lo, ab_lo))
    a_up = np.min((a_up, ab_up))
    b_lo = np.max((b_lo, ab_lo))
    b_up = np.min((b_up, ab_up))
    
    Da = a_up - a_lo
    Db = b_up - b_lo
    if Da < 0 or Db < 0:
        out_str.error = True
    
    out_str.a = np.linspace(a_lo, a_up, opt_str.N)
    out_str.b = np.linspace(b_lo, b_up, opt_str.N)
    
    return out_str


def JointDiagonalization(C):
    """It finds the matrix U that diagonalizes C{1} and C{2} at the same time
    % (it means that U'*C{1}*U and U'*C{2}*U are diagonal)
    % 
    % INPUT
    %      C: [2 x 1] cell of [N x N] Hermitian matrices
    % 
    % OUTPUT
    %       U: diagonalizing matrix"""
    A, V = np.linalg.eig(C[1] @ np.linalg.inv(C[0]))
    U = np.linalg.inv(V.conj().transpose())
    
    return U