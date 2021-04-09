from biopal.ground_cancellation.ground_cancellation import ground_cancellation
import numpy as np


def initial_values():
    Nrg = 378
    Naz = 2656
    data_dummy = 1/np.sqrt(2)*np.random.rand(Nrg,Naz)+ 1j/np.sqrt(2)*np.random.rand(Nrg,Naz)
    pol_dict = {'hh': data_dummy, 'hv': data_dummy, 'vv': data_dummy }
    data_input_dict = {'BSL_00': pol_dict, 'BSL_01': pol_dict, 'BSL_02': pol_dict }
    
    wave_number_dummy = np.ones((Nrg,Naz))
    wave_number_dict = {'BSL_00': 0*wave_number_dummy, 'BSL_01': -0.03*wave_number_dummy, 'BSL_02': 0.08*wave_number_dummy }
    
    multi_master_flag = True
    enhanced_forest_height = 50
    equalization_flag = '1'
    resolution_m_slant_rg = 25
    off_nadir_angle = 0.5*np.ones((Nrg,Naz))
    slope = 0.05*np.ones((Nrg,Naz))
    
    return (
        Nrg, 
        Naz,            
        data_input_dict,
        wave_number_dict,
        multi_master_flag,
        enhanced_forest_height,
        equalization_flag,
        resolution_m_slant_rg,
        off_nadir_angle,
        slope,
        )
    
    
def test_ground_cancellation_1():
    
    (
     Nrg, 
     Naz,            
     data_input_dict,
     wave_number_dict,
     multi_master_flag,
     enhanced_forest_height,
     equalization_flag,
     resolution_m_slant_rg,
     off_nadir_angle,
     slope,
     ) = initial_values()
    
    notched_pol_dict = ground_cancellation(
        data_input_dict,
        wave_number_dict,
        multi_master_flag,
        enhanced_forest_height,
        equalization_flag,
        resolution_m_slant_rg,
        off_nadir_angle,
        slope,
    )
    
    pols = notched_pol_dict.keys()
    assert ('hh' in pols and 'hv' in pols and 'vv' in pols) and len(pols)==3, 'wrong output number of polarizations, in test with space varying OFF'
        
    for pol, data in notched_pol_dict.items():
        
        assert data.shape==(Nrg,Naz), 'wrong output shape, in test with space varying OFF'


def test_ground_cancellation_2():
    
    (
      Nrg, 
      Naz,            
      data_input_dict,
      wave_number_dict,
      multi_master_flag,
      enhanced_forest_height_scalar,
      equalization_flag,
      resolution_m_slant_rg,
      off_nadir_angle,
      slope,
      ) = initial_values()
    
    enhanced_forest_height_map = slope.copy()
    enhanced_forest_height_map.fill(enhanced_forest_height_scalar)

           
    notched_pol_dict = ground_cancellation(
        data_input_dict,
        wave_number_dict,
        multi_master_flag,
        enhanced_forest_height_map,
        equalization_flag,
        resolution_m_slant_rg,
        off_nadir_angle,
        slope,
    )
    
    pols = notched_pol_dict.keys()
    assert ('hh' in pols and 'hv' in pols and 'vv' in pols) and len(pols)==3, 'wrong output number of polarizations, in test with space varying ON'
        
    for pol, data in notched_pol_dict.items():
        
        assert data.shape==(Nrg,Naz), 'wrong output shape, in test with space varying ON'