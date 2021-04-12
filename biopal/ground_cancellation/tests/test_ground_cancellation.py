from biopal.ground_cancellation.ground_cancellation import ground_cancellation
import numpy as np
import warnings


class TestGroundCancellation():
    """
    "ground_cancellation" APP functional tests
    performed with pytest framework: https://docs.pytest.org/en/stable/getting-started.html
    
    "ground_cancellation" APP is executed twice:
        
        first time with space varying feature disabled and dummy inputs.
        
        second time with space varying feature enabled and dummy inputs: 
               to enable the space varying feature, the scalar enhanced forest height 
               coming from the "initial_values" function is converted to a 2D map.
        
    The generated output is tested to be a dictionary with the three expected 
    polarizations, each containing a 2D map with expected Nrg, Naz dimensions 
    and with expected data type.
   
    Test execution:
    "pytest" needs to be installed in the environment:
        
        >>>pip install -U pytest
        
    From command line go in the "BioPAL/" directory and type the following instruction:
        
        >>>pytest
        
    This will automatically execute all the pytest-formatted tests in BioPAL
    
    Parameters
    ----------
    None, all the dummy parameters are loaded from the "initial_values" function
    
    """

    def test_ground_cancellation(self):
        """
        ground_cancellation functional test #1
        Executes ground_cancellation APP, with space varying feature disabled and 
        dummy inputs
        """
        
        (
         Nrg, 
         Naz,            
         data_stack_dict,
         wave_number_dict,
         multi_master_flag,
         enhanced_forest_height,
         equalization_flag,
         resolution_m_slant_rg,
         off_nadir_angle,
         slope,
         ) = self.initial_values()
        
        notched_pol_dict = ground_cancellation(
            data_stack_dict,
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
            
            if multi_master_flag:
                assert isinstance(data[0,0], np.float64), 'wrong output data type; when multi_master_flag is True it should be float64, in test with space varying OFF'
            else:
                assert isinstance(data[0,0], np.complex128), 'wrong output data type; when multi_master_flag is False it is supposed to be complex128, in test with space varying OFF'
    
    
    def test_ground_cancellation_space_varying(self):
        """
        ground_cancellation functional test #2
        Executes ground_cancellation APP, with space varying feature enabled and 
        dummy inputs
        """
        
        (
          Nrg, 
          Naz,            
          data_stack_dict,
          wave_number_dict,
          multi_master_flag,
          enhanced_forest_height_scalar,
          equalization_flag,
          resolution_m_slant_rg,
          off_nadir_angle,
          slope,
          ) = self.initial_values()
        
        # converting to 2D map, in order to enable the space varying feature
        enhanced_forest_height_map = slope.copy()
        enhanced_forest_height_map.fill(enhanced_forest_height_scalar)
        
        warnings.filterwarnings('ignore') # warnings silenced for convenience
        notched_pol_dict = ground_cancellation(
            data_stack_dict,
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
            
            if multi_master_flag:
                assert isinstance(data[0,0], np.float64), 'wrong output data type; when multi_master_flag is True it should be float64, in test with space varying ON'
            else:
                assert isinstance(data[0,0], np.complex128), 'wrong output data type; when multi_master_flag is False it is supposed to be complex128, in test with space varying ON'
    
    
    def initial_values(self):
        """
        Generates dummy variables for the ground cancellation APP, for functional testing
        
        Returns
        -------
        Nrg
            integer, number of range samples for the dummy maps
        Naz
            integer, number of azimuth lines for the dummy maps
        data_stack_dict
            dictionary with three acquisitions, each containing a sub dictionary with 
            three polarizations, each containing  a raster with random complex values
        wave_number_dict
            dictionary with three acquisitions, each containing a raster with 
            constant dummy real values
        multi_master_flag
            boolean (see ground cancellation documentation)
        enhanced_forest_height
                positive number (see ground cancellation documentation)
        equalization_flag
            string, enumeration flag (see ground cancellation documentation)
        resolution_m_slant_rg
           positive number, slant range resolutuon 
        off_nadir_angle
            map with dummy constant values for off nadir angle (supposed radiants)
        slope
            map with dummy constant values for slope angle (supposed radiants)
        """
        
        Nrg = 378
        Naz = 2656
        data_dummy = 1/np.sqrt(2)*np.random.rand(Nrg,Naz)+ 1j/np.sqrt(2)*np.random.rand(Nrg,Naz)
        pol_dict = {'hh': data_dummy, 'hv': data_dummy, 'vv': data_dummy }
        data_stack_dict = {'BSL_00': pol_dict, 'BSL_01': pol_dict, 'BSL_02': pol_dict }
        
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
            data_stack_dict,
            wave_number_dict,
            multi_master_flag,
            enhanced_forest_height,
            equalization_flag,
            resolution_m_slant_rg,
            off_nadir_angle,
            slope,
            )        
        