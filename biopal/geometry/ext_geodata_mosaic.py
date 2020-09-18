class ext_geodata_mosaic:
    def __init__(self, latmin, latmax, lonmin, lonmax, database_dir, data_type, output_dir, geoid_dir=None):

        raise NotImplementedError 

    def generate_csv_with_coordinates(self):
        raise NotImplementedError     
        
    def corr_geoid(self,lat,lon):
        raise NotImplementedError   
        
    def get_mosaic(self):
        raise NotImplementedError   
        
    def save_mosaic(self,plot=False,forcePermission=True):
        raise NotImplementedError    
 
    def save_biodempp_metadata(self):
        raise NotImplementedError

    def Run(self):
        raise NotImplementedError