# Constants

# maimum value of the earth rasius (at equator)
MAX_EARTH_RADIUS_METERS = 6378000

# Earth centred, earth fixed, righthanded 3D coordinate system,
# consisting of 3 orthogonal axes with X and Y axes in the equatorial plane,
# positive Z-axis parallel to mean earth rotation axis and pointing towards
# North Pole. UoM: m.
EPSG_CODE_ECEF = 'EPSG:4978'

# Ellipsoidal 2D CS. Axes: latitude, longitude. Orientations: north, east. UoM: degree
EPSG_CODE_LLA = 'EPSG:4326'

# Oversampling factor to be aplied to input data and auxiliary data when
# resolution / pixel_spacing is < 2
OVERSAMPLING_FACTOR = 2

# Range bandwidth of the Biomass data
RANGE_BANDWIDTH_HZ = 6000000

CARRIER_FREQUENCY_HZ = 435000000

# Velocity of the satellite in [m/s]
SATELLITE_VELOCITY = 7.595823208681849e03
