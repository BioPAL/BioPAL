# SPDX-FileCopyrightText: Aresys S.r.l. <info@aresys.it>
# SPDX-License-Identifier: MIT

import os
import enum
import numpy as np
from scipy import constants
from progressbar import progressbar
from matplotlib import pyplot as plt
from arepytools.io.metadata import StateVectors, RasterInfo
from arepytools.io.productfolder import ProductFolder
from arepytools.timing import precisedatetime
from arepytools.geometry.generalsarorbit import create_general_sar_orbit
from arepytools.geometry import conversions
from arepytools import constants as cst

from biopal.io.xml_io import XmlIO
from biopal.geometry.ext_geodata_mosaic import ext_geodata_mosaic

_NEWTON_TOLERANCE = 1.0e-8
_MAX_ITERATIONS = 8


class EProductType(enum.Enum):
    """
    Products types
    """

    ecef_grid = "ECEFGRID"
    reference_height = "REFERENCEHEIGHT"
    distance = "DISTANCE"
    off_nadir = "OFFNADIR"
    slope = "SLOPE"
    wave_number = "WAVENUMBER"


def shift_image(data, shift_x, shift_y):
    x_frequency_axis = np.fft.fftfreq(data.shape[1])[np.newaxis, :]
    y_frequency_axis = np.fft.fftfreq(data.shape[0])[:, np.newaxis]

    data_shifted = np.fft.ifft(np.fft.fft(data, axis=1) * np.exp(-1j * 2 * np.pi * x_frequency_axis * shift_x), axis=1)
    data_shifted = np.fft.ifft(np.fft.fft(data_shifted, axis=0) * np.exp(-1j * 2 * np.pi * y_frequency_axis * shift_y), axis=0)
    return data_shifted.astype(data.dtype)


def get_zero_doppler_plane(velocity_vector, position_vector):
    a = velocity_vector[0]
    b = velocity_vector[1]
    c = velocity_vector[2]
    d = -(velocity_vector.T @ position_vector).flatten()
    return np.concatenate([a, b, c, d])


def get_distance_from_zero_doppler(plane_coefficients, target_coords):
    a = plane_coefficients[0]
    b = plane_coefficients[1]
    c = plane_coefficients[2]
    d = plane_coefficients[3]
    distance = a * target_coords[0, :] + b * target_coords[1, :] + c * target_coords[2, :] + d
    distance /= np.sqrt(a ** 2 + b ** 2 + c ** 2)
    return distance


def lerp(x_in, y_in):
    return (y_in[0] * x_in[1] - y_in[1] * x_in[0]) / (x_in[1] - x_in[0])


def get_min_indexes(vector):
    index_min_dist = np.zeros(2)
    index_min_dist[0] = np.where(np.diff(vector < 0))[0][0]
    index_min_dist[1] = index_min_dist[0] + 1
    return index_min_dist


def cut_state_vectors(sv_in: StateVectors, azimuth_start, azimuth_stop, guard_sv=5):
    first_sv = np.max([int(np.round((azimuth_start - sv_in.reference_time) / sv_in.time_step)) - guard_sv, 0])
    last_sv = np.min(
        [
            int(np.round((azimuth_stop - sv_in.reference_time) / sv_in.time_step)) + guard_sv,
            sv_in.number_of_state_vectors,
        ]
    )

    t_ref_utc = sv_in.reference_time + first_sv * sv_in.time_step
    position_vector = sv_in.position_vector[first_sv : (last_sv + 1)]
    velocity_vector = sv_in.velocity_vector[first_sv : (last_sv + 1)]
    dt_sv_s = sv_in.time_step

    return StateVectors(position_vector, velocity_vector, t_ref_utc, dt_sv_s)


def compute_interferometric_baseline(
    state_vectors_master,
    state_vectors_slave,
    azimuth_time_master,
    range_axis_master,
    look_direction,
    geodetic_height=0.0,
):
    gso_master = create_general_sar_orbit(state_vectors_master)
    gso_slave = create_general_sar_orbit(state_vectors_slave)

    # Direct geocoding with the master orbit
    targets_coords = gso_master.sat2earth(azimuth_time_master, range_axis_master, look_direction, geodetic_height)

    # Inverse geocoding with the slave orbit
    azimuth_time_slave, _ = gso_slave.earth2sat(targets_coords[:, 0])
    sensor_slave_coords = gso_slave.get_position(azimuth_time_slave)

    # Computation of the baseline vector
    sensor_master_coords = gso_master.get_position(azimuth_time_master)
    baseline_vector = sensor_slave_coords - sensor_master_coords

    # Projection of the baseline vector to the SAR directions
    los_versors = sensor_master_coords - targets_coords
    los_versors = los_versors / np.linalg.norm(los_versors, axis=0)
    velocity_versor = gso_slave.get_velocity(azimuth_time_master)
    velocity_versor = velocity_versor / np.linalg.norm(velocity_versor)
    norm_versor = np.cross(los_versors, velocity_versor, axisa=0, axisb=0).T

    normal_baseline = baseline_vector.T @ norm_versor
    parallel_baseline = baseline_vector.T @ los_versors
    along_track_baseline = baseline_vector.T @ velocity_versor

    return normal_baseline, parallel_baseline, along_track_baseline


def compute_vertical_wavenumber(normal_baseline, wavelength, distance, angles):
    return 4 * np.pi * normal_baseline / (wavelength * distance * np.sin(angles))


def compute_vertical_wavenumber_angles(wavelength, off_nadir_angles_master, off_nadir_angles_slave):
    return (
        4 * np.pi * (off_nadir_angles_slave - off_nadir_angles_master) / (wavelength * np.sin(off_nadir_angles_master))
    )


def save_product(data_matrix_list, azimuth_axis, range_axis, output_path, product_type: EProductType):
    import shutil

    if not os.path.exists(output_path):
        os.mkdir(output_path)
    path_out_product = os.path.join(output_path, "geometric_product_" + product_type.value.lower())
    if os.path.exists(path_out_product):
        shutil.rmtree(path_out_product)
    pf = ProductFolder(path_out_product, "w")

    for channel_index, data_matrix in enumerate(data_matrix_list):
        pf.append_channel(len(azimuth_axis), len(range_axis), "FLOAT64")
        ch = pf.get_channel(channel_index)
        md = ch.metadata
        md_ch = md.get_metadata_channels(0)
        ri = md_ch.get_element("RasterInfo")
        ri.set_lines_axis(azimuth_axis[0], "Utc", np.diff(azimuth_axis[:2]), "s")
        ri.set_samples_axis(range_axis[0], "s", np.diff(range_axis[:2]), "s")
        pf.write_data(channel_index=channel_index, data=data_matrix, start_point=(0, 0))
        pf.write_metadata(channel_index=channel_index)


def cut_raster_info(raster_info_input: RasterInfo, roi_pixels) -> RasterInfo:
    """This function returns a raster info portion according to the roi in pixel in input

    Arguments:
        raster_info_input {RasterInfo} -- [description]
        roi_pixels {[type]} -- [description]

    Returns:
        RasterInfo -- [description]
    """
    lines_start = raster_info_input.lines_start + roi_pixels[0] * raster_info_input.lines_step
    samples_start = raster_info_input.samples_start + roi_pixels[1] * raster_info_input.samples_step

    raster_info_out = RasterInfo(roi_pixels[2], roi_pixels[3], raster_info_input.cell_type.value)
    raster_info_out.set_lines_axis(
        lines_start,
        raster_info_input.lines_start_unit,
        raster_info_input.lines_step,
        raster_info_input.lines_step_unit,
    )
    raster_info_out.set_samples_axis(
        samples_start,
        raster_info_input.samples_start_unit,
        raster_info_input.samples_step,
        raster_info_input.samples_step_unit,
    )
    return raster_info_out


def get_geo_corners(
    raster_info: RasterInfo, state_vectors: StateVectors, side_looking, max_height=10.0e3, min_height=-1.0e3,
):
    """This function returns the maximum extension of the lat/lon corners for the raster info
    and the orbit in input. The range of altitude used to complute the geographic extension is
    from -1 to 10 Km.

    Arguments:
        raster_info {RasterInfo} -- [description]
        state_vectors {StateVectors} -- [description]
        side_looking {[type]} -- [description]

    Keyword Arguments:
        max_height {[type]} -- [description] (default: {10.0e3})
        min_height {[type]} -- [description] (default: {-1.0e3})

    Returns:
        tuple -- (lat_min, lat_max, lon_min, lon_max) [rad]
    """
    gso = create_general_sar_orbit(state_vectors)

    # get data corners
    lat_min = np.nan
    lon_min = np.nan
    lat_max = np.nan
    lon_max = np.nan
    lines_axis = [
        raster_info.lines_start,
        raster_info.lines_start + (raster_info.lines - 1) * raster_info.lines_step,
    ]
    samples_axis = [
        raster_info.samples_start,
        raster_info.samples_start + (raster_info.samples - 1) * raster_info.samples_step,
    ]

    for height in (max_height, min_height):
        for az in lines_axis:
            for rg in samples_axis:
                coord = conversions.xyz2llh(gso.sat2earth(az, rg, side_looking, height))
                lat_min = np.nanmin([lat_min, coord[0]])
                lat_max = np.nanmax([lat_max, coord[0]])
                lon_min = np.nanmin([lon_min, coord[1]])
                lon_max = np.nanmax([lon_max, coord[1]])

    return lat_min, lat_max, lon_min, lon_max


def check_geometric_corners(
    raster_info: RasterInfo, state_vector: StateVectors, side_looking, geometric_corners, tolerance=1e-6,
):
    geometric_corners = np.array(geometric_corners) // tolerance * tolerance
    raster_geocorners = get_geo_corners(raster_info, state_vector, side_looking)
    raster_geocorners = np.array(raster_geocorners) // tolerance * tolerance
    if geometric_corners[0] > raster_geocorners[0]:
        return False
    if geometric_corners[1] < raster_geocorners[1]:
        return False
    if geometric_corners[2] > raster_geocorners[2]:
        return False
    if geometric_corners[3] < raster_geocorners[3]:
        return False
    return True


def check_roi(roi, ri):
    if (roi[0] < 0) or (roi[1] < 0) or (roi[0] + roi[2] > ri.lines) or (roi[1] + roi[3] > ri.samples):
        raise Exception(f"The ROI {roi} is not correct")
    return roi


def get_axis_from_roi(ri, roi=None):
    check_roi(roi, ri)
    if roi is None:
        roi = [0, 0, ri.lines, ri.samples]

    azimuth_axis_rel = np.arange(roi[0], roi[0] + roi[2]) * ri.lines_step
    range_axis = np.arange(roi[1], roi[1] + roi[3]) * ri.samples_step + ri.samples_start

    return azimuth_axis_rel, range_axis


def _initialize_inverse_geocoding(sensor_pos, sensor_vel, point, guess):
    los = point - sensor_pos[:, guess]
    previous_doppler = np.dot(los, sensor_vel[:, guess])
    num_sat_pos = sensor_pos.shape[1]
    for k in range(guess, num_sat_pos):
        los = point - sensor_pos[:, k]
        current_doppler = np.dot(los, sensor_vel[:, k])
        if current_doppler * previous_doppler <= 0 and current_doppler < previous_doppler:
            return max(int(k - 0.5), 0)
        previous_doppler = current_doppler

    raise RuntimeError("Cannot find initialization point")


def perform_inverse_geocoding(gso, points):
    num_points = points.shape[1]

    # Output preallocation
    az = np.full((num_points,), precisedatetime.PreciseDateTime())
    rg = np.zeros((num_points,))

    # Orbit data extraction
    t0, dt, num_sat_pos = gso.t0, gso.dt, gso.n
    sensor_pos = gso._state_vectors
    sensor_vel = gso.get_velocity(gso._time_axis.get_array())
    interpolator = gso.interpolator
    duration = interpolator._POLYNOMIAL_ORDER * dt

    # Orbit local polynomial coefficients
    coeff_pos = np.zeros((3, interpolator._POLYNOMIAL_ORDER + 1))
    coeff_vel = np.zeros((3, interpolator._POLYNOMIAL_ORDER + 1))
    coeff_acc = np.zeros((3, interpolator._POLYNOMIAL_ORDER + 1))

    sat_position = np.zeros((3,))
    sat_velocity = np.zeros((3,))
    sat_acceleration = np.zeros((3,))

    interval_index_guess = 0
    for point_index in progressbar(range(num_points), max_value=num_points):

        interval_index = _initialize_inverse_geocoding(
            sensor_pos, sensor_vel, points[:, point_index], interval_index_guess
        )
        polynomial_index = max(interval_index - interpolator._POLY_CENTER_REL_INTERVAL_INDEX, 0)

        for coord in range(3):
            coeff_pos[coord, :5] = interpolator._pos[polynomial_index, coord].c
            coeff_vel[coord, :4] = interpolator._vel[polynomial_index, coord].c
            coeff_acc[coord, :3] = interpolator._acc[polynomial_index, coord].c

        counter = 0
        time_old = 0
        t = duration / 2

        while counter < _MAX_ITERATIONS and abs(t - time_old) > _NEWTON_TOLERANCE:
            if t < 0 or t > duration:
                raise RuntimeError("Cannot perform inverse geocoding")

            t2 = t * t
            t3 = t2 * t
            t4 = t3 * t

            for coord in range(3):
                sat_position[coord] = (
                    coeff_pos[coord, 0] * t4
                    + coeff_pos[coord, 1] * t3
                    + coeff_pos[coord, 2] * t2
                    + coeff_pos[coord, 3] * t
                    + coeff_pos[coord, 4]
                )
                sat_velocity[coord] = (
                    coeff_vel[coord, 0] * t3 + coeff_vel[coord, 1] * t2 + coeff_vel[coord, 2] * t + coeff_vel[coord, 3]
                )
                sat_acceleration[coord] = coeff_acc[coord, 0] * t2 + coeff_acc[coord, 1] * t + coeff_acc[coord, 2]

            sat2point = sat_position - points[:, point_index]

            f = np.dot(sat2point, sat_velocity)
            df = np.dot(sat_velocity, sat_velocity) + np.dot(sat_acceleration, sat2point)

            time_old = t
            t -= f / df

            counter = counter + 1

        if counter == 8:
            raise RuntimeError("Exceeded max number of iterations during inverse geocoding")

        slant_range = np.linalg.norm(sat2point)

        rg[point_index] = slant_range * 2 / cst.LIGHT_SPEED
        az[point_index] = interpolator._data_axis_array[polynomial_index] + t

        interval_index_guess = max(int((az[point_index] - t0) / dt), 0)

    return az, rg


def is_north_heading(velocity):
    velocity_llh = conversions.xyz2llh(velocity)
    index_max = np.argmax(np.abs(velocity_llh[:2]))
    return index_max == 0


class SARGeometry:
    def __init__(
        self, raster_info: RasterInfo, state_vectors: StateVectors, lat_axis, lon_axis, dem, flag_check_dem=True,
    ):
        self._raster_info = raster_info
        self._state_vector = cut_state_vectors(
            state_vectors,
            self._raster_info.lines_start,
            self._raster_info.lines_start + self._raster_info.lines * self._raster_info.lines_step,
        )
        self._general_sar_orbit = create_general_sar_orbit(self._state_vector)
        self._lat_axis = lat_axis
        self._lon_axis = lon_axis
        self._dem = dem
        self._roi = [0, 0, self._raster_info.lines, self._raster_info.samples]

        # Check the DEM validity
        if flag_check_dem:
            geometric_corners = (
                np.min(lat_axis),
                np.max(lat_axis),
                np.min(lon_axis),
                np.max(lon_axis),
            )
            if not (
                check_geometric_corners(self._raster_info, self._state_vector, "RIGHT", geometric_corners)
                or check_geometric_corners(self._raster_info, self._state_vector, "LEFT", geometric_corners)
            ):
                raise Exception("the DEM in input is not compliant with the ROI of SAR data.")

        self._dem_sar = None
        self._x_sar_coords = None
        self._y_sar_coords = None
        self._z_sar_coords = None
        self._terrain_slope = None

        self.flag_display = False

    @property
    def roi(self):
        return self._roi

    @property
    def lat_axis(self):
        return self._lat_axis

    @property
    def lon_axis(self):
        return self._lon_axis

    @property
    def dem(self):
        return self._dem

    @property
    def dem_sar(self):
        return self._dem_sar

    @property
    def x_sar_coords(self):
        return self._x_sar_coords

    @property
    def y_sar_coords(self):
        return self._y_sar_coords

    @property
    def z_sar_coords(self):
        return self._z_sar_coords

    @property
    def terrain_slope(self):
        return self._terrain_slope

    def get_axis(self):
        azimuth_axis, range_axis = get_axis_from_roi(self._raster_info, self.roi)
        azimuth_axis = azimuth_axis + self._raster_info.lines_start
        return azimuth_axis, range_axis

    def get_targets_coordinates(self, lon_index=None, lat_index=None):
        lon_matrix, lat_matrix = np.meshgrid(self.lon_axis[lon_index], self.lat_axis[lat_index])

        lat_indexes = np.arange(self.lat_axis.size)
        lon_indexes = np.arange(self.lon_axis.size)
        targets_coord = np.concatenate(
            [
                lat_matrix.flatten()[np.newaxis, :],
                lon_matrix.flatten()[np.newaxis, :],
                self.dem[lat_indexes[lat_index].T, lon_indexes[lon_index]].flatten()[np.newaxis, :],
            ]
        )
        targets_coord = conversions.llh2xyz(targets_coord)
        return targets_coord

    def get_indexes_dem_portion(
        self, zero_doppler_plane, targets_coord_first, targets_coord_last, flag_north_heading=True
    ):
        lat_axis = np.arange(self.lat_axis.size, dtype=np.int)
        lon_axis = np.arange(self.lon_axis.size, dtype=np.int)
        indexes = np.zeros(2, dtype=np.int)
        indexes[0] = np.nanargmin(
            np.abs(get_distance_from_zero_doppler(zero_doppler_plane.astype(np.float32), targets_coord_first))
        )
        indexes[1] = np.nanargmin(
            np.abs(get_distance_from_zero_doppler(zero_doppler_plane.astype(np.float32), targets_coord_last))
        )
        if flag_north_heading:
            lat_axis = np.arange(
                np.max([np.min(indexes) - 10, 0]), np.min([np.max(indexes) + 10, self.lat_axis.size - 1]), dtype=np.int,
            )
            lon_matrix_current, lat_matrix_current = np.meshgrid(np.arange(self.lon_axis.size, dtype=np.int), lat_axis)
            current_indexes = np.ravel_multi_index(
                np.concatenate(
                    [lat_matrix_current.flatten()[np.newaxis, :], lon_matrix_current.flatten()[np.newaxis, :],]
                ),
                [targets_coord_first.shape[1], self.lon_axis.size],
            )
        else:
            lon_axis = np.arange(
                np.max([np.min(indexes) - 10, 0]), np.min([np.max(indexes) + 10, self.lon_axis.size - 1]), dtype=np.int,
            )
            lon_matrix_current, lat_matrix_current = np.meshgrid(lon_axis, np.arange(self.lat_axis.size, dtype=np.int))
            current_indexes = np.ravel_multi_index(
                np.concatenate(
                    [lat_matrix_current.flatten()[np.newaxis, :], lon_matrix_current.flatten()[np.newaxis, :],]
                ),
                [self.lat_axis.size, targets_coord_first.shape[1]],
            )

        return lat_axis, lon_axis, current_indexes

    def compute_distance_master(self):
        sensor_positions = self.get_sensor_positions_master()
        distances = np.sqrt(
            (sensor_positions[0, :] - self.x_sar_coords) ** 2
            + (sensor_positions[1, :] - self.y_sar_coords) ** 2
            + (sensor_positions[2, :] - self.z_sar_coords) ** 2
        )
        return distances

    def compute_distance_slave(self, state_vector_ext):
        gso = create_general_sar_orbit(state_vector_ext)

        targets_coords = np.vstack(
            (self.x_sar_coords.flatten(), self.y_sar_coords.flatten(), self.z_sar_coords.flatten())
        )

        _, range_points_master_slave = perform_inverse_geocoding(gso, targets_coords)
        range_points_master_slave.shape = self.x_sar_coords.shape
        distances = range_points_master_slave * constants.c / 2

        return distances

    def compute_distance(self, state_vector_ext: StateVectors = None):
        if state_vector_ext is None:
            return self.compute_distance_master()
        else:
            return self.compute_distance_slave(state_vector_ext)

    def get_sensor_positions_master(self):
        azimuth_time_axis, _ = self.get_axis()
        gso = self._general_sar_orbit
        sensor_position = gso.get_position(azimuth_time_axis)
        sensor_position = np.repeat(sensor_position[:, :, np.newaxis], tuple(self.x_sar_coords.shape)[1], axis=2)

        return sensor_position

    def get_sensor_positions_slave(self, state_vector_ext):
        azimuth_time_axis, _ = self.get_axis()

        gso = create_general_sar_orbit(state_vector_ext)

        targets_coords = np.vstack(
            (self.x_sar_coords.flatten(), self.y_sar_coords.flatten(), self.z_sar_coords.flatten())
        )

        azimuth_points_master_slave, _ = perform_inverse_geocoding(gso, targets_coords)

        sensor_position = gso.get_position(azimuth_points_master_slave)

        sensor_position = np.reshape(
            sensor_position, [3, tuple(self.x_sar_coords.shape)[0], tuple(self.x_sar_coords.shape)[1]],
        )

        return sensor_position

    def get_sensor_positions(self, state_vector_ext: StateVectors = None):
        azimuth_time_axis, _ = self.get_axis()
        if state_vector_ext is None:
            return self.get_sensor_positions_master()
        else:
            return self.get_sensor_positions_slave(state_vector_ext)

    def get_sensor_velocity_master(self):
        azimuth_time_axis, _ = self.get_axis()
        gso = self._general_sar_orbit
        sensor_velocity = gso.get_velocity(azimuth_time_axis)
        sensor_velocity = np.repeat(sensor_velocity[:, :, np.newaxis], tuple(self.x_sar_coords.shape)[1], axis=2)

        return sensor_velocity

    def get_sensor_velocity_slave(self, state_vector_ext):
        gso = create_general_sar_orbit(state_vector_ext)
        sensor_velocity = np.zeros([3, self.x_sar_coords.size])
        targets_coords = np.concatenate(
            [
                self.x_sar_coords.flatten()[:, np.newaxis],
                self.y_sar_coords.flatten()[:, np.newaxis],
                self.z_sar_coords.flatten()[:, np.newaxis],
            ],
            axis=1,
        )
        # Inverse geocoding
        for index_target, target_coord in progressbar(enumerate(targets_coords), max_value=self.x_sar_coords.size):
            azimuth_time, _ = gso.earth2sat(target_coord)
            sensor_velocity[:, index_target] = gso.get_velocity(azimuth_time)
        sensor_velocity = np.reshape(
            sensor_velocity, [3, tuple(self.x_sar_coords.shape)[0], tuple(self.x_sar_coords.shape)[1]],
        )

        return sensor_velocity

    def get_sensor_velocity(self, state_vector_ext: StateVectors = None):
        if state_vector_ext is None:
            return self.get_sensor_velocity_master()
        else:
            return self.get_sensor_velocity_slave(state_vector_ext)

    def compute_off_nadir_angles(self, state_vector_ext: StateVectors = None):
        sensor_positions_master = self.get_sensor_positions_master()
        if state_vector_ext is None:
            sensor_positions = sensor_positions_master
        else:
            sensor_positions = self.get_sensor_positions_slave(state_vector_ext)

        return self.compute_off_nadir_angles_from_positions(sensor_positions_master, sensor_positions)

    def compute_off_nadir_angles_from_positions(self, sensor_positions_master, sensor_positions):
        # Compute nadir versors
        nadir_versors = sensor_positions_master[:, :, 0] / np.linalg.norm(sensor_positions_master[:, :, 0], axis=0)
        # Compute LOS versors
        los_versors_x = sensor_positions[0, :] - self.x_sar_coords
        los_versors_y = sensor_positions[1, :] - self.y_sar_coords
        los_versors_z = sensor_positions[2, :] - self.z_sar_coords
        distances = np.sqrt(np.power(los_versors_x, 2) + np.power(los_versors_y, 2) + np.power(los_versors_z, 2))
        los_versors_x /= distances
        los_versors_y /= distances
        los_versors_z /= distances
        # Compute off-nadir angle
        scalar_product = (
            los_versors_x * nadir_versors[0, :][:, np.newaxis]
            + los_versors_y * nadir_versors[1, :][:, np.newaxis]
            + los_versors_z * nadir_versors[2, :][:, np.newaxis]
        )
        return np.arccos(scalar_product)

    def compute_terrain_slope(self):
        sensor_positions = self.get_sensor_positions_master()[:, :, 0]
        satellite_velocity = self.get_sensor_velocity_master()[:, :, 0]
        satellite_velocity_versors = satellite_velocity / np.linalg.norm(satellite_velocity, axis=0)
        # Compute ground versors
        nadir_versors = -sensor_positions / np.linalg.norm(sensor_positions, axis=0)
        ground_versors = np.zeros_like(nadir_versors)
        ground_versors[0, :] = (
            nadir_versors[1, :] * satellite_velocity_versors[2, :]
            - nadir_versors[2, :] * satellite_velocity_versors[1, :]
        )
        ground_versors[1, :] = (
            nadir_versors[2, :] * satellite_velocity_versors[0, :]
            - nadir_versors[0, :] * satellite_velocity_versors[2, :]
        )
        ground_versors[2, :] = (
            nadir_versors[0, :] * satellite_velocity_versors[1, :]
            - nadir_versors[1, :] * satellite_velocity_versors[0, :]
        )

        # Compute terrain versors
        diff_x = np.diff(self.x_sar_coords, axis=1)
        diff_y = np.diff(self.y_sar_coords, axis=1)
        diff_z = np.diff(self.z_sar_coords, axis=1)
        diff_norm = np.sqrt(np.power(diff_x, 2) + np.power(diff_y, 2) + np.power(diff_z, 2))
        diff_x /= diff_norm
        diff_y /= diff_norm
        diff_z /= diff_norm

        terrain_vector_x = (
            diff_y * satellite_velocity_versors[2, :][:, np.newaxis]
            - diff_z * satellite_velocity_versors[1, :][:, np.newaxis]
        )
        terrain_vector_y = (
            diff_z * satellite_velocity_versors[0, :][:, np.newaxis]
            - diff_x * satellite_velocity_versors[2, :][:, np.newaxis]
        )
        terrain_vector_z = (
            diff_x * satellite_velocity_versors[1, :][:, np.newaxis]
            - diff_y * satellite_velocity_versors[0, :][:, np.newaxis]
        )

        # Compute terrain slope
        scalar_product = (
            terrain_vector_x * ground_versors[0, :][:, np.newaxis]
            + terrain_vector_y * ground_versors[1, :][:, np.newaxis]
            + terrain_vector_z * ground_versors[2, :][:, np.newaxis]
        )
        self._terrain_slope = -np.arcsin(scalar_product)
        self._terrain_slope = np.append(self._terrain_slope, self._terrain_slope[:, -1][:, np.newaxis], axis=1)
        # Interpolate to the data grid
        self._terrain_slope = shift_image(self._terrain_slope, 0, 1 / 2)

    def compute_xyz(self):
        # Get target coordinates
        targets_coord = self.get_targets_coordinates()

        # Display DEM
        if self.flag_display:
            plt.figure()
            plt.imshow(
                self.dem,
                extent=[self.lon_axis[0], self.lon_axis[-1], self.lat_axis[0], self.lat_axis[-1]],
                origin="lower",
            )

            plt.show()

        # Initialization
        azimuth_axis_rel, range_axis_out = get_axis_from_roi(self._raster_info, self.roi)

        range_axis_out = range_axis_out * constants.c / 2

        self._x_sar_coords = np.zeros([self.roi[2], self.roi[3]])
        self._y_sar_coords = np.zeros([self.roi[2], self.roi[3]])
        self._z_sar_coords = np.zeros([self.roi[2], self.roi[3]])

        self._dem_sar = np.zeros([self.roi[2], self.roi[3]])

        range_line_coords_out_grid = np.zeros([3, self.roi[3]])

        sat_positions = self.get_sensor_positions_master()
        sat_velocities = self.get_sensor_velocity_master()

        # Get first and last coords
        mean_velocity = np.nanmean(np.reshape(sat_velocities, [3, self.roi[2] * self.roi[3]]), axis=1)
        if is_north_heading(mean_velocity):
            targets_coord_first = self.get_targets_coordinates(lon_index=0)
            targets_coord_last = self.get_targets_coordinates(lon_index=-1)
        else:
            targets_coord_first = self.get_targets_coordinates(lat_index=0)
            targets_coord_last = self.get_targets_coordinates(lat_index=-1)

        # Slow time loop
        for index_slow_time in progressbar(range(azimuth_axis_rel.size), max_value=azimuth_axis_rel.size):

            current_sat_position = sat_positions[:, index_slow_time, 0][:, np.newaxis]
            current_sat_velocity = sat_velocities[:, index_slow_time, 0][:, np.newaxis]

            zero_doppler_plane = get_zero_doppler_plane(current_sat_velocity, current_sat_position)

            lat_indexes, lon_indexes, current_indexes = self.get_indexes_dem_portion(
                zero_doppler_plane, targets_coord_first, targets_coord_last, is_north_heading(mean_velocity),
            )
            lat_axis_ = self.lat_axis[lat_indexes]
            lon_axis_ = self.lon_axis[lon_indexes]

            # Distance from zero doppler plane
            distances = get_distance_from_zero_doppler(zero_doppler_plane, targets_coord[:, current_indexes])
            distances = np.reshape(distances, [lat_axis_.size, lon_axis_.size])

            # Display distances
            if self.flag_display:
                plt.figure()
                plt.imshow(np.abs(distances), aspect="auto")
                plt.colorbar()
                plt.show()
                import sys

                sys.exit()

            # Interpolation along latitude
            if is_north_heading(mean_velocity):
                index_min_dist_vect = np.zeros([2, self.lon_axis.size], dtype=int)
                range_line_coords = np.zeros([3, self.lon_axis.size])
                for current_lon_index in range(self.lon_axis.size):
                    current_distances = distances[:, current_lon_index]
                    if np.any(np.sign(current_distances) < 0) and np.any(np.sign(current_distances) > 0):
                        index_min_dist_vect[:, current_lon_index] = get_min_indexes(current_distances)
                    else:
                        index_min_dist_vect[:, current_lon_index] = [0, 1]

                distances_vect = np.concatenate(
                    [
                        distances.flatten()[
                            np.ravel_multi_index(
                                [index_min_dist_vect[0], np.arange(self.lon_axis.size, dtype=int)], distances.shape,
                            )
                        ][np.newaxis, :],
                        distances.flatten()[
                            np.ravel_multi_index(
                                [index_min_dist_vect[1], np.arange(self.lon_axis.size, dtype=int)], distances.shape,
                            )
                        ][np.newaxis, :],
                    ]
                )
                dem_vect = np.concatenate(
                    [
                        self.dem.flatten()[
                            np.ravel_multi_index(
                                [lat_indexes[index_min_dist_vect[0]], np.arange(self.lon_axis.size, dtype=int),],
                                self.dem.shape,
                            )
                        ][np.newaxis, :],
                        self.dem.flatten()[
                            np.ravel_multi_index(
                                [lat_indexes[index_min_dist_vect[1]], np.arange(self.lon_axis.size, dtype=int),],
                                self.dem.shape,
                            )
                        ][np.newaxis, :],
                    ]
                )

                range_line_coords[0, :] = lerp(distances_vect, lat_axis_[index_min_dist_vect])
                range_line_coords[1, :] = self.lon_axis
                range_line_coords[2, :] = lerp(distances_vect, dem_vect)
            else:
                index_min_dist_vect = np.zeros([2, self.lat_axis.size], dtype=int)
                range_line_coords = np.zeros([3, self.lat_axis.size])
                for current_lat_index in range(self.lat_axis.size):
                    current_distances = distances[current_lat_index, :]
                    if np.any(np.sign(current_distances) < 0) and np.any(np.sign(current_distances) > 0):
                        index_min_dist_vect[:, current_lat_index] = get_min_indexes(current_distances)
                    else:
                        index_min_dist_vect[:, current_lat_index] = [0, 1]
                complete_lat_indexes = np.arange(self.lat_axis.size, dtype=int)
                distances_vect = np.concatenate(
                    [
                        distances.flatten()[
                            np.ravel_multi_index([complete_lat_indexes, index_min_dist_vect[0]], distances.shape)
                        ][np.newaxis, :],
                        distances.flatten()[
                            np.ravel_multi_index([complete_lat_indexes, index_min_dist_vect[1]], distances.shape)
                        ][np.newaxis, :],
                    ],
                )
                dem_vect = np.concatenate(
                    [
                        self.dem.flatten()[
                            np.ravel_multi_index(
                                [complete_lat_indexes, lon_indexes[index_min_dist_vect[0]]], self.dem.shape,
                            )
                        ][np.newaxis, :],
                        self.dem.flatten()[
                            np.ravel_multi_index(
                                [complete_lat_indexes, lon_indexes[index_min_dist_vect[1]]], self.dem.shape,
                            )
                        ][np.newaxis, :],
                    ]
                )

                range_line_coords[0, :] = self.lat_axis
                range_line_coords[1, :] = lerp(distances_vect, lon_axis_[index_min_dist_vect])
                range_line_coords[2, :] = lerp(distances_vect, dem_vect)

            range_line_coords_xyz = conversions.llh2xyz(range_line_coords)

            target_dist = np.sqrt(np.sum(np.power((range_line_coords_xyz - current_sat_position), 2), 0))

            # Interpolate along range
            range_line_coords = range_line_coords[:, ~np.isnan(target_dist)]
            target_dist = target_dist[~np.isnan(target_dist)]
            if len(target_dist):
                range_line_coords_out_grid[0, :] = np.interp(range_axis_out, target_dist, range_line_coords[0, :])
                range_line_coords_out_grid[1, :] = np.interp(range_axis_out, target_dist, range_line_coords[1, :])
                range_line_coords_out_grid[2, :] = np.interp(range_axis_out, target_dist, range_line_coords[2, :])

                self._dem_sar[index_slow_time, :] = range_line_coords_out_grid[2, :]

                range_line_coords_out_grid = conversions.llh2xyz(range_line_coords_out_grid)

                self._x_sar_coords[index_slow_time, :] = range_line_coords_out_grid[0, :]
                self._y_sar_coords[index_slow_time, :] = range_line_coords_out_grid[1, :]
                self._z_sar_coords[index_slow_time, :] = range_line_coords_out_grid[2, :]


def get_slave_state_vectors(path_slave_pf, channel_index):
    pf_slave = ProductFolder(path_slave_pf)
    ch_slave = pf_slave.get_channel(channel_index)
    index_last_md_channel = ch_slave.metadata.get_number_of_channels() - 1
    ri_slave = ch_slave.get_raster_info(index_last_md_channel)
    sv_slave = cut_state_vectors(
        ch_slave.get_state_vectors(index_last_md_channel),
        ri_slave.lines_start,
        ri_slave.lines_start + ri_slave.lines * ri_slave.lines_step,
    )
    return sv_slave


def main(input_file_path):
    # from src.inout.XmlIO import XmlIO
    # from src.miscellanea.ext_geodata_mosaic import ext_geodata_mosaic
    input_obj = XmlIO(input_file_path)
    configuration_obj = XmlIO(input_obj.configuration_file)

    channel_index = 0  # TODO it is assumed that the product contains 1 data channel
    pf_master = ProductFolder(input_obj.master_path)
    ch_master = pf_master.get_channel(channel_index)
    raster_info = ch_master.get_raster_info()
    state_vectors = ch_master.get_state_vectors()
    side_looking = ch_master.get_dataset_info().side_looking.value

    if "roi_master" in input_obj:
        roi_master = [int(n) for n in np.array(input_obj.roi_master.split(","), dtype="int")]
        raster_info = cut_raster_info(raster_info, roi_master)

    if "geoid_repository" in configuration_obj:
        geoid_dir = configuration_obj.geoid_repository
    else:
        geoid_dir = None

    if not "dem_path" in input_obj:
        # Get DEM ROI

        geo_corners = get_geo_corners(raster_info, state_vectors, side_looking)
        geo_corners = np.rad2deg(geo_corners)

        # Extract DEM
        external_dem_obj = ext_geodata_mosaic(
            *geo_corners,
            configuration_obj.dem_repository,
            configuration_obj.dem_type,
            input_obj.output_path,
            geoid_dir,
        )
        external_dem_obj.Run()
        dem_path = external_dem_obj.geodata_output_dir
    else:
        dem_path = input_obj.dem_path

    pf_dem = ProductFolder(dem_path)
    ch_dem = pf_dem.get_channel(0)
    ri_dem = ch_dem.get_raster_info(0)

    dem_samples = ri_dem.samples
    lat_axis = (np.arange(0, ri_dem.lines) * ri_dem.lines_step + ri_dem.lines_start) * np.pi / 180
    lon_axis = (np.arange(0, dem_samples) * ri_dem.samples_step + ri_dem.samples_start) * np.pi / 180

    dem_data = pf_dem.read_data(0, [0, 0, ri_dem.lines, dem_samples])

    # Run SARGeometry
    sar_geometry = SARGeometry(raster_info, state_vectors, lat_axis, lon_axis, dem_data)

    # XYZ computation
    if (
        int(configuration_obj.save_xyz)
        or int(configuration_obj.save_reference_height)
        or int(configuration_obj.save_distance)
        or int(configuration_obj.save_off_nadir)
        or int(configuration_obj.save_terrain_slope)
        or int(configuration_obj.save_wavenumber)
    ):
        sar_geometry.compute_xyz()
        azimuth_axis, range_axis = sar_geometry.get_axis()

    # Various products computation
    if (
        int(configuration_obj.save_distance)
        or int(configuration_obj.save_off_nadir)
        or int(configuration_obj.save_wavenumber)
    ):
        sv_slave = get_slave_state_vectors(input_obj.slave_path, channel_index)
        gso_slave = create_general_sar_orbit(sv_slave)

        targets_coords = np.vstack(
            (
                sar_geometry.x_sar_coords.flatten(),
                sar_geometry.y_sar_coords.flatten(),
                sar_geometry.z_sar_coords.flatten(),
            )
        )

        azimuth_points_master_slave, range_points_master_slave = perform_inverse_geocoding(gso_slave, targets_coords)

    if int(configuration_obj.save_distance):
        distance_master = sar_geometry.compute_distance_master()
        distance_slave = range_points_master_slave * constants.c / 2
        distance_slave.shape = sar_geometry.x_sar_coords.shape
    else:
        distance_master = None
        distance_slave = None

    if int(configuration_obj.save_off_nadir) or int(configuration_obj.save_wavenumber):
        sensor_positions_master = sar_geometry.get_sensor_positions_master()
        off_nadir_angles_master = sar_geometry.compute_off_nadir_angles_from_positions(
            sensor_positions_master, sensor_positions_master
        )

        sensor_positions_slave = gso_slave.get_position(azimuth_points_master_slave)
        sensor_positions_slave.shape = (3, *sar_geometry.x_sar_coords.shape)
        off_nadir_angles_slave = sar_geometry.compute_off_nadir_angles_from_positions(
            sensor_positions_master, sensor_positions_slave
        )
    else:
        off_nadir_angles_master = None
        off_nadir_angles_slave = None

    if int(configuration_obj.save_terrain_slope) or int(configuration_obj.save_wavenumber):
        sar_geometry.compute_terrain_slope()

    if int(configuration_obj.save_wavenumber):
        wavelength = constants.c / ch_master.get_dataset_info(0).fc_hz
        incidence_angles = off_nadir_angles_master - sar_geometry.terrain_slope
        vertical_wave_number = compute_vertical_wavenumber_angles(
            wavelength, off_nadir_angles_master, off_nadir_angles_slave,
        )

    # Save output products
    if int(configuration_obj.save_xyz):
        save_product(
            [sar_geometry.x_sar_coords, sar_geometry.y_sar_coords, sar_geometry.z_sar_coords],
            azimuth_axis,
            range_axis,
            input_obj.output_path,
            EProductType("ECEFGRID"),
        )

    if int(configuration_obj.save_reference_height):
        save_product(
            [sar_geometry.dem_sar], azimuth_axis, range_axis, input_obj.output_path, EProductType("REFERENCEHEIGHT"),
        )

    if int(configuration_obj.save_distance):
        save_product(
            [distance_master, distance_slave],
            azimuth_axis,
            range_axis,
            input_obj.output_path,
            EProductType("DISTANCE"),
        )

    if int(configuration_obj.save_off_nadir):
        save_product(
            [off_nadir_angles_master, off_nadir_angles_slave],
            azimuth_axis,
            range_axis,
            input_obj.output_path,
            EProductType("OFFNADIR"),
        )

    if int(configuration_obj.save_terrain_slope):
        save_product(
            [sar_geometry.terrain_slope], azimuth_axis, range_axis, input_obj.output_path, EProductType("SLOPE"),
        )

    if int(configuration_obj.save_wavenumber):
        save_product(
            [vertical_wave_number], azimuth_axis, range_axis, input_obj.output_path, EProductType("WAVENUMBER"),
        )


if __name__ == "__main__":
    import sys

    main(sys.argv[1])
