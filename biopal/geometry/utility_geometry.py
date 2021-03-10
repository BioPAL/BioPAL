import os
import logging
import shutil
import numpy as np
from arepytools.io.productfolder import ProductFolder
from arepytools.io.metadata import DataSetInfo
from arepytools.timing.precisedatetime import PreciseDateTime
from arepytools.geometry.generalsarorbit import create_general_sar_orbit

from biopal.utility.utility_functions import (
    convert_rasterinfo_meters_to_seconds,
    convert_rasterinfo_seconds_to_meters,
)
from biopal.geometry import geometric_lib
from biopal.data_operations.data_operations import data_oversample
from biopal.utility.constants import (
    OVERSAMPLING_FACTOR,
    LIGHTSPEED,
)
from biopal.io.data_io import readBiomassHeader


def compute_and_oversample_geometry_auxiliaries(
    L1c_repository,
    aux_file_names,
    stack_id,
    acquisition_pf_names,
    master_id,
    enable_resampling,
    comp_ref_h=True,
    comp_dist=True,
    force_ellipsoid=False,
    sar_geometry_master=None,
):
    ecef_grid = None
    off_nadir_angle_rad = None
    kz = None
    reference_height = None
    R = None
    slope = None

    logging.info("Geometry library: loading auxiliary data: DEM...")
    pf_name = os.path.basename(aux_file_names.dem_folder)
    folder = os.path.dirname(aux_file_names.dem_folder)
    path_dem = os.path.join(folder, pf_name)

    try:
        pf_dem = ProductFolder(path_dem, "r")
    except Exception as e:
        logging.error("Geometry library:  error during auxiliary DEM reading: " + str(e), exc_info=True)
        raise

    ch_dem = pf_dem.get_channel(0)
    ri_dem = ch_dem.get_raster_info(0)

    dem_samples = ri_dem.samples
    lat_axis = (np.arange(0, ri_dem.lines) * ri_dem.lines_step + ri_dem.lines_start) * np.pi / 180
    lon_axis = (np.arange(0, dem_samples) * ri_dem.samples_step + ri_dem.samples_start) * np.pi / 180

    dem = pf_dem.read_data(0, [0, 0, ri_dem.lines, dem_samples])

    if force_ellipsoid:
        logging.info(
            "Geometry library: the ellipsoid altitude has been forced by input flag, only the slope over ellipsoid will be computed, after used to correcting geometry global / local reference."
        )
        dem = 0 * dem

    logging.info("...done\n")

    logging.info("Geometry library: loading orbits from stack metadata...")
    ch_list = []
    dsi_list = []
    swath_id_list = []
    swath_info_list = []

    # read master alone to absure it will be first in the list
    pf_master = ProductFolder(os.path.join(L1c_repository, master_id), "r")
    ch_master = pf_master.get_channel(0)
    ch_list.append(ch_master)
    ri_temp = ch_master.get_raster_info(0)
    
    (
    _, 
    _, 
    _, 
    _, 
    _, 
    _, 
    _, 
    _, 
    sensor_velocity,
    )= readBiomassHeader( ProductFolder( os.path.join(L1c_repository, master_id), 'r'), 0)

    if ri_temp.samples_step_unit == "m":
        ri_master = convert_rasterinfo_meters_to_seconds(ri_temp, sensor_velocity)

        # ri_list.append( ri_master )
        num_samples = ri_temp.samples
        num_lines = ri_temp.lines
        lines_start = ri_temp.lines_start
        lines_start_unit = ri_temp.lines_start_unit
        lines_step = ri_temp.lines_step
        lines_step_unit = ri_temp.lines_step_unit
        samples_start = ri_temp.samples_start
        samples_start_unit = ri_temp.samples_start_unit
        samples_step = ri_temp.samples_step
        samples_step_unit = ri_temp.samples_step_unit

    else:
        ri_master = ri_temp

        ri_meters = convert_rasterinfo_seconds_to_meters(ri_temp, sensor_velocity)
        # ri_list.append( ri_master )
        num_samples = ri_meters.samples
        num_lines = ri_meters.lines
        lines_start = ri_meters.lines_start
        lines_start_unit = ri_meters.lines_start_unit
        lines_step = ri_meters.lines_step
        lines_step_unit = ri_meters.lines_step_unit
        samples_start = ri_meters.samples_start
        samples_start_unit = ri_meters.samples_start_unit
        samples_step = ri_meters.samples_step
        samples_step_unit = ri_meters.samples_step_unit

    sv_master = geometric_lib.cut_state_vectors(
        ch_master.get_state_vectors(0),
        ri_master.lines_start,
        ri_master.lines_start + ri_master.lines * ri_master.lines_step,
    )

    # Initiaize sar geometry
    flag_compute_sar_geometry = False
    if not sar_geometry_master:
        flag_compute_sar_geometry = True
        sar_geometry_master = geometric_lib.SARGeometry(ri_master, sv_master, lat_axis, lon_axis, dem)
        sar_geometry_master.compute_xyz()

    dsi_master = ch_master.get_dataset_info(0)
    dsi_list.append(dsi_master)

    swath_id_list.append(ch_master.get_swath_info(0).swath)
    swath_info_list.append(ch_master.get_swath_info(0))

    for acq_id in acquisition_pf_names:
        if not acq_id == master_id:
            pf_slave = ProductFolder(os.path.join(L1c_repository, acq_id), "r")
            ch_slave = pf_slave.get_channel(0)
            ch_list.append(ch_slave)

            ri_temp = ch_slave.get_raster_info(0)
            if ri_temp.samples_step_unit == "m":
                ri_slave = convert_rasterinfo_meters_to_seconds(ri_temp, sensor_velocity)
            else:
                ri_slave = ri_temp

            # ri_list. append( ri_slave )
            dsi_list.append(ch_slave.get_dataset_info(0))
            swath_id_list.append(ch_slave.get_swath_info(0).swath)
            swath_info_list.append(ch_slave.get_swath_info(0))

    logging.info("...done\n")

    logging.info("Geometry library: computing auxiliary data: slope...")

    if flag_compute_sar_geometry:
        sar_geometry_master.compute_terrain_slope()

    slope = sar_geometry_master.terrain_slope

    if enable_resampling:
        slope = data_oversample(slope, OVERSAMPLING_FACTOR, ri_master)[0]

    if force_ellipsoid:

        slope = slope.transpose()

        return ecef_grid, off_nadir_angle_rad, slope, kz, reference_height, R, sar_geometry_master

    else:
        pf_name = os.path.basename(aux_file_names.slope_file_names[stack_id])
        folder_name = os.path.dirname(aux_file_names.slope_file_names[stack_id])
        out_pf_name = os.path.join(folder_name, stack_id)

        if os.path.exists(out_pf_name):
            logging.warning("Overwriting of Auxiliary slope for stack " + stack_id)
            shutil.rmtree(out_pf_name)
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        # slope = sar_geometry_master.terrain_slope
        pf = ProductFolder(out_pf_name, "w")

        # append a channel with the correct dimensions
        pf.append_channel(num_lines, num_samples, "FLOAT32", header_offset=0)

        # write the ECEF coordinate to the channel
        pf.write_data(0, slope)

        # prepare the metadata elements
        data_channel_obj = pf.get_channel(0)
        metadata_obj = data_channel_obj.metadata
        metadatachannel_obj = metadata_obj.get_metadata_channels(0)

        # Raster Info
        ri = metadatachannel_obj.get_element("RasterInfo")
        ri.set_lines_axis(lines_start, lines_start_unit, lines_step, lines_step_unit)
        ri.set_samples_axis(samples_start, samples_start_unit, samples_step, samples_step_unit)

        # DataSetInfo: to insert the X, Y or Z description
        dsi_list[0].description = "Auxiliary data: Slope [rad]"
        metadatachannel_obj.insert_element(dsi_list[0])

        pf.write_metadata(0)

        pf = None

    print("...done.\n")

    logging.info("Geometry library: computing auxiliary data: ECEFGRID...")

    ecef_grid = {
        "X": sar_geometry_master.x_sar_coords,
        "Y": sar_geometry_master.y_sar_coords,
        "Z": sar_geometry_master.z_sar_coords,
    }

    if enable_resampling:
        ecef_grid = data_oversample(ecef_grid, OVERSAMPLING_FACTOR, ri_master)[0]

    pf_name = os.path.basename(aux_file_names.ECEF_grid_file_names[stack_id])
    folder_name = os.path.dirname(aux_file_names.ECEF_grid_file_names[stack_id])
    out_pf_name = os.path.join(folder_name, stack_id)

    if os.path.exists(out_pf_name):
        logging.warning("Overwriting of Auxiliary ECEFGRID for stack " + stack_id)
        shutil.rmtree(out_pf_name)
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    pf = ProductFolder(out_pf_name, "w")
    dsi_ecefgrid = DataSetInfo("STRIPMAP", 435000000)
    dsi_ecefgrid.description = ""
    dsi_ecefgrid.image_type = "SLC"
    dsi_ecefgrid.projection = "SLANT RANGE"
    dsi_ecefgrid.side_looking = "RIGHT"
    dsi_ecefgrid.sensor_name = "BIOMASS"
    dsi_ecefgrid.sense_date = PreciseDateTime()
    dsi_ecefgrid.acquisition_station = "NOT_AVAILABLE"
    dsi_ecefgrid.processing_center = "NOT_AVAILABLE"
    dsi_ecefgrid.processing_date = PreciseDateTime()
    dsi_ecefgrid.processing_software = "NOT_AVAILABLE"

    coord_names = ["X", "Y", "Z"]
    for channel_idx, coord_name in enumerate(coord_names):

        pf.append_channel(num_lines, num_samples, "FLOAT32", header_offset=0)
        pf.write_data(channel_idx, ecef_grid[coord_name])

        # prepare the metadata elements
        data_channel_obj = pf.get_channel(channel_idx)
        metadata_obj = data_channel_obj.metadata
        metadatachannel_obj = metadata_obj.get_metadata_channels(0)

        # Raster Info
        ri = metadatachannel_obj.get_element("RasterInfo")
        ri.set_lines_axis(lines_start, lines_start_unit, lines_step, lines_step_unit)
        ri.set_samples_axis(samples_start, samples_start_unit, samples_step, samples_step_unit)

        # DataSetInfo: to insert the X, Y or Z description
        dsi_ecefgrid.description = "Auxiliary data: " + coord_name + " ECEF GRID [m]"
        metadatachannel_obj.insert_element(dsi_ecefgrid)

        pf.write_metadata(channel_idx)

    pf = None
    logging.info("...done\n")

    targets_coords = np.vstack(
        (
            sar_geometry_master.x_sar_coords.flatten(),
            sar_geometry_master.y_sar_coords.flatten(),
            sar_geometry_master.z_sar_coords.flatten(),
        )
    )

    if comp_ref_h:
        logging.info("Geometry library: computing auxiliary data: reference height...")
        reference_height = sar_geometry_master.dem_sar

        if enable_resampling:
            reference_height = data_oversample(reference_height, OVERSAMPLING_FACTOR, ri_master)[0]

        pf_name = os.path.basename(aux_file_names.reference_height_file_names[stack_id])
        folder_name = os.path.dirname(aux_file_names.reference_height_file_names[stack_id])
        out_pf_name = os.path.join(folder_name, stack_id)

        if os.path.exists(out_pf_name):
            logging.warning("Overwriting of Auxiliary reference height for stack " + stack_id)
            shutil.rmtree(out_pf_name)
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        # create the Product Folder in writing mode:
        pf = ProductFolder(out_pf_name, "w")
        # append a channel with the correct dimensions
        pf.append_channel(num_lines, num_samples, "FLOAT32", header_offset=0)
        # write the DTM to PF (reference height)
        pf.write_data(0, reference_height)

        # prepare the metadata elements
        data_channel_obj = pf.get_channel(0)
        metadata_obj = data_channel_obj.metadata
        metadatachannel_obj = metadata_obj.get_metadata_channels(0)

        # Raster Info
        ri = metadatachannel_obj.get_element("RasterInfo")
        ri.set_lines_axis(lines_start, lines_start_unit, lines_step, lines_step_unit)

        ri.set_samples_axis(samples_start, samples_start_unit, samples_step, samples_step_unit)

        # DataSetInfo
        dsi = DataSetInfo("STRIPMAP", 435000000)
        dsi.description = "Auxiliary data: Reference height (DTM)"
        dsi.image_type = "SLC"
        dsi.projection = "SLANT RANGE"
        dsi.side_looking = "RIGHT"
        dsi.sensor_name = "BIOMASS"
        dsi.sense_date = PreciseDateTime()
        dsi.acquisition_station = "NOT_AVAILABLE"
        dsi.processing_center = "NOT_AVAILABLE"
        dsi.processing_date = PreciseDateTime()
        dsi.processing_software = "NOT_AVAILABLE"

        metadatachannel_obj.insert_element(dsi)

        pf.write_metadata(0)
        pf = None

        print("...done.\n")

    if comp_dist:
        logging.info("Geometry library: computing auxiliary data: slant range distances...")

        pf_name = os.path.basename(aux_file_names.slant_range_distances_file_names[stack_id])
        folder_name = os.path.dirname(aux_file_names.slant_range_distances_file_names[stack_id])
        out_pf_name = os.path.join(folder_name, stack_id)

        if os.path.exists(out_pf_name):
            logging.warning("Overwriting of Auxiliary slant range distances for stack " + stack_id)
            shutil.rmtree(out_pf_name)
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        pf = ProductFolder(out_pf_name, "w")
        R = {}
        for channel_idx, channel_curr in enumerate(ch_list):

            if channel_idx == 0:
                # distance between targets and master orbit
                R[swath_id_list[channel_idx]] = sar_geometry_master.compute_distance_master()
            else:

                # distance between targets and slave orbit
                ri_temp = channel_curr.get_raster_info(0)
                if ri_temp.samples_step_unit == "m":
                    ri_slave = convert_rasterinfo_meters_to_seconds(ri_temp, sensor_velocity)
                else:
                    ri_slave = ri_temp

                sv_slave = geometric_lib.cut_state_vectors(
                    channel_curr.get_state_vectors(0),
                    ri_slave.lines_start,
                    ri_slave.lines_start + ri_slave.lines * ri_slave.lines_step,
                )

                gso_slave = create_general_sar_orbit(sv_slave)

                _, range_points_master_slave = geometric_lib.perform_inverse_geocoding(gso_slave, targets_coords)

                distance_slave = range_points_master_slave * LIGHTSPEED / 2
                distance_slave.shape = sar_geometry_master.x_sar_coords.shape

                R[swath_id_list[channel_idx]] = distance_slave

            if enable_resampling:
                R[swath_id_list[channel_idx]] = data_oversample(
                    R[swath_id_list[channel_idx]], OVERSAMPLING_FACTOR, ri_master
                )[0]

            # append a channel with the correct dimensions
            pf.append_channel(num_lines, num_samples, "FLOAT32", header_offset=0)

            # write the R coordinate to the channel
            pf.write_data(channel_idx, R[swath_id_list[channel_idx]])

            # prepare the metadata elements
            data_channel_obj = pf.get_channel(channel_idx)
            metadata_obj = data_channel_obj.metadata
            metadatachannel_obj = metadata_obj.get_metadata_channels(0)

            # Raster Info
            ri = metadatachannel_obj.get_element("RasterInfo")
            ri.set_lines_axis(lines_start, lines_start_unit, lines_step, lines_step_unit)
            ri.set_samples_axis(samples_start, samples_start_unit, samples_step, samples_step_unit)

            # DataSetInfo: to insert the X, Y or Z description
            dsi_list[channel_idx].description = "Auxiliary data: slant range distances [m]"
            metadatachannel_obj.insert_element(dsi_list[channel_idx])

            metadatachannel_obj.insert_element(swath_info_list[channel_idx])

            pf.write_metadata(channel_idx)
        pf = None

        print("...done.\n")

    logging.info("Geometry library: computing auxiliary data: off nadir...")

    pf_name = os.path.basename(aux_file_names.off_nadir_angle_file_names[stack_id])
    folder_name = os.path.dirname(aux_file_names.off_nadir_angle_file_names[stack_id])
    out_pf_name = os.path.join(folder_name, stack_id)

    if os.path.exists(out_pf_name):
        logging.warning("Overwriting of Auxiliary off nadir for stack " + stack_id)
        shutil.rmtree(out_pf_name)
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    pf = ProductFolder(out_pf_name, "w")
    off_nadir_angle_rad = {}
    for channel_idx, channel_curr in enumerate(ch_list):

        if channel_idx == 0:
            # master
            sensor_positions_master = sar_geometry_master.get_sensor_positions_master()
            off_nadir_angle_rad[
                swath_id_list[channel_idx]
            ] = sar_geometry_master.compute_off_nadir_angles_from_positions(
                sensor_positions_master, sensor_positions_master
            )
        else:
            # slaves
            ri_temp = channel_curr.get_raster_info(0)
            if ri_temp.samples_step_unit == "m":
                ri_slave = convert_rasterinfo_meters_to_seconds(ri_temp, sensor_velocity)
            else:
                ri_slave = ri_temp

            sv_slave = geometric_lib.cut_state_vectors(
                channel_curr.get_state_vectors(0),
                ri_slave.lines_start,
                ri_slave.lines_start + ri_slave.lines * ri_slave.lines_step,
            )

            gso_slave = create_general_sar_orbit(sv_slave)

            azimuth_points_master_slave, _ = geometric_lib.perform_inverse_geocoding(gso_slave, targets_coords)

            sensor_positions_slave = gso_slave.get_position(azimuth_points_master_slave)
            sensor_positions_slave.shape = (3, *sar_geometry_master.x_sar_coords.shape)
            off_nadir_angle_rad[
                swath_id_list[channel_idx]
            ] = sar_geometry_master.compute_off_nadir_angles_from_positions(
                sensor_positions_master, sensor_positions_slave
            )

        if enable_resampling:
            off_nadir_angle_rad[swath_id_list[channel_idx]] = data_oversample(
                off_nadir_angle_rad[swath_id_list[channel_idx]], OVERSAMPLING_FACTOR, ri_master
            )[0]

        # append a channel with the correct dimensions
        pf.append_channel(num_lines, num_samples, "FLOAT32", header_offset=0)

        # write the ECEF coordinate to the channel
        pf.write_data(channel_idx, off_nadir_angle_rad[swath_id_list[channel_idx]])

        # prepare the metadata elements
        data_channel_obj = pf.get_channel(channel_idx)
        metadata_obj = data_channel_obj.metadata
        metadatachannel_obj = metadata_obj.get_metadata_channels(0)

        # Raster Info
        ri = metadatachannel_obj.get_element("RasterInfo")
        ri.set_lines_axis(lines_start, lines_start_unit, lines_step, lines_step_unit)
        ri.set_samples_axis(samples_start, samples_start_unit, samples_step, samples_step_unit)

        # DataSetInfo: to insert the X, Y or Z description
        dsi_list[channel_idx].description = "Auxiliary data: Incidence Angle [rad]"
        metadatachannel_obj.insert_element(dsi_list[channel_idx])

        metadatachannel_obj.insert_element(swath_info_list[channel_idx])

        pf.write_metadata(channel_idx)

    pf = None
    print("...done.\n)")

    logging.info("Geometry library: computing auxiliary data: KZ (wavenumbers)...")

    pf_name = os.path.basename(aux_file_names.kz_file_names[stack_id])
    folder_name = os.path.dirname(aux_file_names.kz_file_names[stack_id])
    out_pf_name = os.path.join(folder_name, stack_id)

    if os.path.exists(out_pf_name):
        logging.warning("Overwriting of Auxiliary KZ (wavenumbers) for stack " + stack_id)
        shutil.rmtree(out_pf_name)
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    pf = ProductFolder(out_pf_name, "w")
    kz = {}
    for channel_idx, channel_curr in enumerate(ch_list):

        if channel_idx == 0:
            # Master
            wavelength = LIGHTSPEED / ch_master.get_dataset_info(0).fc_hz

            kz[swath_id_list[channel_idx]] = np.zeros((ri_master.lines, ri_master.samples))
        else:
            kz[swath_id_list[channel_idx]] = geometric_lib.compute_vertical_wavenumber_angles(
                wavelength, off_nadir_angle_rad[swath_id_list[0]], off_nadir_angle_rad[swath_id_list[channel_idx]],
            )

        if enable_resampling:
            kz[swath_id_list[channel_idx]] = data_oversample(
                kz[swath_id_list[channel_idx]], OVERSAMPLING_FACTOR, ri_master
            )[0]

        # append a channel with the correct dimensions
        pf.append_channel(num_lines, num_samples, "FLOAT32", header_offset=0)

        # write the KZ  to the channel
        pf.write_data(channel_idx, kz[swath_id_list[channel_idx]])

        # prepare the metadata elements
        data_channel_obj = pf.get_channel(channel_idx)
        metadata_obj = data_channel_obj.metadata
        metadatachannel_obj = metadata_obj.get_metadata_channels(0)

        # Raster Info
        ri = metadatachannel_obj.get_element("RasterInfo")
        ri.set_lines_axis(lines_start, lines_start_unit, lines_step, lines_step_unit)
        ri.set_samples_axis(samples_start, samples_start_unit, samples_step, samples_step_unit)

        # DataSetInfo: to insert the X, Y or Z description
        dsi_list[channel_idx].description = "Auxiliary data: KZ (Interferometric wavenumber matrix)"
        metadatachannel_obj.insert_element(dsi_list[channel_idx])

        metadatachannel_obj.insert_element(swath_info_list[channel_idx])

        pf.write_metadata(channel_idx)

    pf = None
    logging.info("...done\n")

    # Transpose all the outputs as done in the case of auxiliaries directly read (and not computed), see "read_and_oversample_aux_data"
    if not ecef_grid is None:
        for coord_name in coord_names:
            ecef_grid[coord_name] = ecef_grid[coord_name].transpose()

    for swath_id in swath_id_list:
        if not off_nadir_angle_rad is None:
            off_nadir_angle_rad[swath_id] = off_nadir_angle_rad[swath_id].transpose()
        if not kz is None:
            kz[swath_id] = kz[swath_id].transpose()
        if not R is None:
            R[swath_id] = R[swath_id].transpose()

    if not slope is None:
        slope = slope.transpose()

    if not reference_height is None:
        reference_height = reference_height.transpose()

    return ecef_grid, off_nadir_angle_rad, slope, kz, reference_height, R, sar_geometry_master
