import numpy as np
from . import pyxb_metadata_types, pyxb_metadata
from .. import metadata
from ...timing.precisedatetime import PreciseDateTime
from arepytools import constants
import logging
import pyxb


def convert_sampling_constants_to_aresys(i_sampling_constants: pyxb_metadata_types.SamplingConstantsType) -> metadata.SamplingConstants:
    o_sampling_constants = metadata.SamplingConstants(i_sampling_constants.frg_hz.value(),
                                                      i_sampling_constants.Brg_hz.value(),
                                                      i_sampling_constants.faz_hz.value(),
                                                      i_sampling_constants.Baz_hz.value())

    return o_sampling_constants


def convert_pulse_to_aresys(i_pulse: pyxb_metadata_types.PulseType) -> metadata.Pulse:
    o_pulse = metadata.Pulse(i_pulse.PulseLength.value(),
                             i_pulse.Bandwidth.value(),
                             i_pulse.PulseSamplingRate.value(),
                             i_pulse.PulseEnergy.value(),
                             )
    o_pulse.pulse_direction = i_pulse.Direction
    if i_pulse.PulseStartFrequency is not None:
        o_pulse.pulse_start_frequency = i_pulse.PulseStartFrequency.value()
    if i_pulse.PulseStartPhase is not None:
        o_pulse.pulse_start_phase = i_pulse.PulseStartPhase.value()

    return o_pulse


def convert_swath_info_to_aresys(i_swath_info: pyxb_metadata_types.SwathInfoType) -> metadata.SwathInfo:
    o_swath_info = metadata.SwathInfo(i_swath_info.Swath.value(),
                                      i_swath_info.Polarization,
                                      i_swath_info.AcquisitionPRF)

    o_swath_info.swath_acquisition_order = i_swath_info.SwathAcquisitionOrder.value()
    o_swath_info.rank = i_swath_info.Rank.value()
    o_swath_info.range_delay_bias = i_swath_info.RangeDelayBias.value()
    o_swath_info.acquisition_start_time = PreciseDateTime().set_from_utc_string(i_swath_info.AcquisitionStartTime.value())
    o_swath_info.azimuth_steering_rate_reference_time = i_swath_info.AzimuthSteeringRateReferenceTime.value()
    o_swath_info.azimuth_steering_rate_pol = (i_swath_info.AzimuthSteeringRatePol.val[0].value(),
                                              i_swath_info.AzimuthSteeringRatePol.val[1].value(),
                                              i_swath_info.AzimuthSteeringRatePol.val[2].value())
    o_swath_info.echoes_per_burst = i_swath_info.EchoesPerBurst
    if i_swath_info.ChannelDelay is not None:
        o_swath_info.channel_delay = i_swath_info.ChannelDelay.value()
    if i_swath_info.RxGain is not None:
        o_swath_info.rx_gain = i_swath_info.RxGain
    return o_swath_info


def convert_rasterinfo_to_aresys(i_rasterinfo: pyxb_metadata_types.RasterInfoType):
    samples = i_rasterinfo.Samples
    samples_start = float(i_rasterinfo.SamplesStart.value())
    samples_start_unit = i_rasterinfo.SamplesStart.unit
    samples_step = float(i_rasterinfo.SamplesStep.value())
    samples_step_unit = i_rasterinfo.SamplesStep.unit
    lines = i_rasterinfo.Lines
    lines_start_unit = i_rasterinfo.LinesStart.unit
    if lines_start_unit.lower() == 'utc':
        lines_start = PreciseDateTime()
        lines_start.set_from_utc_string(i_rasterinfo.LinesStart.value())
    else:
        lines_start = float(i_rasterinfo.LinesStart.value())
    lines_step = float(i_rasterinfo.LinesStep.value())
    lines_step_unit = i_rasterinfo.LinesStep.unit
    celltype = i_rasterinfo.CellType
    filename = i_rasterinfo.FileName
    header_offset_bytes = i_rasterinfo.HeaderOffsetBytes
    row_prefix_bytes = i_rasterinfo.RowPrefixBytes
    byteorder = i_rasterinfo.ByteOrder

    o_rasterinfo = metadata.RasterInfo(lines=lines, samples=samples, celltype=celltype,
                                       filename=filename,
                                       header_offset_bytes=header_offset_bytes,
                                       row_prefix_bytes=row_prefix_bytes,
                                       byteorder=byteorder)

    o_rasterinfo.set_lines_axis(lines_start, lines_start_unit, lines_step, lines_step_unit)
    o_rasterinfo.set_samples_axis(samples_start, samples_start_unit, samples_step, samples_step_unit)

    return o_rasterinfo


def convert_dataset_info_to_aresys(i_dataset_info: pyxb_metadata_types.DataSetInfoType):
    o_dataset_info = metadata.DataSetInfo(i_dataset_info.AcquisitionMode.value(), i_dataset_info.fc_hz.value())

    o_dataset_info.sensor_name = i_dataset_info.SensorName
    o_dataset_info.description = i_dataset_info.Description.value()
    sense_date = i_dataset_info.SenseDate.value()
    if sense_date == 'NOT_AVAILABLE' or sense_date is None:
        o_dataset_info.sense_date = None
    else:
        o_dataset_info.sense_date = \
            PreciseDateTime().set_from_utc_string(sense_date)
    o_dataset_info.image_type = i_dataset_info.ImageType.value()
    o_dataset_info.projection = i_dataset_info.Projection.value()
    o_dataset_info.acquisition_station = i_dataset_info.AcquisitionStation.value()
    o_dataset_info.processing_center = i_dataset_info.ProcessingCenter.value()
    process_date = i_dataset_info.ProcessingDate.value()
    if process_date == 'NOT_AVAILABLE' or process_date is None:
        o_dataset_info.sense_date = None
    else:
        o_dataset_info.processing_date = \
            PreciseDateTime().set_from_utc_string(process_date)
    o_dataset_info.processing_software = i_dataset_info.ProcessingSoftware.value()
    o_dataset_info.side_looking = i_dataset_info.SideLooking
    o_dataset_info.external_calibration_factor = i_dataset_info.ExternalCalibrationFactor
    o_dataset_info.data_take_id = i_dataset_info.DataTakeID

    return o_dataset_info


def convert_state_vectors_to_aresys(i_state_vectors: pyxb_metadata_types.StateVectorDataType):
    number_of_state_vectors = i_state_vectors.nSV_n.value()

    position_vector = np.zeros((number_of_state_vectors, 3))
    velocity_vector = np.zeros((number_of_state_vectors, 3))
    for ecef_coord_index, ecef_coord in enumerate(i_state_vectors.pSV_m.val):
        position_vector[ecef_coord_index // 3][ecef_coord_index % 3] = \
            ecef_coord.value()
        velocity_vector[ecef_coord_index // 3][ecef_coord_index % 3] = \
            i_state_vectors.vSV_mOs.val[ecef_coord_index].value()

    reference_time = PreciseDateTime()
    reference_time.set_from_utc_string(i_state_vectors.t_ref_Utc)
    delta_time = i_state_vectors.dtSV_s.value()

    o_state_vectors = metadata.StateVectors(position_vector,
                                            velocity_vector,
                                            reference_time,
                                            delta_time)

    if i_state_vectors.OrbitNumber != 'NOT_AVAILABLE':
        o_state_vectors.orbit_number = int(i_state_vectors.OrbitNumber)

    if i_state_vectors.Track != 'NOT_AVAILABLE':
        o_state_vectors.track_number = int(i_state_vectors.Track)

    return o_state_vectors


def convert_attitude_info_to_aresys(i_attitude_info: pyxb_metadata_types.AttitudeInfoType):
    o_attitude_info = metadata.AttitudeInfo(yaw=_pyxb_sequence_to_list(i_attitude_info.yaw_deg),
                                            pitch=_pyxb_sequence_to_list(i_attitude_info.pitch_deg),
                                            roll=_pyxb_sequence_to_list(i_attitude_info.roll_deg),
                                            t0=PreciseDateTime().set_from_utc_string(i_attitude_info.t_ref_Utc),
                                            delta_t=i_attitude_info.dtYPR_s.value(),
                                            ref_frame=i_attitude_info.referenceFrame,
                                            rot_order=i_attitude_info.rotationOrder)

    o_attitude_info.attitude_type = i_attitude_info.AttitudeType

    return o_attitude_info


def _pyxb_sequence_to_list(pyxb_seq):
    if pyxb_seq is None:
        return []
    return [v.value() for v in pyxb_seq.val]


def convert_acquisition_time_line_to_aresys(i_acquisition_time_line: pyxb_metadata_types.AcquisitionTimelineType):
    swl_changes_number = i_acquisition_time_line.Swl_changes_number if i_acquisition_time_line.Swl_changes_number else 0
    o_acquisition_time_line = metadata.AcquisitionTimeLine(
        missing_lines_number_i=i_acquisition_time_line.MissingLines_number.value(),
        missing_lines_azimuth_times_i=_pyxb_sequence_to_list(i_acquisition_time_line.MissingLines_azimuthtimes),
        swst_changes_number_i=i_acquisition_time_line.Swst_changes_number.value(),
        swst_changes_azimuth_times_i=_pyxb_sequence_to_list(i_acquisition_time_line.Swst_changes_azimuthtimes),
        swst_changes_values_i=_pyxb_sequence_to_list(i_acquisition_time_line.Swst_changes_values),
        noise_packets_number_i=i_acquisition_time_line.noise_packets_number,
        noise_packets_azimuth_times_i=_pyxb_sequence_to_list(i_acquisition_time_line.noise_packets_azimuthtimes),
        internal_calibration_number_i=i_acquisition_time_line.Internal_calibration_number,
        internal_calibration_azimuth_times_i=_pyxb_sequence_to_list(
            i_acquisition_time_line.Internal_calibration_azimuthtimes),
        swl_changes_number_i=swl_changes_number,
        swl_changes_azimuth_times_i=_pyxb_sequence_to_list(i_acquisition_time_line.Swl_changes_azimuthtimes),
        swl_changes_values_i=_pyxb_sequence_to_list(i_acquisition_time_line.Swl_changes_values))

    if i_acquisition_time_line.DuplicatedLines_number is not None:
        o_acquisition_time_line.duplicated_lines = (
            i_acquisition_time_line.DuplicatedLines_number.value(),
            _pyxb_sequence_to_list(i_acquisition_time_line.DuplicatedLines_azimuthtimes))

    return o_acquisition_time_line


def convert_geo_point_to_aresys(i_point: pyxb_metadata_types.PointType) -> metadata.GeoPoint:
    return metadata.GeoPoint(*_pyxb_sequence_to_list(i_point))


def convert_ground_corner_points_to_aresys(
        i_ground_corner_points: pyxb_metadata_types.GroundCornersPointsType) -> metadata.GroundCornerPoints:
    o_ground_corner_points = metadata.GroundCornerPoints()
    o_ground_corner_points.easting_grid_size = i_ground_corner_points.EastingGridSize.value()
    o_ground_corner_points.northing_grid_size = i_ground_corner_points.NorthingGridSize.value()
    o_ground_corner_points.center_point = convert_geo_point_to_aresys(i_ground_corner_points.Center.Point)
    o_ground_corner_points.nw_point = convert_geo_point_to_aresys(i_ground_corner_points.NorthWest.Point)
    o_ground_corner_points.ne_point = convert_geo_point_to_aresys(i_ground_corner_points.NorthEast.Point)
    o_ground_corner_points.sw_point = convert_geo_point_to_aresys(i_ground_corner_points.SouthWest.Point)
    o_ground_corner_points.se_point = convert_geo_point_to_aresys(i_ground_corner_points.SouthEast.Point)

    return o_ground_corner_points


def convert_burst_to_aresys(i_burst: pyxb_metadata_types.BurstType, lines) -> metadata.Burst:
    o_burst = metadata.Burst(i_burst.RangeStartTime.value(),
                             PreciseDateTime().set_from_utc_string(i_burst.AzimuthStartTime.value()), lines,
                             i_burst.BurstCenterAzimuthShift)
    return o_burst


def convert_antenna_info_to_aresys(i_antenna_info: pyxb_metadata_types.AntennaInfoType) -> metadata.AntennaInfo:
    o_antenna_info = metadata.AntennaInfo(i_antenna_info.SensorName,
                                   i_antenna_info.Polarization,
                                   i_antenna_info.AcquisitionMode,
                                   i_antenna_info.BeamName)

    if i_antenna_info.LinesPerPattern is not None:
        o_antenna_info.lines_per_pattern = i_antenna_info.LinesPerPattern

    return o_antenna_info


def convert_data_statistics_to_aresys(
        i_data_statistics: pyxb_metadata_types.DataStatisticsType) -> metadata.DataStatistics:
    o_data_statistics = metadata.DataStatistics(i_data_statistics.NumSamples.value(),
                                                i_data_statistics.MaxI.value(),
                                                i_data_statistics.MaxQ.value(),
                                                i_data_statistics.MinI.value(),
                                                i_data_statistics.MinQ.value(),
                                                i_data_statistics.SumI.value(),
                                                i_data_statistics.SumQ.value(),
                                                i_data_statistics.Sum2I.value(),
                                                i_data_statistics.Sum2Q.value(),
                                                i_data_statistics.StdDevI.value(),
                                                i_data_statistics.StdDevQ.value())

    if i_data_statistics.StatisticsList is not None:
        for data_block_stat_pyxb in i_data_statistics.StatisticsList.DataBlockStatistic:
            data_block_stat = metadata.DataBlockStatistic(data_block_stat_pyxb.lineStart,
                                                          data_block_stat_pyxb.lineStop,
                                                          data_block_stat_pyxb.NumSamples.value(),
                                                          data_block_stat_pyxb.MaxI.value(),
                                                          data_block_stat_pyxb.MaxQ.value(),
                                                          data_block_stat_pyxb.MinI.value(),
                                                          data_block_stat_pyxb.MinQ.value(),
                                                          data_block_stat_pyxb.SumI.value(),
                                                          data_block_stat_pyxb.SumQ.value(),
                                                          data_block_stat_pyxb.Sum2I.value(),
                                                          data_block_stat_pyxb.Sum2Q.value())
            o_data_statistics.add_data_block_statistic(data_block_stat)

    return o_data_statistics


def convert_poly_2d_to_aresys(i_poly_2d, i_poly_type=None):

    if i_poly_type is None:
        poly_type = metadata._Poly2DVector
    else:
        poly_type = i_poly_type

    o_poly_2d_vector = poly_type()
    for current_poly in i_poly_2d:
        coefficients = np.zeros(len(current_poly.pol.val))
        for index_pol, value in enumerate(current_poly.pol.val):
            coefficients[index_pol] = value.value()
        o_poly_2d_vector.add_poly(
            poly_type.get_single_poly_type()(current_poly.taz0_Utc.value(),
                                             current_poly.trg0_s.value(),
                                             coefficients))

    return o_poly_2d_vector


def convert_doppler_centroid_to_aresys(
    i_poly: pyxb_metadata_types.polyType) -> metadata.DopplerCentroidVector:
    return convert_poly_2d_to_aresys(i_poly, metadata.DopplerCentroidVector)


def convert_doppler_rate_to_aresys(
    i_poly: pyxb_metadata_types.polyType) -> metadata.DopplerRateVector:
    return convert_poly_2d_to_aresys(i_poly, metadata.DopplerRateVector)


def convert_tops_azimuth_modulation_rate_to_aresys(
    i_poly: pyxb_metadata_types.polyType) -> metadata.TopsAzimuthModulationRateVector:
    return convert_poly_2d_to_aresys(i_poly, metadata.TopsAzimuthModulationRateVector)


def convert_slant_to_ground_to_aresys(
    i_poly: pyxb_metadata_types.polyType) -> metadata.SlantToGroundVector:
    return convert_poly_2d_to_aresys(i_poly, metadata.SlantToGroundVector)


def convert_ground_to_slant_to_aresys(
    i_poly: pyxb_metadata_types.polyType) -> metadata.GroundToSlantVector:
    return convert_poly_2d_to_aresys(i_poly, metadata.GroundToSlantVector)


def convert_slant_to_incidence_to_aresys(
    i_poly: pyxb_metadata_types.polyType) -> metadata.SlantToIncidenceVector:
    return convert_poly_2d_to_aresys(i_poly, metadata.SlantToIncidenceVector)


def convert_slant_to_elevation_to_aresys(
    i_poly: pyxb_metadata_types.polyType) -> metadata.SlantToElevationVector:
    return convert_poly_2d_to_aresys(i_poly, metadata.SlantToElevationVector)


_COREG_POLY_COEF_PYXB_TO_ARESYS = {0: 4,
                                   1: 1,
                                   2: 2,
                                   3: 0,
                                   }


def convert_coreg_poly_to_aresys(
    i_coreg_poly: pyxb_metadata_types.polyCoregType) -> metadata.CoregPolyVector:
    o_coreg_poly_vector = metadata.CoregPolyVector()
    for current_poly in i_coreg_poly:
        coefficients_rg = np.zeros(max(_COREG_POLY_COEF_PYXB_TO_ARESYS.values())+1)
        coefficients_az = np.zeros(max(_COREG_POLY_COEF_PYXB_TO_ARESYS.values())+1)
        for index_pyxb, index_aresys in _COREG_POLY_COEF_PYXB_TO_ARESYS.items():
            coefficients_rg[index_aresys] = \
                current_poly.polRg.val[index_pyxb].value()
            coefficients_az[index_aresys] = \
                current_poly.polAz.val[index_pyxb].value()
        o_coreg_poly_vector.add_poly(
            metadata.CoregPoly(current_poly.taz0_Utc.value(),
                               current_poly.trg0_s.value(),
                               coefficients_az,
                               coefficients_rg))

    return o_coreg_poly_vector


def _expand_lines_per_burst_dict(lines_per_burst, num_bursts):
    out = list()
    for i_burst in range(num_bursts):
        for i_burst_change in sorted(lines_per_burst.keys(), reverse=True):
            if i_burst + 1 >= i_burst_change:
                out.append(lines_per_burst[i_burst_change])
                break
    return out


def _reduce_to_lines_per_burst_dict(lines_per_burst):
    out = dict()
    last_num_lines = 0
    for i_burst, lines_per_burst in enumerate(lines_per_burst):
        if lines_per_burst != last_num_lines:
            last_num_lines = lines_per_burst
            out[i_burst] = last_num_lines
    return out


def convert_burst_info_to_aresys(i_burst_info: pyxb_metadata_types.BurstInfoType) -> metadata.BurstInfo:
    number_of_bursts = i_burst_info.NumberOfBursts

    o_burst_info = metadata.BurstInfo(burst_repetition_frequency=i_burst_info.BurstRepetitionFrequency.value())
    if i_burst_info.LinesPerBurst is not None:
        lines_per_burst = [i_burst_info.LinesPerBurst for _ in range(number_of_bursts)]
    else:
        lines_per_burst = dict()
        for lines in i_burst_info.LinesPerBurstChangeList.Lines:
            lines_per_burst[lines.FromBurst] = lines.value()
        lines_per_burst = _expand_lines_per_burst_dict(lines_per_burst, number_of_bursts)

    for i_burst in range(number_of_bursts):
        burst = convert_burst_to_aresys(i_burst_info.Burst[i_burst], lines_per_burst[i_burst])
        o_burst_info.add_burst(burst.range_start_time, burst.azimuth_start_time, burst.lines, burst.burst_center_azimuth_shift)

    return o_burst_info


_SECTION_TO_CONVERTER_TO_ARESYS = {
    'SamplingConstants': convert_sampling_constants_to_aresys,
    'Pulse': convert_pulse_to_aresys,
    'SwathInfo': convert_swath_info_to_aresys,
    'RasterInfo': convert_rasterinfo_to_aresys,
    'DataSetInfo': convert_dataset_info_to_aresys,
    'StateVectorData': convert_state_vectors_to_aresys,
    'AttitudeInfo': convert_attitude_info_to_aresys,
    'AcquisitionTimeLine': convert_acquisition_time_line_to_aresys,
    'GroundCornerPoints': convert_ground_corner_points_to_aresys,
    'BurstInfo': convert_burst_info_to_aresys,
    'DopplerCentroid': convert_doppler_centroid_to_aresys,
    'DopplerRate': convert_doppler_rate_to_aresys,
    'TopsAzimuthModulationRate': convert_tops_azimuth_modulation_rate_to_aresys,
    'SlantToGround': convert_slant_to_ground_to_aresys,
    'GroundToSlant': convert_ground_to_slant_to_aresys,
    'SlantToIncidence': convert_slant_to_incidence_to_aresys,
    'SlantToElevation': convert_slant_to_elevation_to_aresys,
    'AntennaInfo': convert_antenna_info_to_aresys,
    'DataStatistics': convert_data_statistics_to_aresys,
    'CoregPoly': convert_coreg_poly_to_aresys,
}

_SECTION_ARESYS_TO_PYXB = {
    'SamplingConstants': 'SamplingConstants',
    'Pulse': 'Pulse',
    'SwathInfo': 'SwathInfo',
    'RasterInfo': 'RasterInfo',
    'DataSetInfo': 'DataSetInfo',
    'StateVectors': 'StateVectorData',
    'AttitudeInfo': 'AttitudeInfo',
    'AcquisitionTimeLine': 'AcquisitionTimeLine',
    'GroundCornerPoints': 'GroundCornerPoints',
    'BurstInfo': 'BurstInfo',
    'DopplerCentroidVector': 'DopplerCentroid',
    'DopplerRateVector': 'DopplerRate',
    'TopsAzimuthModulationRateVector': 'TopsAzimuthModulationRate',
    'SlantToGroundVector': 'SlantToGround',
    'GroundToSlantVector': 'GroundToSlant',
    'SlantToIncidenceVector': 'SlantToIncidence',
    'SlantToElevationVector': 'SlantToElevation',
    'AntennaInfo': 'AntennaInfo',
    'DataStatistics': 'DataStatistics',
    'CoregPolyVector': 'CoregPoly',
}

supported_elements = metadata.MetaDataChannel.get_supported_metadata_elements()
supported_elements.append('StateVectorData')  # in MetaDataChannel is StateVectors
readable_elements = list(_SECTION_TO_CONVERTER_TO_ARESYS.keys())
readable_elements.append('StateVectors')  # in MetaDataChannel is StateVectors
for supported_element in supported_elements:
    if supported_element not in readable_elements:
        logging.getLogger(__name__).warning("{} not readable".format(supported_element))
for readable_element in readable_elements:
    if readable_element not in supported_elements:
        logging.getLogger(__name__).warning("{} not supported".format(readable_element))


def run_conversion_to_aresys(metadata_i: pyxb_metadata.AresysXmlDocType):
    metadata_o = metadata.MetaData()

    for ch in metadata_i.Channel:
        mdc = metadata.MetaDataChannel()

        # RasterInfo
        if not (ch.RasterInfo is None):
            raster_info = convert_rasterinfo_to_aresys(ch.RasterInfo)
            mdc.insert_element(raster_info)
        else:
            logger = logging.getLogger(__name__)
            logger.debug("Raster info not available in Channel #{}".format(ch.Number))

        for section, converter in _SECTION_TO_CONVERTER_TO_ARESYS.items():
            if section != "RasterInfo":
                section_pyxb = getattr(ch, section)
                if section_pyxb is not None:
                    mdc.insert_element(converter(section_pyxb))

        metadata_o.append_channel(mdc)

    return metadata_o


def convert_sampling_constants_to_pyxb(
        i_sampling_constants: metadata.SamplingConstants) -> pyxb_metadata_types.SamplingConstantsType:
    o_sampling_constants = pyxb_metadata_types.SamplingConstantsType(frg_hz=i_sampling_constants.frg_hz,
                                                                     Brg_hz=i_sampling_constants.brg_hz,
                                                                     faz_hz=i_sampling_constants.faz_hz,
                                                                     Baz_hz=i_sampling_constants.baz_hz)

    o_sampling_constants.frg_hz.unit = str(i_sampling_constants.frg_hz_unit)
    o_sampling_constants.Brg_hz.unit = str(i_sampling_constants.brg_hz_unit)
    o_sampling_constants.faz_hz.unit = str(i_sampling_constants.faz_hz_unit)
    o_sampling_constants.Baz_hz.unit = str(i_sampling_constants.baz_hz_unit)

    return o_sampling_constants


def convert_pulse_to_pyxb(i_pulse: metadata.Pulse):
    o_pulse = pyxb_metadata_types.PulseType(PulseLength=i_pulse.pulse_length,
                                            Bandwidth=i_pulse.bandwidth,
                                            PulseSamplingRate=i_pulse.pulse_sampling_rate,
                                            PulseEnergy=i_pulse.pulse_energy)

    if i_pulse.pulse_direction is not None:
        o_pulse.Direction = i_pulse.pulse_direction.value
    if i_pulse.pulse_start_frequency is not None:
        o_pulse.PulseStartFrequency = i_pulse.pulse_start_frequency
        o_pulse.PulseStartFrequency.unit = str(i_pulse.pulse_start_frequency_unit)
    if i_pulse.pulse_start_phase is not None:
        o_pulse.PulseStartPhase = i_pulse.pulse_start_phase
        o_pulse.PulseStartPhase.unit = str(i_pulse.pulse_start_phase_unit)

    o_pulse.PulseLength.unit = str(i_pulse.pulse_length_unit)
    o_pulse.Bandwidth.unit = str(i_pulse.bandwidth_unit)
    o_pulse.PulseEnergy.unit = str(i_pulse.pulse_energy_unit)
    o_pulse.PulseSamplingRate.unit = str(i_pulse.pulse_sampling_rate_unit)

    return o_pulse


def convert_swath_info_to_pyxb(i_swath_info: metadata.SwathInfo) -> pyxb_metadata_types.SwathInfoType:
    o_swath_info = pyxb_metadata_types.SwathInfoType(Swath=i_swath_info.swath,
                                                     SwathAcquisitionOrder=i_swath_info.swath_acquisition_order,
                                                     Polarization=i_swath_info.polarization.value,
                                                     Rank=i_swath_info.rank,
                                                     RangeDelayBias=i_swath_info.range_delay_bias,
                                                     AcquisitionStartTime=str(i_swath_info.acquisition_start_time),
                                                     AzimuthSteeringRateReferenceTime=
                                                     i_swath_info.azimuth_steering_rate_reference_time,
                                                     AzimuthSteeringRatePol=pyxb.BIND(),
                                                     AcquisitionPRF=i_swath_info.acquisition_prf,
                                                     EchoesPerBurst=i_swath_info.echoes_per_burst,
                                                     )

    for index, val in enumerate(i_swath_info.azimuth_steering_rate_pol):
        o_swath_info.AzimuthSteeringRatePol.append(val)
        o_swath_info.AzimuthSteeringRatePol.val[index].N = index+1

    if i_swath_info.channel_delay is not None:
        o_swath_info.ChannelDelay = i_swath_info.channel_delay
    if i_swath_info.rx_gain is not None:
        o_swath_info.RxGain = i_swath_info.rx_gain

    o_swath_info.AcquisitionPRF.unit = str(i_swath_info.acquisition_prf_unit)
    o_swath_info.RangeDelayBias.unit = str(i_swath_info.range_delay_bias_unit)
    o_swath_info.AcquisitionStartTime.unit = str(i_swath_info.acquisition_start_time_unit)
    o_swath_info.AzimuthSteeringRateReferenceTime.unit = str(i_swath_info.az_steering_rate_ref_time_unit)

    return o_swath_info


def convert_rasterinfo_to_pyxb(i_rasterinfo: metadata.RasterInfo):
    o_rasterinfo = pyxb_metadata_types.RasterInfoType(Lines=i_rasterinfo.lines,
                                                      LinesStep=i_rasterinfo.lines_step,
                                                      LinesStart=str(i_rasterinfo.lines_start),
                                                      Samples=i_rasterinfo.samples,
                                                      SamplesStep=i_rasterinfo.samples_step,
                                                      SamplesStart=str(i_rasterinfo.samples_start),
                                                      HeaderOffsetBytes=i_rasterinfo.header_offset_bytes,
                                                      RowPrefixBytes=i_rasterinfo.row_prefix_bytes,
                                                      FileName=i_rasterinfo.file_name,
                                                      ByteOrder=i_rasterinfo.byte_order.value,
                                                      CellType=i_rasterinfo.cell_type.value)

    o_rasterinfo.SamplesStart.unit = str(i_rasterinfo.samples_start_unit)
    o_rasterinfo.SamplesStep.unit = str(i_rasterinfo.samples_step_unit)
    o_rasterinfo.LinesStart.unit = str(i_rasterinfo.lines_start_unit)
    o_rasterinfo.LinesStep.unit = str(i_rasterinfo.lines_step_unit)

    return o_rasterinfo


def convert_dataset_info_to_pyxb(i_dataset_info: metadata.DataSetInfo) -> pyxb_metadata_types.DataSetInfoType:
    if i_dataset_info.sense_date is not None:
        sense_date_str = str(i_dataset_info.sense_date)
    else:
        sense_date_str = 'NOT_AVAILABLE'

    if i_dataset_info.processing_date is not None:
        processing_date_str = str(i_dataset_info.processing_date)
    else:
        processing_date_str = 'NOT_AVAILABLE'

    o_dataset_info = pyxb_metadata_types.DataSetInfoType(SensorName=i_dataset_info.sensor_name,
                                                         Description=i_dataset_info.description,
                                                         SenseDate=sense_date_str,
                                                         AcquisitionMode=i_dataset_info.acquisition_mode,
                                                         ImageType=i_dataset_info.image_type,
                                                         Projection=i_dataset_info.projection,
                                                         AcquisitionStation=i_dataset_info.acquisition_station,
                                                         ProcessingCenter=i_dataset_info.processing_center,
                                                         ProcessingDate=processing_date_str,
                                                         ProcessingSoftware=i_dataset_info.processing_software,
                                                         fc_hz=i_dataset_info.fc_hz,
                                                         SideLooking=i_dataset_info.side_looking.value,
                                                         ExternalCalibrationFactor=
                                                         i_dataset_info.external_calibration_factor,
                                                         DataTakeID=i_dataset_info.data_take_id)
    return o_dataset_info


def convert_state_vectors_to_pyxb(i_state_vectors: metadata.StateVectors) -> pyxb_metadata_types.StateVectorDataType:
    o_state_vectors = pyxb_metadata_types.StateVectorDataType(
        pSV_m=pyxb.BIND(),
        vSV_mOs=pyxb.BIND(),
        OrbitNumber=str(i_state_vectors.orbit_number),
        Track=str(i_state_vectors.track_number),
        OrbitDirection=i_state_vectors.orbit_direction.value,
        t_ref_Utc=str(i_state_vectors.reference_time),
        dtSV_s=str(i_state_vectors.time_step),
        nSV_n=i_state_vectors.number_of_state_vectors)

    for index_position, position in enumerate(i_state_vectors.position_vector):
        for index_sv_vect, coord in enumerate(position):
            index_current = index_position * 3 + index_sv_vect + 1

            o_state_vectors.pSV_m.append(coord)
            o_state_vectors.pSV_m.val[index_current - 1].N = str(index_current)

            o_state_vectors.vSV_mOs.append(i_state_vectors.velocity_vector[index_position][index_sv_vect])
            o_state_vectors.vSV_mOs.val[index_current - 1].N = str(index_current)

    return o_state_vectors


def _list_to_pyxb_sequence(input_list, pyxb_obj):
    if input_list is None:
        return pyxb_obj
    if len(input_list) > 0:
        for i_val, val in enumerate(input_list):
            pyxb_obj.append(val)
            pyxb_obj.val[i_val].N = i_val + 1


def convert_attitude_info_to_pyxb(i_attitude_info: metadata.AttitudeInfo) -> pyxb_metadata_types.AttitudeInfoType:
    o_attitude_info = pyxb_metadata_types.AttitudeInfoType(
        yaw_deg=pyxb.BIND(),
        pitch_deg=pyxb.BIND(),
        roll_deg=pyxb.BIND(),
        t_ref_Utc=str(i_attitude_info.reference_time),
        dtYPR_s=str(i_attitude_info.time_step),
        nYPR_n=i_attitude_info.attitude_records_number,
        referenceFrame=i_attitude_info.reference_frame.value.upper(),
        rotationOrder=i_attitude_info.rotation_order.value.upper(),
        AttitudeType=i_attitude_info.attitude_type.value.upper())

    _list_to_pyxb_sequence(i_attitude_info.yaw_vector, o_attitude_info.yaw_deg)
    _list_to_pyxb_sequence(i_attitude_info.pitch_vector, o_attitude_info.pitch_deg),
    _list_to_pyxb_sequence(i_attitude_info.roll_vector, o_attitude_info.roll_deg),

    return o_attitude_info


def convert_acquisition_time_line_to_pyxb(
        i_acquisition_time_line: metadata.AcquisitionTimeLine) -> pyxb_metadata_types.AcquisitionTimelineType:
    o_acquisition_time_line = pyxb_metadata_types.AcquisitionTimelineType(
        MissingLines_number=i_acquisition_time_line.missing_lines[0],
        MissingLines_azimuthtimes=pyxb.BIND(),
        Swst_changes_number=i_acquisition_time_line.swst_changes[0],
        Swst_changes_azimuthtimes=pyxb.BIND(),
        Swst_changes_values=pyxb.BIND(),
        noise_packets_number=i_acquisition_time_line.noise_packet[0],
        noise_packets_azimuthtimes=pyxb.BIND(),
        Internal_calibration_number=i_acquisition_time_line.internal_calibration[0],
        Internal_calibration_azimuthtimes=pyxb.BIND(),
        Swl_changes_number=i_acquisition_time_line.swl_changes[0],
        Swl_changes_azimuthtimes=pyxb.BIND(),
        Swl_changes_values=pyxb.BIND(),
        DuplicatedLines_number=i_acquisition_time_line.duplicated_lines[0],
        DuplicatedLines_azimuthtimes=pyxb.BIND())

    _list_to_pyxb_sequence(i_acquisition_time_line.missing_lines[1], o_acquisition_time_line.MissingLines_azimuthtimes)
    _list_to_pyxb_sequence(i_acquisition_time_line.swst_changes[1], o_acquisition_time_line.Swst_changes_azimuthtimes)
    _list_to_pyxb_sequence(i_acquisition_time_line.swst_changes[2], o_acquisition_time_line.Swst_changes_values)
    _list_to_pyxb_sequence(i_acquisition_time_line.noise_packet[1], o_acquisition_time_line.noise_packets_azimuthtimes)
    _list_to_pyxb_sequence(i_acquisition_time_line.internal_calibration[1],
                           o_acquisition_time_line.Internal_calibration_azimuthtimes)
    _list_to_pyxb_sequence(i_acquisition_time_line.swl_changes[1], o_acquisition_time_line.Swl_changes_azimuthtimes)
    _list_to_pyxb_sequence(i_acquisition_time_line.swl_changes[2], o_acquisition_time_line.Swl_changes_values)
    _list_to_pyxb_sequence(i_acquisition_time_line.duplicated_lines[1],
                           o_acquisition_time_line.DuplicatedLines_azimuthtimes)

    return o_acquisition_time_line


def convert_geo_point_to_pyxb(i_geo_point: metadata.GeoPoint) -> pyxb_metadata_types.PointType:
    o_point = pyxb_metadata_types.PointType()

    for i_val, val in enumerate(i_geo_point.to_list()):
        o_point.append(val)
        o_point.val[i_val].N = i_val + 1

    return o_point


def convert_ground_corner_points_to_pyxb(
        i_ground_corner_points: metadata.GroundCornerPoints) -> pyxb_metadata_types.GroundCornersPointsType:
    o_ground_corner_points = pyxb_metadata_types.GroundCornersPointsType(
        EastingGridSize=i_ground_corner_points.easting_grid_size,
        NorthingGridSize=i_ground_corner_points.northing_grid_size,
        NorthWest=pyxb.BIND(convert_geo_point_to_pyxb(i_ground_corner_points.nw_point)),
        NorthEast=pyxb.BIND(convert_geo_point_to_pyxb(i_ground_corner_points.ne_point)),
        SouthWest=pyxb.BIND(convert_geo_point_to_pyxb(i_ground_corner_points.sw_point)),
        SouthEast=pyxb.BIND(convert_geo_point_to_pyxb(i_ground_corner_points.se_point)),
        Center=pyxb.BIND(convert_geo_point_to_pyxb(i_ground_corner_points.center_point))
    )
    return o_ground_corner_points


def convert_burst_to_pyxb(i_burst: metadata.Burst) -> pyxb_metadata_types.BurstType:
    o_burst = pyxb_metadata_types.BurstType(RangeStartTime=i_burst.range_start_time,
                                            AzimuthStartTime=str(i_burst.azimuth_start_time),
                                            BurstCenterAzimuthShift=i_burst.burst_center_azimuth_shift)
    o_burst.RangeStartTime.unit = constants.SECOND_STR
    o_burst.AzimuthStartTime.unit = constants.UTC_STR
    if o_burst.BurstCenterAzimuthShift is not None:
        o_burst.BurstCenterAzimuthShift.unit = constants.SECOND_STR
    return o_burst


def convert_burst_info_to_pyxb(i_burst_info: metadata.BurstInfo) -> pyxb_metadata_types.BurstInfoType:
    o_burst_info = pyxb_metadata_types.BurstInfoType(NumberOfBursts=i_burst_info.get_number_of_bursts(),
                                                     BurstRepetitionFrequency=i_burst_info.burst_repetition_frequency)

    o_burst_info.BurstRepetitionFrequency.unit = constants.HERTZ_STR

    lines_per_burst = list()
    for i_burst in range(i_burst_info.get_number_of_bursts()):
        lines_per_burst.append(i_burst_info.get_lines(i_burst))
    from_burst_dict = _reduce_to_lines_per_burst_dict(lines_per_burst)

    o_burst_info.LinesPerBurstChangeList = pyxb.BIND()

    for i_from_burst, (from_burst, lines) in enumerate(from_burst_dict.items()):
        o_burst_info.LinesPerBurstChangeList.Lines.append(lines)
        o_burst_info.LinesPerBurstChangeList.Lines[i_from_burst].FromBurst = from_burst + 1

    for i_burst in range(i_burst_info.get_number_of_bursts()):
        burst = i_burst_info.get_burst(i_burst)
        o_burst_info.Burst.append(convert_burst_to_pyxb(burst))
        o_burst_info.Burst[i_burst].N = i_burst + 1

    return o_burst_info


def convert_antenna_info_to_pyxb(i_antenna_info: metadata.AntennaInfo) -> pyxb_metadata_types.AntennaInfoType:
    o_antenna_info = pyxb_metadata_types.AntennaInfoType(SensorName=i_antenna_info.sensor_name,
                                                         AcquisitionMode=i_antenna_info.acquisition_mode,
                                                         BeamName=i_antenna_info.acquisition_beam,
                                                         Polarization=i_antenna_info.polarization.value,
                                                         )
    if i_antenna_info.lines_per_pattern != 0:
        o_antenna_info.LinesPerPattern = i_antenna_info.lines_per_pattern

    return o_antenna_info


def convert_data_block_statistics_to_pyxb(
        i_data_block_statistics: metadata.DataBlockStatistic) -> pyxb_metadata_types.DataBlockStatisticsType:
    o_data_block_statistics = pyxb_metadata_types.DataBlockStatisticsType(
        NumSamples=i_data_block_statistics.num_samples,
        MaxI=i_data_block_statistics.max_i,
        MinI=i_data_block_statistics.min_i,
        MaxQ=i_data_block_statistics.max_q,
        MinQ=i_data_block_statistics.min_q,
        SumI=i_data_block_statistics.sum_i,
        SumQ=i_data_block_statistics.sum_q,
        Sum2I=i_data_block_statistics.sum_2_i,
        Sum2Q=i_data_block_statistics.sum_2_q,
        )
    o_data_block_statistics.lineStart = i_data_block_statistics.line_start
    o_data_block_statistics.lineStop = i_data_block_statistics.line_stop

    return o_data_block_statistics


def convert_data_statistics_to_pyxb(
        i_data_statistics: metadata.DataStatistics) -> pyxb_metadata_types.DataStatisticsType:
    o_data_block_statistics = pyxb_metadata_types.DataStatisticsType(NumSamples=i_data_statistics.num_samples,
                                                                     MaxI=i_data_statistics.max_i,
                                                                     MinI=i_data_statistics.min_i,
                                                                     MaxQ=i_data_statistics.max_q,
                                                                     MinQ=i_data_statistics.min_q,
                                                                     SumI=i_data_statistics.sum_i,
                                                                     SumQ=i_data_statistics.sum_q,
                                                                     Sum2I=i_data_statistics.sum_2_i,
                                                                     Sum2Q=i_data_statistics.sum_2_q,
                                                                     StdDevI=i_data_statistics.std_dev_i,
                                                                     StdDevQ=i_data_statistics.std_dev_q,
                                                                     )
    if i_data_statistics.get_number_of_data_block_statistic() != 0:
        for index in range(0, i_data_statistics.get_number_of_data_block_statistic()-1):
            o_data_block_statistics.StatisticsList.DataBlockStatistic[index] = \
                convert_data_block_statistics_to_pyxb(i_data_statistics.get_data_block_statistic(index))

    return o_data_block_statistics


def convert_poly_2d_to_pyxb(i_poly_2d_vector: metadata._Poly2DVector) -> pyxb_metadata_types.polyType:
    o_poly_2d = list()
    units = metadata._Poly2D.get_units()

    for index_poly, poly_2d in enumerate(i_poly_2d_vector):
        current_poly = pyxb_metadata_types.polyType(trg0_s=poly_2d.t_ref_rg,
                                                      taz0_Utc=poly_2d.t_ref_az,
                                                      pol=pyxb.BIND())
        for index_coef, coef in enumerate(poly_2d.coefficients):
            current_poly.pol.val.append(coef)
            current_poly.pol.val[index_coef].N = index_coef + 1
            current_poly.pol.val[index_coef].unit = units[index_coef]

        current_poly.Number = index_poly + 1
        current_poly.Total = i_poly_2d_vector.get_number_of_poly()

        o_poly_2d.append(current_poly)

    return o_poly_2d


def convert_coreg_poly_to_pyxb(
    i_coreg_poly_vector: metadata.CoregPolyVector) -> pyxb_metadata_types.polyCoregType:
    o_coreg_poly = list()

    for index_poly, coreg_poly in enumerate(i_coreg_poly_vector):
        current_poly = pyxb_metadata_types.polyCoregType(trg0_s=coreg_poly.ref_range_time,
                                                         taz0_Utc=coreg_poly.ref_azimuth_time,
                                                         polRg=pyxb.BIND(),
                                                         polAz=pyxb.BIND())
        index_coef = 0
        for index_pyxb, index_aresys in _COREG_POLY_COEF_PYXB_TO_ARESYS.items():
            current_poly.polAz.val.append(coreg_poly.azimuth_poly.coefficients[index_aresys])
            current_poly.polRg.val.append(coreg_poly.range_poly.coefficients[index_aresys])
            current_poly.polAz.val[index_coef].N = index_coef + 1
            current_poly.polRg.val[index_coef].N = index_coef + 1
            index_coef += 1

        current_poly.Number = index_poly + 1
        current_poly.Total = i_coreg_poly_vector.get_number_of_poly()

        o_coreg_poly.append(current_poly)

    return o_coreg_poly


_SECTION_TO_CONVERTER_TO_PYXB = {
    'SamplingConstants': convert_sampling_constants_to_pyxb,
    'Pulse': convert_pulse_to_pyxb,
    'SwathInfo': convert_swath_info_to_pyxb,
    'RasterInfo': convert_rasterinfo_to_pyxb,
    'DataSetInfo': convert_dataset_info_to_pyxb,
    'StateVectors': convert_state_vectors_to_pyxb,
    'AttitudeInfo': convert_attitude_info_to_pyxb,
    'AcquisitionTimeLine': convert_acquisition_time_line_to_pyxb,
    'GroundCornerPoints': convert_ground_corner_points_to_pyxb,
    'BurstInfo': convert_burst_info_to_pyxb,
    'DopplerCentroidVector': convert_poly_2d_to_pyxb,
    'DopplerRateVector': convert_poly_2d_to_pyxb,
    'TopsAzimuthModulationRateVector': convert_poly_2d_to_pyxb,
    'SlantToGroundVector': convert_poly_2d_to_pyxb,
    'GroundToSlantVector': convert_poly_2d_to_pyxb,
    'SlantToIncidenceVector': convert_poly_2d_to_pyxb,
    'SlantToElevationVector': convert_poly_2d_to_pyxb,
    'AntennaInfo': convert_antenna_info_to_pyxb,
    'DataStatistics': convert_data_statistics_to_pyxb,
    'CoregPolyVector': convert_coreg_poly_to_pyxb,
}

supported_elements = metadata.MetaDataChannel.get_supported_metadata_elements()
writable_elements = list(_SECTION_TO_CONVERTER_TO_PYXB.keys())
for supported_element in supported_elements:
    if supported_element not in writable_elements:
        logging.getLogger(__name__).warning("{} not writable".format(supported_element))
for writable_element in writable_elements:
    if writable_element not in supported_elements:
        logging.getLogger(__name__).warning("{} not supported".format(writable_element))


def run_conversion_pyxb(metadata_i: metadata.MetaData) -> pyxb_metadata.AresysXmlDocType:
    version = 2.1
    o_metadata = pyxb_metadata.AresysXmlDocType(metadata_i.get_number_of_channels(), version, '')

    for index_channel in range(0, metadata_i.get_number_of_channels()):
        o_metadata.Channel.append(pyxb.BIND())
        o_metadata.Channel[index_channel].Number = index_channel + 1
        o_metadata.Channel[index_channel].Total = metadata_i.get_number_of_channels()

        # RasterInfo is mandatory
        raster_info_aresys = metadata_i.get_metadata_channels(index_channel).get_element('RasterInfo')
        o_metadata.Channel[index_channel].RasterInfo = convert_rasterinfo_to_pyxb(raster_info_aresys)

        # Other sections
        for section, converter in _SECTION_TO_CONVERTER_TO_PYXB.items():
            if section != "RasterInfo":
                section_aresys = metadata_i.get_metadata_channels(index_channel).get_element(section)
                if section_aresys is not None:
                    setattr(o_metadata.Channel[index_channel], _SECTION_ARESYS_TO_PYXB[section], converter(section_aresys))

    return o_metadata
