# SPDX-FileCopyrightText: BioPAL <biopal@esa.int>
# SPDX-License-Identifier: MIT

import numpy as np
import logging


def check_pm_kz0(kz_stack, opt_str):

    # convert input kz from dict to list in order to perform ne max mmin search:
    kz_list = []
    [kz_list.append(v) for v in kz_stack.values()]
    del kz_stack

    condition = np.nansum(
        (np.nanmax(kz_list, axis=0) >= opt_str.kz0) & (np.nanmin(kz_list, axis=0) <= opt_str.kz0)
    ) < np.nansum((np.nanmax(kz_list, axis=0) >= -opt_str.kz0) & (np.nanmin(kz_list, axis=0) <= -opt_str.kz0))

    if condition:
        return -opt_str.kz0, kz_list

    else:
        return opt_str.kz0, kz_list


def kz0_crit(kz_stack, opt_str):

    # Check if it is better to generate +kz0 or -kz0
    opt_str.kz0, kz_list = check_pm_kz0(kz_stack, opt_str)

    # get the minimums and maximums of the modified kz stack, return as a 2D array
    min_array = np.nanmin(kz_list, axis=0)
    max_array = np.nanmax(kz_list, axis=0)

    # check for how many pixels the kz0 is in between
    num_array = (min_array < opt_str.kz0) & (opt_str.kz0 < max_array)
    pixelcount = np.sum(num_array)

    return pixelcount


def ground_cancellation_core(data_stack, pol_name, kz_stack, opt_str):
    """
     It computes the ground cancellation from input data stack, for the desired polarization.
    
     INPUT
          data_stack: stack of calibrated, ground steered slc images.
                      It is a dictionary of two nested dictionaries where:
                      each data_stack[ acquisition_name_string ][ polarization_name_string ] is an array of shape [Nrg x Naz]
          kz_stack:   array of phase-to-height conversion
                      factors (kz) OR [num_acq x 1] vector of kz for constant geometries. Needed if num_acq > 2.
                      It is a dictionary where kz_stack[acquisition_name_string] is an array of shape [Nrg x Naz]
          pol_name:   string, polarization to be processed from the ones into the data_stack dictionary
          opt_str.
                      kz0: desired phase-to-height
                      master: index of the master image
                      z_demod: (optional. Default: 0) vertical spectrum
                               demodulation height to perform interpolation
                      model_sign: (optional. Default: +1) model for the
                                  interferometric phase: model_sign*kz*z
                      z_demod_map: (optional. [Nrg x Naz] numpy array) vertical
                                   spectrum demodulation height to perform
                                   interpolation
                      kz0_map: (optional. [Nrg x Naz] numpy array) desired
                              phase-to-height map
                      space_varying_cancellation: (optional. Default: false)
                               if true opt_str.z_demod_map and opt_str.kz0_map are used
    
     OUTPUTkz0_crit
           GroundNotchedSLC: [Nrg x Naz] ground notche slc
    
     DEPENDENCIES
                 kzInterp
    """

    num_acq = len(data_stack)
    acq_names = list(data_stack.keys())
    first_acq_dict = data_stack[acq_names[0]]
    Nrg, Naz = first_acq_dict[pol_name].shape

    # Checking whether space-varying cancellation should run. It has already
    # been checked by ground_cancellation method, trust it.
    if hasattr(opt_str, "space_varying_cancellation") and opt_str.space_varying_cancellation:
        logging.info("AGB: (core) opt_str.space_varying_cancellation = True...:")
        space_varying_cancellation = True
    else:
        logging.info("AGB: (core) opt_str.space_varying_cancellation = False...:")
        space_varying_cancellation = False

    if num_acq == 2:

        master_acq_id = acq_names[0]
        slave_acq_id = acq_names[1]

        GroundNotchedSLC = data_stack[master_acq_id][pol_name] - data_stack[slave_acq_id][pol_name]
        mask_extrap = np.zeros((Nrg, Naz), dtype=np.int8)

    else:
        # Checking opt_str
        if not hasattr(opt_str, "kz0"):
            logging.error("Ground cancellation module: missing opt_str.kz0 input.")
            raise

        if not hasattr(opt_str, "z_demod"):
            logging.error("Ground cancellation module: missing opt_str.z_demod input.")
            raise

        if not hasattr(opt_str, "master_id"):
            logging.error("Ground cancellation module: missing opt_str.master_id input.")
            raise

        # Check if it is better to generate +kz0 or -kz0
        # Just one single kz0 to perform this check even in case of space-varying cancellation
        opt_str.kz0, _ = check_pm_kz0(kz_stack, opt_str)

        # Demodulation
        for acq_id in data_stack.keys():
            if space_varying_cancellation:
                data_stack[acq_id][pol_name] = data_stack[acq_id][pol_name] * np.exp(
                    -1j * kz_stack[acq_id] * opt_str.z_demod_map
                )
            else:
                data_stack[acq_id][pol_name] = data_stack[acq_id][pol_name] * np.exp(
                    -1j * kz_stack[acq_id] * opt_str.z_demod
                )

        # Generating synthetic SLC through interpolation
        if space_varying_cancellation:
            Ikz0, mask_extrap = kzInterp(data_stack, kz_stack, opt_str.kz0_map, pol_name)
        else:
            Ikz0, mask_extrap = kzInterp(data_stack, kz_stack, opt_str.kz0, pol_name)

        # Modulation
        for acq_id in data_stack.keys():
            if space_varying_cancellation:
                data_stack[acq_id][pol_name] = data_stack[acq_id][pol_name] * np.exp(
                    1j * kz_stack[acq_id] * opt_str.z_demod_map
                )
            else:
                data_stack[acq_id][pol_name] = data_stack[acq_id][pol_name] * np.exp(
                    1j * kz_stack[acq_id] * opt_str.z_demod
                )

        if space_varying_cancellation:
            Ikz0 = Ikz0 * np.exp(1j * opt_str.kz0 * opt_str.z_demod_map)
        else:
            Ikz0 = Ikz0 * np.exp(1j * opt_str.kz0 * opt_str.z_demod)

        GroundNotchedSLC = data_stack[opt_str.master_id][pol_name] - Ikz0

    return GroundNotchedSLC, mask_extrap


def ground_cancellation(data_stack, kz_stack, multi_master_flag, z_emph, eq_flag, pho_r, theta_look, ground_slope):

    #      data_stack: stack of calibrated, ground steered slc images.
    #                  It is a dictionary of two nested dictionaries where:
    #                    each data_stack[ acquisition_name_string ][ polarization_name_string ] is an array of shape [Nrg x Naz]
    #      kz_stack:   array of phase-to-height conversion
    #                  factors (kz) OR [num_acq x 1] vector of kz for constant geometries. Needed if num_acq > 2.
    #                  It is a dictionary where kz_stack[acquisition_name_string] is an array of shape [Nrg x Naz]
    if not eq_flag == "1" and not eq_flag == "2" and not eq_flag == "3":
        error_str = 'Configuration flag "ModelBasedEqualization" value "{}" not valid, choose among "1", "2" or "3", where "1"->always OFF; "2"->always ON; "3"->ON only if two acquisitions are present'.format(
            eq_flag
        )
        logging.error(error_str)
        raise

    # Nrg, Naz, num_acq = data_stack[0].shape

    # Initializing local paramters
    num_acq = len(data_stack)
    acq_names = list(data_stack.keys())
    first_acq_dict = data_stack[acq_names[0]]
    pol_names = list(first_acq_dict.keys())
    Nrg, Naz = first_acq_dict[pol_names[0]].shape

    # Ground notching
    GroundNotchedSLC = {}

    # Desired elevation for the peak of the ground notching processing
    class opt_str:
        pass

    # Checking whether z_emph is a suitable map to drive the ground cancellation algorithm
    opt_str.space_varying_cancellation = is_suitable_z_emph_map(z_emph, Nrg, Naz)

    # If it is not a suitable map and it is not a scalar, an exception should be raised

    # Ground cancellation parameters
    if opt_str.space_varying_cancellation:
        logging.info("AGB: (ground_cancellation) z_emph suitable 2D map...:")
        opt_str.z_demod_map = z_emph / 2  # [m]
        opt_str.kz0_map = np.pi / opt_str.z_demod_map / 2
        opt_str.z_demod = np.mean(opt_str.z_demod_map)
        opt_str.kz0 = np.pi / opt_str.z_demod / 2
    else:
        logging.info("AGB: (ground_cancellation) z_emph not suitable 2D map...:")
        opt_str.z_demod = z_emph / 2  # [m]
        opt_str.kz0 = np.pi / opt_str.z_demod / 2

    # Polarimetric channel
    for pol_name in pol_names:

        if num_acq == 2:
            opt_str.master_id = acq_names[0]
            notch_final, mask_final = ground_cancellation_core(data_stack, pol_name, kz_stack, opt_str)
            GroundNotchedSLC[pol_name] = notch_final

        elif not multi_master_flag:

            # choose optimal master from a stack. This is based on whether the kz0 is
            # within the kz interval of the respective master-slaves combination
            num_pixels_kz = {}

            for acq_id in kz_stack.keys():
                TEMP = kz_stack[acq_id].reshape((Nrg, Naz))
                opt_str.master_id = acq_id

                kz_stack_mod = {}
                for acq_id_to_mod in kz_stack.keys():
                    kz_stack_mod[acq_id_to_mod] = kz_stack[acq_id_to_mod] - TEMP

                # check the number of pixels for which the kz0 criterion holds
                pixelcount = kz0_crit(kz_stack_mod, opt_str)
                num_pixels_kz[acq_id] = pixelcount
                print(
                    "For Master {} of pol {} the number of pixels that satisfy the criterion is {}".format(
                        acq_id, pol_name, pixelcount
                    )
                )

            # get the index of the max
            acq_id_max_pixels = max(num_pixels_kz, key=num_pixels_kz.get)
            print("Choosing {} as optimal master SLC for pol {}".format(acq_id_max_pixels, pol_name))

            # compute ground notch
            opt_str.master_id = acq_id_max_pixels
            TEMP = kz_stack[opt_str.master_id].reshape((Nrg, Naz))
            kz_stack_mod = {}
            for acq_id_to_mod in kz_stack.keys():
                kz_stack_mod[acq_id_to_mod] = kz_stack[acq_id_to_mod] - TEMP

            notch_final, mask_final = ground_cancellation_core(data_stack, pol_name, kz_stack_mod, opt_str)
            notch_final[mask_final == 1] = 0

            GroundNotchedSLC[pol_name] = notch_final

        else:  # multi-master

            notch_final = np.zeros((Nrg, Naz))
            mask_final = np.zeros((Nrg, Naz))

            for acq_id in kz_stack.keys():
                TEMP = kz_stack[acq_id].reshape((Nrg, Naz))
                opt_str.master_id = acq_id

                kz_stack_mod = {}
                for acq_id_to_mod in kz_stack.keys():
                    kz_stack_mod[acq_id_to_mod] = kz_stack[acq_id_to_mod] - TEMP

                notch_temp, notch_mask = ground_cancellation_core(data_stack, pol_name, kz_stack_mod, opt_str)
                notch_temp[notch_mask == 1] = 0
                notch_final = notch_final + np.abs(notch_temp) ** 2
                mask_final = mask_final + 1 - notch_mask
            mask_final[mask_final == 0] = 1
            GroundNotchedSLC[pol_name] = np.sqrt(notch_final / mask_final)

        if eq_flag == "2":  # 2 = always active

            if num_acq == 2:
                normalization_factor = SingleBaselineEqualization(
                    z_emph, pho_r, kz_stack[acq_names[1]], theta_look, ground_slope
                )
                GroundNotchedSLC[pol_name] = GroundNotchedSLC[pol_name] / np.sqrt(normalization_factor)
            else:
                normalization_factor = SingleBaselineEqualization(z_emph, pho_r, opt_str.kz0, theta_look, ground_slope)
                GroundNotchedSLC[pol_name] = GroundNotchedSLC[pol_name] / np.sqrt(normalization_factor)

        elif eq_flag == "3" and num_acq == 2:  # 3 = active only when nyum_acq == 2
            normalization_factor = SingleBaselineEqualization(
                z_emph, pho_r, kz_stack[acq_names[1]], theta_look, ground_slope
            )
            GroundNotchedSLC[pol_name] = GroundNotchedSLC[pol_name] / np.sqrt(normalization_factor)

    return GroundNotchedSLC


def SingleBaselineEqualization(eq_height, pho_r, kz, theta_look, ground_slope):
    """It generates an model based equalization map for ground notched data from geometry.

    INPUT
      eq_height: the desired height for which the power is equalized
      pho_r: slant range resolution
      kz: [Nrg x Naz]  the phase-toheight conversion factor for the interferometric pair
      theta_look: [Nrg x Naz]  off-nadir angle in slant range geometry
      ground_slope: [Nrg x Naz]  terrain slope in slant range geometry

    OUTPUT
      normalization_factor: [Nrg x Naz] power normalization factor"""

    # Phase-to-crossrange conversion factor
    kxr = kz * np.sin(theta_look)

    # Cross-range coordinate of the forest top height in the middle of the
    # resolution cell
    xr0 = eq_height / (np.cos(theta_look)) / (np.tan(theta_look) - np.tan(ground_slope))

    # Spread along the cross-range direction of a surface due to acquisition
    # geometry (pho_r = range resolution)
    delta_xr = pho_r / np.tan(theta_look - ground_slope)

    # Model power (up to the absolute forest reflectivity density)
    P = 2 * (
        (xr0 + delta_xr / 2) * (1 - np.sin(kxr * (xr0 + delta_xr / 2)) / (kxr * (xr0 + delta_xr / 2)))
        + delta_xr / 2 * (1 - np.sin(kxr * delta_xr / 2) / (kxr * delta_xr / 2))
    )

    # Surface normalization
    normalization_factor = P / np.sin(theta_look - ground_slope)

    return normalization_factor


def kzInterp(data_stack_in, kz_stack_in, kz0, pol_name):
    """It generates a synthetic SLC stack of acquisitions by interpolating the original stack of SLC data_stack
    defined over the kz axis specified by kz_stack in correspondence of the desired kz0.

    INPUT
    data_stack: stack of calibrated, ground steered slc images.
                It is a dictionary of two nested dictionaries where:
                 each data_stack[ acquisition_name_string ][ polarization_name_string ] is an array of shape [Nrg x Naz]
     kz_stack:  array of phase-to-height conversion
                factors (kz) OR [num_acq x 1] vector of kz for constant geometries. Needed if num_acq > 2.
                It is a dictionary where kz_stack[acquisition_name_string] is an array of shape [Nrg x Naz]
      kz0:  scalar or [Nrg x Naz] desired phase-to-height

    OUTPUT
      Ikz0: [Nrg x Naz] synthetic SLC
      mask_extrap: [Nrg x Naz] logical mask true if the desired kz0 is out of the available range"""

    num_acq = len(data_stack_in)
    acq_names = list(data_stack_in.keys())
    first_acq_dict = data_stack_in[acq_names[0]]
    Nrg, Naz = first_acq_dict[pol_name].shape

    # convert input dictionaries in lists, in order to have better performances in the computatios that follows
    data_stack = np.zeros((Nrg, Naz, num_acq), dtype=np.complex64)
    kz_stack = np.zeros((Nrg, Naz, num_acq))
    for acq_idx, acq_name in enumerate(acq_names):
        data_stack[:, :, acq_idx] = data_stack_in[acq_name][pol_name]
        kz_stack[:, :, acq_idx] = kz_stack_in[acq_name]
    del data_stack_in, kz_stack_in

    Nr, Nc, N = data_stack.shape

    # Linear interpolation
    pre_kz_ind = np.zeros((Nr, Nc), dtype=np.int8)
    post_kz_ind = np.zeros((Nr, Nc), dtype=np.int8)
    pre_kz_abs_diff = np.zeros((Nr, Nc)) + np.inf
    post_kz_abs_diff = np.zeros((Nr, Nc)) + np.inf

    for n in np.arange(N):
        curr_kz_diff = kz_stack[:, :, n] - kz0
        curr_kz_abs_diff = np.abs(curr_kz_diff)

        pre_kz_mask = curr_kz_diff < 0
        post_kz_mask = np.logical_not(pre_kz_mask)

        # To Be Replaced
        pre_tbr = (np.abs(curr_kz_diff) < pre_kz_abs_diff) & pre_kz_mask
        post_tbr = (np.abs(curr_kz_diff) < post_kz_abs_diff) & post_kz_mask

        pre_kz_ind[pre_tbr] = n
        post_kz_ind[post_tbr] = n

        pre_kz_abs_diff[pre_tbr] = curr_kz_abs_diff[pre_tbr]
        post_kz_abs_diff[post_tbr] = curr_kz_abs_diff[post_tbr]

    # Desired kz_stack out of range (To Be Extrapolated)
    pre_tbe = np.isinf(pre_kz_abs_diff)
    post_tbe = np.isinf(post_kz_abs_diff)

    pre_kz_ind[pre_tbe] = 0
    post_kz_ind[post_tbe] = N - 1

    [C, R] = np.meshgrid(np.arange(Nc), np.arange(Nr))

    kz_pre = kz_stack[R, C, pre_kz_ind]
    kz_post = kz_stack[R, C, post_kz_ind]
    frac_part = (kz0 - kz_pre) / (kz_post - kz_pre)

    Ikz0 = (1 - frac_part) * data_stack[R, C, pre_kz_ind] + frac_part * data_stack[R, C, post_kz_ind]

    mask_extrap = pre_tbe | post_tbe

    Ikz0[mask_extrap] = np.spacing(1)

    return Ikz0, mask_extrap


def is_suitable_z_emph_map(z_emph, Nrg, Naz):
    """It checks whether the input (z_emph) is a suitable 2D map to be used as a
    space-varying height to be emphasized through ground cancellation. It must be:
        - A numpy ndarray instance
        - Two-dimensional
        - Same size (in pixels) as the area currently processed
        - All values greater than zero

    INPUT
    z_emph: variable to be evaluated
    Nrg   : Number of rows currently processed
    Naz   : Number of columns currently processed

    OUTPUT
    is_suitable: """

    is_suitable = False
    if isinstance(z_emph, np.ndarray):
        if np.ndim(z_emph) == 2 and z_emph.shape == (Nrg, Naz) and z_emph.min() > 0:
            is_suitable = True

    return is_suitable
