"""
CALIMA - Component Analysis and Line Identification for MAsers
Author: Pablo Santo-Tomás Ros
Affiliation: Instituto de Astrofísica de Andalucía (IAA-CSIC)
Year: 2025
License: GPL v3 (https://img.shields.io/badge/License-GPL%20v3-blue.svg)

Description:
This script detects spectral components in maser emission and fits Gaussian profiles
to determine their spatial coordinates.
"""

import os
import warnings

import numpy as np

try:
    import pandas as pd
except ImportError:
    import pip

    pip.main(["install", "pandas"])
    import pandas as pd

# //////////////////////////////////////////////////////////
#   DEFINITIONS
# //////////////////////////////////////////////////////////


def read_input_file(file_path):
    """Reads a file with key: value pairs and returns a dictionary."""
    if not os.path.isfile(file_path):
        warnings.warn(f"Error: File '{file_path}' not found.")

        return None
    inputs = {}
    try:
        with open(file_path, "r", encoding="utf-8") as file:
            for line in file:
                line = line.strip()
                if line and not line.startswith("#"):
                    key_value = line.split(":", 1)
                    if len(key_value) == 2:
                        key, value = key_value
                        inputs[key.strip()] = value.strip()
                    else:
                        warnings.warn(f"Invalid line format -> {line}")
    except Exception as e:
        warnings.warn(f"Error reading file: {e}")
        return None

    return inputs


def pix_to_Units(Units, pixUnits, coordsys):
    """Converts pixel coordinates to desired units"""
    if freq_axis == 2:
        new_coords = coordsys.convert(
            coordin=[pixUnits[0], pixUnits[1], 1, 1],
            unitsin=["pix", "pix", "Hz", ""],
            unitsout=[Units, Units, "Hz", ""],
        )
    elif freq_axis == 3:
        new_coords = coordsys.convert(
            coordin=[pixUnits[0], pixUnits[1], 1, 1],
            unitsin=["pix", "pix", "", "Hz"],
            unitsout=[Units, Units, "", "Hz"],
        )
    return new_coords


def is_in_ellipse(CentralMaxCoords, Maxcoords, beamMaj, beamMin, beamPA):
    """Checks whether a given point lies within an ellipse of given axes"""
    x = Maxcoords[0]
    y = Maxcoords[1]
    X = CentralMaxCoords[0]
    Y = CentralMaxCoords[1]
    ellipse_equation = (
        x * np.sin(beamPA)
        - y * np.cos(beamPA)
        - X * np.sin(beamPA)
        + Y * np.cos(beamPA)
    ) ** 2 / (beamMaj**2) + (
        x * np.cos(beamPA)
        + y * np.sin(beamPA)
        - X * np.cos(beamPA)
        - Y * np.sin(beamPA)
    ) ** 2 / (
        beamMin**2
    )
    if ellipse_equation <= 1:
        return True
    else:
        return False


def checkMaxima(
    CentralMaxCoords,
    chanRange,
    chanMax,
    coordsys,
    beamMajUnits,
    beamMaj,
    beamMin,
    beamPA,
    factor,
):
    """Checks if the given peaks are surrounded by peaks in adjacent channels that lie within an ellipse of given dimensions"""
    channelRange = range(chanMax - chanRange, chanMax + chanRange + 1)
    channelRange = np.delete(channelRange, chanRange)
    ellipse_cond = np.zeros(2 * chanRange, dtype=bool)
    aux = 0
    ans = False
    for index in range(2 * chanRange):
        try:
            channel = channelRange[index]
            Max = MaximaData[channel]
            if Max > nsigma_cont * RMSData[channel]:
                MaxPixCoords = imstat(filename, box=box_max, chans=f"{channel}")[
                    "maxpos"
                ]
                MaxCoords = pix_to_Units(beamMajUnits, MaxPixCoords, coordsys)
                mod_beamMaj = factor * beamMaj
                mod_beamMin = factor * beamMin
                ellipse_cond[index] = is_in_ellipse(
                    CentralMaxCoords, MaxCoords, mod_beamMaj, mod_beamMin, beamPA
                )
        except IndexError:
            warnings.warn(f"Consecutive channels out of bounds for maxima at {chanMax}")
            index_chans_to_delete = np.linspace(
                index, 2 * chanRange, 2 * chanRange - index, endpoint=False, dtype=int
            )
            ellipse_cond = np.delete(ellipse_cond, index_chans_to_delete)
            break

        if ellipse_cond[index]:
            aux += 1
        else:
            aux = 0

        if aux == chanRange:
            ans = True
            break

    return ans


# //////////////////////////////////////////////////////////
#   MAIN
# //////////////////////////////////////////////////////////

path_to_file = input("Parameter File: ")
inputs = read_input_file(path_to_file)

if inputs is None:
    warnings.warn("Error reading parameter file, terminating task")
else:
    filename = inputs["File Name"]
    box_max = inputs["Source box"]
    box_rms = inputs["Noise box"]

    try:
        outfile = inputs["Output file"]
    except ValueError:
        warnings.warn("Error reading output file name, continue with default name")
        outfile = f"Results_{filename}.csv"

    try:
        nsigma_cont = int(inputs["Nsigma adjacent"])
        nsigma_peak = int(inputs["Nsigma source"])
        nsigma_val = int(inputs["Nsigma valley"])
    except ValueError:
        warnings.warn("Nsigma inputs must be integers, continue with default values")
        nsigma_cont = 2
        nsigma_peak = 3
        nsigma_val = 2

    if nsigma_peak < nsigma_val:
        warnings.warn(
            "Nsigma source must be bigger than Nsigma valley, continue with default values"
        )
        nsigma_cont = 2
        nsigma_peak = 3
        nsigma_val = 2

    try:
        chanRange = int(inputs["N adjacent chans"])
    except ValueError:
        warnings.warn(
            "N adjacent chans must be an integer, continue with default value"
        )
        chanRange = 2

    summary = inputs["Summary Flag"]
    if summary == "False":
        summary = False
    elif summary == "True":
        summary == True
    else:
        warnings.warn(
            "Summary Flag must be bool (True or False), continue with default value"
        )
        summary = False

    verbose = inputs["Verbose Flag"]
    if not verbose == "True":
        if not verbose == "False":
            warnings.warn("Verbose Flag must be bool, continue with default value")
            verbose = "False"

    try:
        axisfactor = float(inputs["Ellipse size"])
    except ValueError:
        warnings.warn("Ellipse size must be float, continue with default value")
        axisfactor = 0.75
        summary = False

    header = imhead(filename)
    refpix = header["refpix"]
    refval = header["refval"]
    increment = header["incr"]
    axisunits = header["axisunits"]
    axisnames = header["axisnames"]
    freq_axis = np.where(axisnames == "Frequency")[0][0]
    myaxes = [0, 1, 2, 3]
    myaxes = np.delete(myaxes, freq_axis)
    MaximaData = imstat(filename, axes=myaxes, box=box_max)["max"]
    RMSData = imstat(filename, axes=myaxes, box=box_rms)["rms"]

    size = MaximaData.size
    DetectionRecord = np.zeros(size, dtype=bool)
    channel = 1
    while channel < size - 1:
        maximum = MaximaData[channel]
        sigma3 = nsigma_peak * RMSData[channel]
        sigma2 = nsigma_val * RMSData[channel]
        Peak = 0.0
        Valley = 0.0
        PeakFound = False
        while maximum > sigma2 and channel < size - 1:
            if maximum > sigma3 and maximum > Peak and maximum - Valley > sigma3:
                Peak = maximum
                channelPeak = channel
                PeakFound = True
            elif PeakFound and Peak - maximum > sigma3:
                Valley = maximum
                PeakFound = False
                DetectionRecord[channelPeak] = True
                Peak = 0.0
            elif not PeakFound and maximum < Valley:
                Valley = maximum

            channel += 1
            maximum = MaximaData[channel]
            sigma3 = nsigma_peak * RMSData[channel]
            sigma2 = nsigma_val * RMSData[channel]

        if PeakFound:
            DetectionRecord[channelPeak] = True

        channel += 1

    PeaksChannels = np.nonzero(DetectionRecord)[0]
    if verbose == "True":
        print(f"The candidates are in channels: {PeaksChannels}")

    NMax = len(PeaksChannels)
    MaxPos = np.zeros((NMax, 2))

    head = imhead(filename, mode="list")
    try:
        beamMajAx = head["beammajor"]["value"]
        beamMajUnits = head["beammajor"]["unit"]
        beamMinAx = head["beamminor"]["value"]
        beamMinUnits = head["beamminor"]["unit"]
        beamPosAngle = head["beampa"]["value"]
        beamPosAngleUnits = head["beampa"]["unit"]
        MaxUnits = ["pix", "pix"]
    except KeyError:
        beamMajAx = head["perplanebeams"]["*0"]["major"]["value"]
        beamMajUnits = head["perplanebeams"]["*0"]["major"]["unit"]
        beamMinAx = head["perplanebeams"]["*0"]["minor"]["value"]
        beamMinUnits = head["perplanebeams"]["*0"]["minor"]["unit"]
        beamPosAngle = head["perplanebeams"]["*0"]["positionangle"]["value"]
        beamPosAngleUnits = head["perplanebeams"]["*0"]["positionangle"]["unit"]

    # We want the position angle to be always in radians
    if beamPosAngleUnits == "deg":
        beamPosAngleUnits = "rad"
        beamPosAngle = beamPosAngle * np.pi / 180

    if not beamMajUnits == beamMinUnits:
        raise NameError("Beam major and minor axis units do not coincide")

    ia.open(filename)
    coord_sys = ia.coordsys()
    ia.close()

    MaxPosRecord = np.zeros((PeaksChannels[-1] + 1, 2))
    MaxPixPos = np.zeros((NMax, 2))
    iteration = 1
    MaxPosPixFitError = np.zeros(2)
    record = open(outfile, "a")

    # Make header
    refvel = coord_sys.frequencytovelocity(
        value=refval[freq_axis],
        frequnit=axisunits[freq_axis],
        doppler="radio",
        velunit="km/s",
    )[0]
    addedvel = coord_sys.frequencytovelocity(
        value=refval[freq_axis] + increment[freq_axis],
        frequnit=axisunits[freq_axis],
        doppler="radio",
        velunit="km/s",
    )[0]
    incrvel = addedvel - refvel
    for index, unit in enumerate(axisunits):
        if unit == "rad":
            refval[index] = refval[index] * 360 / (2 * np.pi)
            increment[index] = increment[index] * 360 / (2 * np.pi)
            axisunits[index] = "deg"

    line1 = f"axisnames,{axisnames[0]},{axisnames[1]},{axisnames[2]},{axisnames[3]},Velocity,,,,,,,,,\n"
    line2 = (
        f"refpix,{refpix[0]},{refpix[1]},{refpix[2]},{refpix[3]},{refpix[3]},,,,,,,,,\n"
    )
    line3 = f"axisunits,{axisunits[0]},{axisunits[1]},{axisunits[2]},{axisunits[3]},km/s,,,,,,,,\n"
    line4 = (
        f"refval,{refval[0]},{refval[1]},{refval[2]},{refval[3]},{refvel},,,,,,,,,\n"
    )
    line5 = f"increment,{increment[0]},{increment[1]},{increment[2]},{increment[3]},{incrvel},,,,,,,,,\n"
    record.write(line1)
    record.write(line2)
    record.write(line3)
    record.write(line4)
    record.write(line5)

    aux = 1
    for index in range(NMax):
        channel = PeaksChannels[index]
        stats = imstat(filename, chans=str(channel), box=box_max)
        MaxPixPos[index, 0] = stats["maxpos"][0]  # in pixels
        MaxPixPos[index, 1] = stats["maxpos"][1]  # in pixels
        MaxPos = pix_to_Units(beamMajUnits, MaxPixPos[index, :], coord_sys)[0:2]
        DetectionRecord[channel] = checkMaxima(
            MaxPos,
            chanRange,
            channel,
            coord_sys,
            beamMajUnits,
            beamMajAx,
            beamMinAx,
            beamPosAngle,
            axisfactor,
        )
        if DetectionRecord[channel]:
            MaxPosRecord[channel, :] = MaxPos
            fit = imfit(filename, box=box_max, chans=str(channel))
            if fit["converged"][0]:
                MaxPosPixFit = fit["results"]["component0"]["pixelcoords"]
                lat_error = fit["deconvolved"]["component0"]["shape"]["direction"][
                    "error"
                ]["latitude"]["value"]
                lat_error_units = fit["deconvolved"]["component0"]["shape"][
                    "direction"
                ]["error"]["latitude"]["unit"]
                long_error = fit["deconvolved"]["component0"]["shape"]["direction"][
                    "error"
                ]["longitude"]["value"]
                long_error_units = fit["deconvolved"]["component0"]["shape"][
                    "direction"
                ]["error"]["longitude"]["unit"]
                pix_arcsec = fit["pixelsperarcsec"]

                if long_error_units == lat_error_units:
                    error_units = long_error_units
                    if error_units == "arcsec":
                        lat_error_pix = lat_error * pix_arcsec[0]
                        long_error_pix = long_error * pix_arcsec[1]
                    elif error_units == "deg":
                        lat_error_pix = lat_error * 3600 * pix_arcsec[0]
                        long_error_pix = long_error * 3600 * pix_arcsec[1]
                    else:
                        raise ValueError("Error units not in arcsec nor deg")
                else:
                    if lat_error_units == "arcsec" and long_error_units == "deg":
                        lat_error_pix = lat_error * pix_arcsec[0]
                        long_error_pix = long_error * 3600 * pix_arcsec[1]
                    elif lat_error_units == "deg" and long_error_units == "arcsec":
                        lat_error_pix = lat_error * 3600 * pix_arcsec[0]
                        long_error_pix = long_error * pix_arcsec[1]
                    else:
                        raise ValueError("Error units not in arcsec nor deg")

                MaxPosPixFitError[0] = long_error_pix
                MaxPosPixFitError[1] = lat_error_pix
                MaxPosFit = pix_to_Units("deg", MaxPosPixFit, coord_sys)[0:2]
                MaxPosFitError = np.zeros(2)
                MaxPosFitError[0] = (
                    pix_to_Units("deg", MaxPosPixFit + MaxPosPixFitError, coord_sys)[0]
                    - MaxPosFit[0]
                )

                MaxPosFitError[1] = lat_error

                frequency = fit["results"]["component0"]["spectrum"]["frequency"]["m0"][
                    "value"
                ]
                frequency_units = fit["results"]["component0"]["spectrum"]["frequency"][
                    "m0"
                ][
                    "unit"
                ]  # We are assuming all frequencies are in the same units as the reference one
                velocity = coord_sys.frequencytovelocity(
                    value=frequency,
                    frequnit=frequency_units,
                    doppler="radio",
                    velunit="km/s",
                )[0]
                peak_flux = fit["results"]["component0"]["peak"]["value"]
                peak_flux_error = fit["results"]["component0"]["peak"]["error"]
                peak_flux_unit = fit["results"]["component0"]["peak"]["unit"]
                flux_density = fit["results"]["component0"]["flux"]["value"][0]
                flux_density_error = fit["results"]["component0"]["flux"]["error"][0]
                flux_density_unit = fit["results"]["component0"]["flux"]["unit"]

                if aux == 1:
                    firstline = f"Channel,RA (pix),DEC (pix),RAerror (pix),DECerror (pix),RA (deg),DEC (deg),RAerror ({lat_error_units}),DECerror ({lat_error_units}),Freq ({frequency_units}),Vel (km/s),PeakFlux ({peak_flux_unit}),PeakFluxError ({peak_flux_unit}),FluxDens ({flux_density_unit}),FluxDensError ({flux_density_unit})\n"
                    record.write(firstline)
                    aux += 1
                newline = f"{channel},{MaxPosPixFit[0]},{MaxPosPixFit[1]},{MaxPosPixFitError[0]},{MaxPosPixFitError[1]},{MaxPosFit[0]},{MaxPosFit[1]},{MaxPosFitError[0]},{MaxPosFitError[1]},{frequency},{velocity},{peak_flux},{peak_flux_error},{flux_density},{flux_density_error}\n"
                record.write(newline)

                if summary:
                    if iteration == 1:
                        imfit(
                            filename,
                            box=box_max,
                            chans=str(channel),
                            summary=f"imfit_{filename}.log",
                        )
                        file1 = open(f"imfit_{filename}.log", "a")
                        iteration += 1
                    else:
                        imfit(
                            filename,
                            box=box_max,
                            chans=str(channel),
                            summary=f"imfit_{filename}_{index}.log",
                        )
                        file2 = open(f"imfit_{filename}_{index}.log", "r")
                        dataline = file2.readlines()[2]
                        file1.write(dataline)
                        file2.close()
                        os.remove(f"imfit_{filename}_{index}.log")
            else:
                newline = f"{channel},NoCon,NoCon,NoCon,NoCon,NoCon,NoCon,NoCon,NoCon,NoCon,NoCon,NoCon,NoCon,NoCon,NoCon\n"
                record.write(newline)

    if summary:
        file1.close()

    record.close()
    record = pd.read_csv(outfile)
    record.to_csv(outfile, index=False)

    if verbose == "True":
        print(
            f"The maxima have been found in channels: {np.nonzero(DetectionRecord)[0]}"
        )
        # print("The location of the maxima is found in the array: MaxPosRecord")
        # print(f"The units of those coordinates are: {beamMajUnits}")
        if summary:
            print(
                f"The parameters of the Gaussian fit on each maximum are stored in imfit_{filename}.log"
            )

        print(
            f"The results from the Gaussian fit on each maximum are stored in {outfile}"
        )
