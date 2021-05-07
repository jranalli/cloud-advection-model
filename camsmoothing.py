import numpy as np


def cam(input_dt, input_ghi, plant_dx, plant_distr, cloudspeed, ref_pos=0):
    """
    Compute a smoothed time series based on a given plant distribution.

    :param input_dt:
        Time resolution for input measurement in seconds.

    :param input_ghi:
        GHI time series for the input.

    :param plant_dx:
        Spatial resolution of the plant distribution in meters.

    :param plant_distr:
        Plant generation distribution d(x). Pad with trailing zeros to increase
        the frequency resolution of the filter.

    :param cloudspeed:
        Cloud speed in the plant direction in meters per second.

    :param ref_pos:
        The position of the reference within the plant distribution in meters.
        Defaults to zero, corresponding to the leading edge of the plant.

    :return: smoothed, (frequency, filter)
        smoothed - A smoothed version of the input GHI time series.
        frequency - the frequency axis for the filter
        filter - the transfer function representing the plant
    """

    # Compute the FFT and frequency vector of the input time series.
    input_fft = np.fft.fft(input_ghi) * 2 / len(input_ghi)
    input_freq = np.fft.fftfreq(input_ghi.shape[-1], input_dt)

    # Normalize the plant distribution. Then compute the FFT of the plant.
    qint = sum(plant_distr)
    qdist = plant_distr / qint
    plant_fft = np.fft.fft(qdist)  # What does it look like in f domain

    # Compute the effective plant time and frequency vectors
    plant_time = plant_dx / np.abs(cloudspeed)
    plant_freq = np.fft.fftfreq(plant_distr.shape[-1], plant_time)

    # Interpolate so that the frequency axes coincide.
    sortinds = plant_freq.argsort()  # Sort the frequencies for interpolation.
    # Interpolate Magnitude and Phase separately to preserve complex numbers.
    interp_mag = np.interp(input_freq, plant_freq[sortinds], np.abs(plant_fft[sortinds]))
    interp_phase = np.interp(input_freq, plant_freq[sortinds], np.angle(plant_fft[sortinds]))
    interp_plant_fft = interp_mag * np.exp(1j * interp_phase)

    # Shift the phase shift to account for the position of the reference.
    t_delay = ref_pos / cloudspeed
    if cloudspeed > 0:
        interp_plant_fft = interp_plant_fft * np.exp(1j * input_freq * (2 * np.pi) * t_delay)
    else:
        interp_plant_fft = np.conj(interp_plant_fft * np.exp(1j * input_freq * (2 * np.pi) * -t_delay))

    # Compute the filtered frequency response.
    filtresp = input_fft * interp_plant_fft

    # Inverse FFT to yield the modeled time series, only keep the real part.
    invsig = np.fft.ifft(filtresp * len(input_ghi) / 2)
    invsig = np.real(invsig)

    return invsig, (input_freq, interp_plant_fft)
