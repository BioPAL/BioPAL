import numpy as np

def shift_image(data, shift_x, shift_y):
    x_frequency_axis = np.fft.fftfreq(data.shape[1])[np.newaxis, :]
    y_frequency_axis = np.fft.fftfreq(data.shape[0])[:, np.newaxis]

    data_shifted = np.fft.ifft(np.fft.fft(data, axis=1) * np.exp(-1j * 2 * np.pi * x_frequency_axis * shift_x), axis=1)
    data_shifted = np.fft.ifft(np.fft.fft(data_shifted, axis=0) * np.exp(-1j * 2 * np.pi * y_frequency_axis * shift_y), axis=0)
    return data_shifted.astype(data.dtype)