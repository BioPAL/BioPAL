import numpy as np
from scipy import signal


def coherence2d(data_1, data_2, windows_size=(5, 5)):
    """
    Function to compute the coherence 2D
    :param data_1:
    :param data_2:
    :param windows_size:
    :return:
    """
    interferogram = data_1 * np.conj(data_2)
    data_1_square = data_1 * np.conj(data_1)
    data_2_square = data_2 * np.conj(data_2)

    window = np.ones((np.min((windows_size[0], data_1.shape[0])), np.min((windows_size[1], data_1.shape[1]))))
    window = window / window.sum()

    numerator = signal.convolve2d(interferogram, window)
    numerator = numerator * np.conj(numerator) * np.exp(1j * np.angle(numerator))
    denominator = signal.convolve2d(data_1_square, window) * signal.convolve2d(data_2_square, window)

    coherence = numerator / denominator

    return coherence
