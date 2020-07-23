import scipy as sp
import numpy as np

def get_abs_preview(data, min_percentile=1, max_percentile=99):

    data_abs = np.abs(data)
    max_ = sp.percentile(data_abs, max_percentile)
    min_ = sp.percentile(data_abs, min_percentile)

    data_abs[data_abs > max_] = max_
    data_abs[data_abs < min_] = min_
    return data_abs