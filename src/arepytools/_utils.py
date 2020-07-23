"""
Utilities
---------
"""

import os
import numbers
from functools import wraps
import numpy as np

_ERROR_CHECK_TYPE = "Wrong {} type: {} != {}"

_ERROR_CHECK_SIZE_OF_NUMPY_ARRAY = "Vector {} has wrong size: {} != {}"
_ERROR_CHECK_DTYPE_OF_NUMPY_ARRAY = "Vector {} has wrong dtype: {} != {}"
_ERROR_CHECK_SHAPE_OF_NUMPY_ARRAY = "Vector {} has wrong shape: {} != {}"
_ERROR_CHECK_NDIM_OF_NUMPY_ARRAY = "Vector {} has wrong ndim: {} != {}"
_ERROR_CHECK_FIRST_AXIS_SIZE_OF_NUMPY_ARRAY = "Vector {} has wrong first axis size: {} != {}"


def check_type(data, dtype, name="", throw_on_error=True) -> bool:
    """
    Checks that data is of type dtype, raises a RunTimeError

    :param data: input variable
    :param dtype: expected type
    :param name: name of the variable (for the error msg) [optional]
    :param throw_on_error: whether to raise a RunTimeError on check failure [optional, default false]
    """
    if not isinstance(data, dtype):
        if throw_on_error:
            raise RuntimeError(_ERROR_CHECK_TYPE.format(name, type(data), dtype))
        return False
    return True


def check_size_of_numpy_array(vector: np.ndarray, size, name="", throw_on_error=True) -> bool:
    """
    Checks the size of vector, raises a RunTimeError if it is not as expected.

    (It does not check if it is a numpy array)

    :param vector: should be a np.ndarray
    :param size: expected size
    :param name: name of the variable (for the error msg) [optional]
    :param throw_on_error: whether to raise a RunTimeError on check failure [optional, default false]
    """
    if size != vector.size:
        if throw_on_error:
            raise RuntimeError(_ERROR_CHECK_SIZE_OF_NUMPY_ARRAY.format(name, vector.size, size))
        return False
    return True


def check_dtype_of_numpy_array(vector: np.ndarray, dtype, name="", throw_on_error=True) -> bool:
    """
    Checks the data type (dtype) of vector, raises a RunTimeError if not as expected.

    (It does not check if it is a numpy array)

    :param vector: should be a np.ndarray
    :param dtype: expected type of the data inside
    :param name: name of the variable (for the error msg) [optional]
    :param throw_on_error: whether to raise a RunTimeError on check failure [optional, default false]
    """
    if dtype != vector.dtype:
        if throw_on_error:
            raise RuntimeError(_ERROR_CHECK_DTYPE_OF_NUMPY_ARRAY.format(name, vector.dtype, dtype))
        return False
    return True


def check_shape_of_numpy_array(vector: np.ndarray, shape, name="", throw_on_error=True) -> bool:
    """
    Checks the data type (dtype) of vector, raises a RunTimeError if not as expected.

    (It does not check if it is a numpy array)

    :param vector: should be a np.ndarray
    :param shape: expected shape of the data inside
    :param name: name of the variable (for the error msg) [optional]
    :param throw_on_error: whether to raise a RunTimeError on check failure [optional, default false]
    """
    if shape != vector.shape or len(shape) != len(vector.shape):
        if throw_on_error:
            raise RuntimeError(_ERROR_CHECK_SHAPE_OF_NUMPY_ARRAY.format(name, vector.shape, shape))
        return False
    return True


def check_ndim_of_numpy_array(vector: np.ndarray, ndim, name="", throw_on_error=True) -> bool:
    """
    Checks the num of dimensions (ndim) of vector, raises a RunTimeError if not as expected.

    (It does not check if it is a numpy array)

    :param vector: should be a np.ndarray
    :param ndim: expected ndim of the data
    :param name: name of the variable (for the error msg) [optional]
    :param throw_on_error: whether to raise a RunTimeError on check failure [optional, default false]
    """
    if ndim != vector.ndim:
        if throw_on_error:
            raise RuntimeError(_ERROR_CHECK_NDIM_OF_NUMPY_ARRAY.format(name, vector.ndim, ndim))
        return False
    return True


def check_first_axis_size_of_numpy_array(vector: np.ndarray, first_axis_size, name="", throw_on_error=True) -> bool:
    """
    Checks that the size of the first axis (first_axis_size) of vector is equal to first_axis_size,
     raises a RunTimeError if not as expected.

     Example: check first_axis_size=3 is ok for (3,) (3,1,1), but not for (1,3,1)

    (It does not check if it is a numpy array)

    :param vector: should be a np.ndarray
    :param first_axis_size: expected first_axis_size of the data
    :param name: name of the variable (for the error msg) [optional]
    :param throw_on_error: whether to raise a RunTimeError on check failure [optional, default false]
    """
    if first_axis_size != vector.shape[0]:
        if throw_on_error:
            raise RuntimeError(
                _ERROR_CHECK_FIRST_AXIS_SIZE_OF_NUMPY_ARRAY.format(name, vector.shape[0], first_axis_size))
        return False
    return True


def input_data_to_numpy_array_with_checks(data, dtype=None, size=None, shape=None, ndim=None, first_axis_size=None,
                                          name="") -> np.ndarray:
    """
    It takes a data and returns data converted to a vector, after some checks

    If data is a np.ndarray it checks dtype, size, shap, ndim (if provided) and returns the data as it is.
    If data is a list it converts it to a np.ndarray and then performs the same checks as stated above
    If data is a single object, the object is converted to a numpy array (following, shape and ndim if given
    and then the checks are performed.
    be (1, ) if ndim is not provided, (1,1,...) of length ndim if ndim was provided.

    :param data: the data
    :param dtype: expected type
    :param size: expected size
    :param ndim: expected ndim
    :param shape: expected shape
    :param first_axis_size: expected value for shape[0]
    :param name: name for error msg

    :return: return a numpy vector with dtype size as shape as expect, or it raises a RunTimeError
    """
    if isinstance(data, list):
        data_array = np.array(data)
    elif isinstance(data, np.ndarray):
        data_array = data
    else:
        if shape is not None:
            data_array = np.full(shape, data)
        else:
            if ndim is not None:
                shape_out = tuple(1 for _ in range(ndim))
            else:
                shape_out = (1,)
            data_array = np.full(shape_out, data)

    check_numpy_array(data_array, dtype=dtype, size=size, shape=shape, ndim=ndim, first_axis_size=first_axis_size,
                      name=name)
    return data_array


def check_numpy_array(data, dtype=None, size=None, shape=None, ndim=None, first_axis_size=None, name="",
                      throw_on_error=True) -> bool:
    """
    Verify that data is a numpy.ndarray and that it satisfies the specified constraints.
    If some constraints is not as expected an exception is raised.
    To disable this behavior and just get a boolean specify throw_on_error=False.
    Specify a name to customize the error messages.

    :param data: the data
    :param dtype: expected type
    :param size: expected size
    :param ndim: expected ndim
    :param first_axis_size: shape[0]
    :param shape: expected shape
    :param name: name for error msg
    :param throw_on_error: whether to raise or not an exception on check failure
    """
    correctness = check_type(data, np.ndarray, name=name, throw_on_error=throw_on_error)
    if dtype is not None:
        correctness = correctness and check_dtype_of_numpy_array(data, dtype, name=name, throw_on_error=throw_on_error)
    if size is not None:
        correctness = correctness and check_size_of_numpy_array(data, size, name=name, throw_on_error=throw_on_error)
    if shape is not None:
        correctness = correctness and check_shape_of_numpy_array(data, shape, name=name, throw_on_error=throw_on_error)
    if ndim is not None:
        correctness = correctness and check_ndim_of_numpy_array(data, ndim, name=name, throw_on_error=throw_on_error)
    if first_axis_size is not None:
        correctness = correctness and check_first_axis_size_of_numpy_array(data, first_axis_size, name=name,
                                                                           throw_on_error=throw_on_error)
    return correctness


def _size_from_shape(shape: tuple) -> int:
    """
    From a tuple reprensenting a shape of a vector it retrives the size.

    For example (5, 1, 8) --> 40

    :param shape: the shape
    :return: the corresponding size
    """
    if not shape:
        return 0

    out = 1
    for k in shape:
        out *= k
    return out


def check_exists(path_to_file, tag=None):
    """
    Check if path exist, if not raise a runtime error
    :param path_to_file: the path to be checked
    :param tag: tag for error msg
    """
    if not os.path.exists(path_to_file):
        raise RuntimeError("{}: {} not found".format("Requested file" if tag is None else tag, path_to_file))


def check_file_exists(*positions):
    """
    Parametrized decorator.
    :param positions: position(s) of the arguments to be checked.
    Example:

    @check_file_exists(2, 4)
    def fun(a, b , path_in, c, path_in2, path_out):
        pass

    when calling fun the path_in and path_in2 will be validated
    the default arguments will be ignored
    """
    for pos in positions:
        if not isinstance(pos, int):
            raise RuntimeError("Decorator parameter should be one or more int")

    def decorator(fun):
        @wraps(fun)
        def decorated_fun(*args_fun):
            for position in positions:
                if position < len(args_fun):
                    check_exists(args_fun[position], "Argument in position {}".format(position))
            return fun(*args_fun)

        return decorated_fun

    return decorator


def check_args_are_numeric(*positions):
    """
    Parametrized decorator to check if the arguments in certain positions are numeric.
    :param positions:  position(s) of the arguments to be checked.

    """
    for pos in positions:
        if not isinstance(pos, int):
            raise RuntimeError("Decorator parameter should be one or more int")

    def decorator(fun):
        @wraps(fun)
        def decorated_fun(*args):

            for position in positions:
                if position < len(args):
                    if not isinstance(args[position], numbers.Number):
                        raise ValueError("Argument {} of {} has to be Numeric".format(position, fun.__name__))
            return fun(*args)

        return decorated_fun

    return decorator
