"""
Newton method module
--------------------
"""
import numpy as np


def newton_for_geocoding(function_with_jacobian, initial_guess, max_iter=20, tolerance=1.e-10):
    """3D Newton method for geocoding"""
    x = initial_guess

    for _ in range(max_iter):
        f_eval, j_eval = function_with_jacobian(x)
        dx = -_inv_3x3_transpose(j_eval, f_eval)
        x += dx

        stop_condition = np.dot(dx, dx)
        if stop_condition <= tolerance:
            break
    else:
        raise RuntimeError("Newton did not converge: "
                           "maximum number of iterations reached.")

    return x


def _inv_3x3_transpose(mat, f):
    det = + mat[0][0] * (mat[2][2] * mat[1][1] - mat[2][1] * mat[1][2]) \
          - mat[1][0] * (mat[2][2] * mat[0][1] - mat[2][1] * mat[0][2]) \
          + mat[2][0] * (mat[1][2] * mat[0][1] - mat[1][1] * mat[0][2])

    x = f[0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]) + \
        f[1] * (mat[2][0] * mat[1][2] - mat[1][0] * mat[2][2]) + \
        f[2] * (mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1])

    y = f[0] * (mat[2][1] * mat[0][2] - mat[0][1] * mat[2][2]) + \
        f[1] * (mat[0][0] * mat[2][2] - mat[2][0] * mat[0][2]) + \
        f[2] * (mat[2][0] * mat[0][1] - mat[0][0] * mat[2][1])

    z = f[0] * (mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2]) + \
        f[1] * (mat[1][0] * mat[0][2] - mat[0][0] * mat[1][2]) + \
        f[2] * (mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1])

    return np.asarray([x, y, z]) / det
