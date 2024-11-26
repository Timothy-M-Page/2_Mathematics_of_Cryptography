import math
import numpy as np

A = np.array([[1, 2], [3, 4]])          # Use np.array to write matrices


def det(matrix):
    det = np.linalg.det(matrix)
    if math.isclose(det, round(det)):   # Avoid floating point precision errors.
        return round(det)
    else:
        return det


def transpose(matrix):
    return np.transpose(matrix)


def matrix_product(matrix1, matrix2):
    try:
        if matrix1.shape[1] != matrix2.shape[0]:    # Ensure matrices are the right shape to be multiplied.
            raise ValueError(f"Matrices of shape {matrix1.shape}, and {matrix2.shape} not compatible.")
        return np.dot(matrix1, matrix2)
    except ValueError as e:
        return(f"{e}")


def rank(matrix):
    result = np.linalg.matrix_rank(matrix)
    return result


def matrix_inverse(matrix):
    try:
        return np.linalg.inv(matrix)
    except np.linalg.LinAlgError:    # The case that the determinant is zero.
        return("The matrix is singular and cannot be inverted.")


def adjugate(matrix):
    rows, cols = matrix.shape
    cofactor_matrix = np.zeros_like(matrix, dtype=float)    # Create an empty matrix of zeros as floats
    for i in range(rows):                                   # the same size as the input matrix.
        for j in range(cols):
            minor_matrix = np.delete(np.delete(matrix, i, axis=0), j, axis=1)  # New entries are determinants
            cofactor_matrix[i, j] = ((-1) ** (i + j)) * det(minor_matrix)      # of the minor matrices
    return np.transpose(cofactor_matrix)


def modular_matrix_inverse(A, mod):
    try:
        A = np.array(A)
        det_A = det(A) % mod
        if det_A == 0:
            return "Matrix determinant is zero mod " + str(mod)     # Determinant must be non zero under the modulus
        det_inv = pow(det_A, -1, mod)
        adj = adjugate(A)
        inv_matrix = (det_inv * adj) % mod      # Modular inverse may be written as modular determinant multipled by adjugate.
        return inv_matrix.astype(int)
    except ValueError:
        return "Determinant " + str(det(A)) + " is not invertible modulo " + str(mod) + "."

