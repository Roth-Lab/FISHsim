import time
import numpy as np
import scipy.signal
import sparse


def test_sparse_convolution():
    row, col, z = 100, 100, 101  # 3d matrix dimensions
    indexes = np.floor(
        np.random.uniform(0, 1, size=(1000, 3)) * np.array([row, col, z])
    ).astype(int)
    data = np.random.uniform(0, 1, size=1000)

    matrix = np.zeros((row, col, z))
    count = 0
    for index in indexes:
        matrix[index[0], index[1], index[2]] = data[count]
        count += 1

    sparse_matrix = sparse.SparseMatrix3D(data, indexes, (row, col, z), False)

    kernel = np.random.uniform(0, 1000, (3, 3, z))

    # Compute convolution
    start = time.time()
    dense_output = scipy.signal.convolve(matrix, kernel, mode="valid")[:, :, 0]
    print(f"dense: {time.time()-start}")
    start = time.time()
    sparse_output = sparse.sparse_convolve3d(sparse_matrix, kernel)
    print(f"sparse: {time.time()-start}")

    diff = dense_output.flatten() - sparse_output.flatten()
    is_same = True
    for i in range(diff.size):
        if diff[i] > 1e-8:
            is_same = False
            break

    assert is_same == True


def test_simple_sparse_convolution():
    row, col, z = 10, 10, 1  # 3d matrix dimensions
    indexes = np.array([[0, 0, 0], [1, 3, 0], [8, 7, 0], [0, 1, 0]])
    data = np.array([1, 1, 1, 2])

    matrix = np.zeros((row, col, z))
    count = 0
    for index in indexes:
        matrix[index[0], index[1], index[2]] = data[count]
        count += 1

    sparse_matrix = sparse.SparseMatrix3D(data, indexes, (row, col, z), False)

    kernel = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]).reshape((3, 3, 1))

    # Compute convolution
    dense_output = scipy.signal.convolve(matrix, kernel, mode="valid")[:, :, 0]
    sparse_output = sparse.sparse_convolve3d(sparse_matrix, kernel)

    print(dense_output.shape)
    print(sparse_output.shape)

    diff = dense_output.flatten() - sparse_output.flatten()
    is_same = True
    for i in range(diff.size):
        if diff[i] > 1e-8:
            is_same = False
            break

    assert is_same == True


if __name__ == "__main__":
    test_simple_sparse_convolution()
    test_sparse_convolution()
