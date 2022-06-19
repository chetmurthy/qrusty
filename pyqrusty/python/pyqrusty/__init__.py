from .pyqrusty import Pauli, SparsePauliOp, SpMat

def py_sum(a, b):
    """Returns the sum of two numbers (Python only)"""
    return a + b

def csr_matrix(m):
    (shape, data, indices, indptr) = m.export()
    from scipy.sparse import csr_matrix
    return csr_matrix((data, indices, indptr), shape=shape, dtype=complex)
