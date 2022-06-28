from .pyqrusty import Pauli, SparsePauliOp, SpMat

def csr_matrix(m):
    (shape, data, indices, indptr) = m.export()
    from scipy.sparse import csr_matrix
    return csr_matrix((data, indices, indptr), shape=shape, dtype=complex)
