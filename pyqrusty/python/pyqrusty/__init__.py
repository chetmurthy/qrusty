# (C) Copyright IBM 2020, 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

from .pyqrusty import Pauli, SparsePauliOp, SpMat, spmat_dot_densevec, \
                      axpby, axpy, ax, precond, precond2, a_spmat_p_b_spmat

def csr_matrix(m):
    (shape, data, indices, indptr) = m.export()
    from scipy.sparse import csr_matrix
    return csr_matrix((data, indices, indptr), shape=shape, dtype=complex)
