# (C) Copyright IBM 2020, 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

import pytest

import time
import numpy as np
import scipy.sparse as sps

from pyqrusty import *
from H_fixtures import *

def timer(msg, f):
    t0 = time.time()
    print("START %s" % (msg,))
    rv = f()
    t1 = time.time()
    print('END ELAPSED %s: %.03fms' % (msg, 1000 * (t1-t0)))
    return rv

def test_paulis():
    p = Pauli('IXYZ')
    assert (p.num_qubits() == 4)
    assert p.label() == 'IXYZ'
    assert repr(p) == "Pauli('IXYZ')"

def test_spop_wrong_lengths():
    p = Pauli('I')
    q = Pauli('IXYZ')
    with pytest.raises(Exception):
        spop = SparsePauliOp([p,q], [1.0 + 0.0j, 1.0 + 0.0j])

def test_spop():
    p = Pauli('IIII')
    q = Pauli('IXYZ')
    spop = SparsePauliOp([p,q], [1.0 + 0.0j, 1.0 + 0.0j])
    assert repr(spop) == "SparsePauliOp('IIII','IXYZ', [1+0j, 1+0j])"

def test_spmat():
    dmat = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=complex)
    spmat = sps.csr_matrix(dmat)
    spmat2 = SpMat.new_unchecked(dmat.shape, spmat.data, spmat.indices, spmat.indptr)
    assert repr(spmat2) == "<2x2 sparse matrix of type Complex64\n\twith 4 stored elements in Compressed Sparse Row format>"
    assert np.array_equal(spmat2.diagonal(), np.diagonal(dmat))

def test_spmat3():
    dmat = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], dtype=complex)
    spmat = sps.csr_matrix(dmat)
    spmat2 = SpMat.new_unchecked  (spmat.shape, spmat.data, spmat.indices, spmat.indptr)
    assert repr(spmat2) == "<3x2 sparse matrix of type Complex64\n\twith 6 stored elements in Compressed Sparse Row format>"
    spmat3 = csr_matrix(spmat2)

I = np.array([[1.0, 0.0],[0.0, 1.0]], dtype=complex)
Z = np.array([[1.0, 0.0],[0.0, -1.0]], dtype=complex)
X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
Y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)

def test_pauli_I():
    p = Pauli("I")
    spmat = csr_matrix(p.to_matrix())
    assert np.array_equal(I, spmat.todense())

def test_pauli_X():
    p = Pauli("X")
    spmat = csr_matrix(p.to_matrix())
    assert np.array_equal(X, spmat.todense())

def test_pauli_Z():
    p = Pauli("Z")
    spmat = csr_matrix(p.to_matrix())
    assert np.array_equal(Z, spmat.todense())

def test_pauli_Y():
    p = Pauli("Y")
    spmat = csr_matrix(p.to_matrix())
    assert np.array_equal(Y, spmat.todense())

def test_pauli_I_plus_Y():
    pI = Pauli("I")
    pY = Pauli("Y")
    spop = SparsePauliOp([pI,pY], [1.0 + 0.0j, 2.0 + 0.0j])
    spmat = spop.to_matrix()
    target = "/tmp/I_plus_Y.mtx"
    target2 = "/tmp/I_plus_Y-2.mtx"
    spmat.matrixmarket_write(target)
    spmat2 = csr_matrix(spmat)
    from scipy.io import mmread, mmwrite
    mmwrite(target2, spmat2)
    spmat3 = mmread(target)
    assert np.array_equal(I+2.0 * Y, spmat2.todense())
    assert np.array_equal(I+2.0 * Y, spmat3.todense())
    from os import remove
    remove(target)
    remove(target2)

def test_add():
    pI = Pauli("I")
    pX = Pauli("X")
    matI = pI.to_matrix()
    matX = pX.to_matrix()
    matsum = matI + matX
    l = csr_matrix(matsum)
    assert np.array_equal(l.todense(), (I+X))


def test_sub():
    pI = Pauli("I")
    matI = pI.to_matrix()
    matsub = matI - matI
    l = csr_matrix(matsub)
    assert np.array_equal(l.todense(), (I-I))

def test_nz():
    pI = Pauli("I")
    matI = pI.to_matrix()
    matI.scale(1e-8 + 0j)
    assert matI.count_zeros() == 2
    mat2 = matI.eliminate_zeros()
    assert mat2.count_zeros() == 0
    assert mat2.nnz() == 0
    

def test_copy():
    pI = Pauli("I")
    matI = pI.to_matrix()
    import copy
    mat2 = copy.copy(matI)

from scipy.io import mmread, mmwrite

def write_roundtrip(path, m):
    import copy
    m = copy.copy(m)
    m.matrixmarket_write(str(path))
    m2 = mmread(path)
    assert np.array_equal(csr_matrix(m).todense(), m2.todense())

def write_roundtrip_gzip(path, m):
    gzpath = str(path) + ".gz"
    import copy
    m = copy.copy(m)
    m.matrixmarket_write(gzpath)
    check_call(['gunzip', gzpath])
    m2 = mmread(path)
    assert np.array_equal(csr_matrix(m).todense(), m2.todense())

def write_roundtrip2(path, spmat, symmetry=None):
    mmwrite(str(path), spmat, symmetry=symmetry)
    spmat2 = mmread(str(path))
    assert np.array_equal(spmat2.todense(), spmat.todense())
    m = SpMat.matrixmarket_read(str(path))
    assert np.array_equal(csr_matrix(m).todense(), spmat.todense())

from subprocess import check_call

def write_roundtrip2_gzip(path, spmat, symmetry=None):
    mmwrite(str(path), spmat, symmetry=symmetry)
    spmat2 = mmread(str(path))
    assert np.array_equal(spmat2.todense(), spmat.todense())
    check_call(['gzip', str(path)])
    gzpath = str(path) + ".gz"
    m = SpMat.matrixmarket_read(gzpath)
    assert np.array_equal(csr_matrix(m).todense(), spmat.todense())

def test_write(tmp_path):
    pI = Pauli("I")
    matI = pI.to_matrix()
    write_roundtrip(tmp_path / "I.mtx", matI)
    spmatI = csr_matrix(matI)
    write_roundtrip2(tmp_path / "I2-gz.mtx", spmatI, symmetry="symmetric")
    write_roundtrip2(tmp_path / "I3-gz.mtx", spmatI, symmetry="hermitian")

    pY = Pauli("Y")
    matY = pY.to_matrix()
    write_roundtrip(tmp_path / "Y.mtx", matY)
    write_roundtrip_gzip(tmp_path / "Y-gz.mtx", matY)
    spmatY = csr_matrix(matY)
    write_roundtrip2(tmp_path / "Y2.mtx", spmatY, symmetry="skew-symmetric")
    write_roundtrip2_gzip(tmp_path / "Y2-gz.mtx", spmatY, symmetry="skew-symmetric")
    write_roundtrip2(tmp_path / "Y3.mtx", spmatY, symmetry="hermitian")
    write_roundtrip2_gzip(tmp_path / "Y3-gz.mtx", spmatY, symmetry="hermitian")

def test_slice():
    p = H2[3]
    assert repr(p) == "(Pauli('IIXX'), (-0.00170526+0j))"
    s = H2[3:8:2]
    print(s)
    assert repr(s) == "[(Pauli('IIXX'), (-0.00170526+0j)), (Pauli('IZII'), (0.18388659+0j)), (Pauli('XXII'), (-0.00170526+0j))]"

def test_spmat_dot_densevec():
    pIX = Pauli("IX")
    matIX = pIX.to_matrix()
    csr_matIX = csr_matrix(pIX.to_matrix())
    x = np.array([1.0 + 2j, 2.0 + 3j, 3.0 + 4j, 4.0 + 5j])
    assert np.allclose(spmat_dot_densevec(matIX, x), csr_matIX.dot(x))

def test_spmat_dot_densevec_H6():
    m = H6.to_matrix()
    csrm = csr_matrix(H6.to_matrix())
    rows = csrm.shape[0]
    x = np.random.random(rows) + np.random.random(rows) * 1j
    assert np.allclose(spmat_dot_densevec(m, x), csrm.dot(x))

def test_axpy():
    x = np.array([1.0 + 2j, 2.0 + 3j, 3.0 + 4j, 4.0 + 5j])
    y = np.ones(4, dtype=complex)
    a = 0.0 + 1.0j
    assert np.allclose(axpy(a,x,y), a * x + y)

def test_axpy2():
    rows = 1<<20
    x = np.random.random(rows) + np.random.random(rows) * 1j
    y = np.random.random(rows) + np.random.random(rows) * 1j
    a = 0.0 + 1.0j
    assert np.allclose(axpy(a,x,y), a * x + y)

def test_ax():
    x = np.array([1.0 + 2j, 2.0 + 3j, 3.0 + 4j, 4.0 + 5j])
    a = 0.0 + 1.0j
    assert np.allclose(ax(a,x), a * x)

def test_ax2():
    rows = 1<<20
    x = np.random.random(rows) + np.random.random(rows) * 1j
    a = 0.0 + 1.0j
    assert np.allclose(ax(a,x), a * x)
