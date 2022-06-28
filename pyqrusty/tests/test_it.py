import pytest

import numpy as np
import scipy.sparse as sps

import pyqrusty as pq
#from fixtures import *
from fix2 import *

def test_paulis():
    p = pq.Pauli('IXYZ')
    assert (p.num_qubits() == 4)
    assert p.label() == 'IXYZ'
    assert repr(p) == "Pauli('IXYZ')"

def test_spop_wrong_lengths():
    p = pq.Pauli('I')
    q = pq.Pauli('IXYZ')
    with pytest.raises(Exception):
        spop = pq.SparsePauliOp([p,q], [1.0 + 0.0j, 1.0 + 0.0j])

def test_spop():
    p = pq.Pauli('IIII')
    q = pq.Pauli('IXYZ')
    spop = pq.SparsePauliOp([p,q], [1.0 + 0.0j, 1.0 + 0.0j])
    assert repr(spop) == "SparsePauliOp('IIII','IXYZ', [1+0j, 1+0j])"

def test_spmat():
    dmat = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=complex)
    spmat = sps.csr_matrix(dmat)
    spmat2 = pq.SpMat.new_unchecked(dmat.shape, spmat.data, spmat.indices, spmat.indptr)
    assert repr(spmat2) == "<2x2 sparse matrix of type Complex64\n\twith 4 stored elements in Compressed Sparse Row format>"

def test_spmat3():
    dmat = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], dtype=complex)
    spmat = sps.csr_matrix(dmat)
    spmat2 = pq.SpMat.new_unchecked  (spmat.shape, spmat.data, spmat.indices, spmat.indptr)
    assert repr(spmat2) == "<3x2 sparse matrix of type Complex64\n\twith 6 stored elements in Compressed Sparse Row format>"
    spmat3 = pq.csr_matrix(spmat2)

I = np.array([[1.0, 0.0],[0.0, 1.0]], dtype=complex)
Z = np.array([[1.0, 0.0],[0.0, -1.0]], dtype=complex)
X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
Y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)

def test_pauli_I():
    p = pq.Pauli("I")
    spmat = pq.csr_matrix(p.to_matrix())
    assert np.array_equal(I, spmat.todense())

def test_pauli_X():
    p = pq.Pauli("X")
    spmat = pq.csr_matrix(p.to_matrix())
    assert np.array_equal(X, spmat.todense())

def test_pauli_Z():
    p = pq.Pauli("Z")
    spmat = pq.csr_matrix(p.to_matrix())
    assert np.array_equal(Z, spmat.todense())

def test_pauli_Y():
    p = pq.Pauli("Y")
    spmat = pq.csr_matrix(p.to_matrix())
    assert np.array_equal(Y, spmat.todense())

def test_pauli_I_plus_Y():
    pI = pq.Pauli("I")
    pY = pq.Pauli("Y")
    spop = pq.SparsePauliOp([pI,pY], [1.0 + 0.0j, 2.0 + 0.0j])
    spmat = spop.to_matrix()
    target = "/tmp/I_plus_Y.mtx"
    target2 = "/tmp/I_plus_Y-2.mtx"
    spmat.matrixmarket_write(target)
    spmat2 = pq.csr_matrix(spmat)
    from scipy.io import mmread, mmwrite
    mmwrite(target2, spmat2)
    spmat3 = mmread(target)
    assert np.array_equal(I+2.0 * Y, spmat2.todense())
    assert np.array_equal(I+2.0 * Y, spmat3.todense())
    from os import remove
    remove(target)
    remove(target2)

def test_add():
    pI = pq.Pauli("I")
    pX = pq.Pauli("X")
    matI = pI.to_matrix()
    matX = pX.to_matrix()
    matsum = matI + matX
    l = pq.csr_matrix(matsum)
    assert np.array_equal(l.todense(), (I+X))


def test_sub():
    pI = pq.Pauli("I")
    matI = pI.to_matrix()
    matsub = matI - matI
    l = pq.csr_matrix(matsub)
    assert np.array_equal(l.todense(), (I-I))

def test_nz():
    pI = pq.Pauli("I")
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
    assert np.array_equal(pq.csr_matrix(m).todense(), m2.todense())

def write_roundtrip_gzip(path, m):
    gzpath = str(path) + ".gz"
    import copy
    m = copy.copy(m)
    m.matrixmarket_write(gzpath)
    check_call(['gunzip', gzpath])
    m2 = mmread(path)
    assert np.array_equal(pq.csr_matrix(m).todense(), m2.todense())

def write_roundtrip2(path, spmat, symmetry=None):
    mmwrite(str(path), spmat, symmetry=symmetry)
    spmat2 = mmread(str(path))
    assert np.array_equal(spmat2.todense(), spmat.todense())
    m = pq.SpMat.matrixmarket_read(str(path))
    assert np.array_equal(pq.csr_matrix(m).todense(), spmat.todense())

from subprocess import check_call

def write_roundtrip2_gzip(path, spmat, symmetry=None):
    mmwrite(str(path), spmat, symmetry=symmetry)
    spmat2 = mmread(str(path))
    assert np.array_equal(spmat2.todense(), spmat.todense())
    check_call(['gzip', str(path)])
    gzpath = str(path) + ".gz"
    m = pq.SpMat.matrixmarket_read(gzpath)
    assert np.array_equal(pq.csr_matrix(m).todense(), spmat.todense())

def test_write(tmp_path):
    pI = pq.Pauli("I")
    matI = pI.to_matrix()
    write_roundtrip(tmp_path / "I.mtx", matI)
    spmatI = pq.csr_matrix(matI)
    write_roundtrip2(tmp_path / "I2-gz.mtx", spmatI, symmetry="symmetric")
    write_roundtrip2(tmp_path / "I3-gz.mtx", spmatI, symmetry="hermitian")

    pY = pq.Pauli("Y")
    matY = pY.to_matrix()
    write_roundtrip(tmp_path / "Y.mtx", matY)
    write_roundtrip_gzip(tmp_path / "Y-gz.mtx", matY)
    spmatY = pq.csr_matrix(matY)
    write_roundtrip2(tmp_path / "Y2.mtx", spmatY, symmetry="skew-symmetric")
    write_roundtrip2_gzip(tmp_path / "Y2-gz.mtx", spmatY, symmetry="skew-symmetric")
    write_roundtrip2(tmp_path / "Y3.mtx", spmatY, symmetry="hermitian")
    write_roundtrip2_gzip(tmp_path / "Y3-gz.mtx", spmatY, symmetry="hermitian")
