import numpy as np
import scipy.sparse as sps
import pytest

from pyqrusty import *

# content of test_sample.py
def inc(x):
    return x + 1


def test_answer():
    assert inc(3) == 4

def test_sum():
    assert py_sum(3,4) == 7

def test_paulis():
    p = Pauli('IXYZ')
    assert p.num_qubits() == 4
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
    np.array_equal(I, spmat.todense())

def test_pauli_X():
    p = Pauli("X")
    spmat = csr_matrix(p.to_matrix())
    np.array_equal(X, spmat.todense())

def test_pauli_Z():
    p = Pauli("Z")
    spmat = csr_matrix(p.to_matrix())
    np.array_equal(Z, spmat.todense())

def test_pauli_Y():
    p = Pauli("Y")
    spmat = csr_matrix(p.to_matrix())
    np.array_equal(Y, spmat.todense())

def test_pauli_I_plus_X():
    pI = Pauli("I")
    pX = Pauli("X")
    spop = SparsePauliOp([pI,pX], [1.0 + 0.0j, 2.0 + 0.0j])
    spmat = csr_matrix(spop.to_matrix())
    np.array_equal(I+2.0 * X, spmat.todense())
