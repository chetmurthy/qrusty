import numpy as np
import scipy.sparse as sps
import pytest

from pyqrusty import py_sum, Pauli, SparsePauliOp, SpMat

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
    spmat2 = SpMat.new_unchecked  (dmat.shape, spmat.data, spmat.indices, spmat.indptr)
    assert repr(spmat2) == "<3x2 sparse matrix of type Complex64\n\twith 6 stored elements in Compressed Sparse Row format>"
