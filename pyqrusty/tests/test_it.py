import pytest

from pyqrusty import py_sum, Pauli, SparsePauliOp

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
