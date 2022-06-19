import pytest

from pyqrusty import py_sum, PyPauli, PySparsePauliOp

# content of test_sample.py
def inc(x):
    return x + 1


def test_answer():
    assert inc(3) == 4

def test_sum():
    assert py_sum(3,4) == 7

def test_paulis():
    p = PyPauli('IXYZ')
    assert p.num_qubits() == 4
    assert p.label() == 'IXYZ'

def test_spop_wrong_lengths():
    p = PyPauli('I')
    q = PyPauli('IXYZ')
    with pytest.raises(Exception):
        spop = PySparsePauliOp([p,q], [1.0 + 0.0j, 1.0 + 0.0j])

def test_spop():
    p = PyPauli('IIII')
    q = PyPauli('IXYZ')
    spop = PySparsePauliOp([p,q], [1.0 + 0.0j, 1.0 + 0.0j])
