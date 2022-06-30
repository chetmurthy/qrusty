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

def test_mode():
    with pytest.raises(Exception):
        H2.to_matrix_mode("foo")

def test_H2():
    m1 = csr_matrix(H2.to_matrix())
    assert np.array_equal(m1.todense(), csr_matrix(H2.to_matrix_mode(mode="")).todense())
    m2 = H2.to_matrix_mode(mode="RowwiseUnsafeChunked/100")
    assert np.array_equal(m1.todense(), csr_matrix(m2).todense())

def test_H4():
    m1 = csr_matrix(H4.to_matrix())
    assert np.array_equal(m1.todense(), csr_matrix(H4.to_matrix_mode(mode="")).todense())
    m2 = H4.to_matrix_mode(mode="RowwiseUnsafeChunked/100")
    assert np.array_equal(m1.todense(), csr_matrix(m2).todense())

def test_H6():
    m1 = csr_matrix(H6.to_matrix())
    assert np.array_equal(m1.todense(), csr_matrix(H6.to_matrix_mode(mode="")).todense())
    m2 = H6.to_matrix_mode(mode="RowwiseUnsafeChunked/100")
    assert np.array_equal(m1.todense(), csr_matrix(m2).todense())
