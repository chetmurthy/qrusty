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
from pymisc import *

def timeit(msg, f):
    t0 = time.time()
    print("START %s" % (msg,))
    rv = f()
    t1 = time.time()
    print('END ELAPSED %s: %.03fms' % (msg, 1000 * (t1-t0)))
    return rv

def test_foo():
    assert foo() == 1

def test_axpy():
    x = np.array([1.0 + 2j, 2.0 + 3j, 3.0 + 4j, 4.0 + 5j])
    y = np.ones(4, dtype=complex)
    a = 0.0 + 1.0j
    assert np.allclose(axpy(a,x,y), a * x + y)

def reg(x,tol):
    def f(x):
        if(np.abs(x)<tol): return tol
        return x
    f = np.vectorize(f)
    return f(x)

def precond(spmat, dx,e,x0, tol=1e-14):
        print("tol=%s, dx=%s, e=%s, x0=%s" % (tol, dx,e,x0))
        rv = timeit("precond", lambda: dx/reg(spmat.diagonal()-e,tol))
        print("rv=%s" % (rv,))
        return rv


tol=1e-14
dx=[ 0.00000000e+00+0.j,  0.00000000e+00+0.j,  0.00000000e+00+0.j,
  0.00000000e+00+0.j,  0.00000000e+00+0.j,  0.00000000e+00+0.j,
 -9.57567359e-14+0.j,  0.00000000e+00+0.j,  0.00000000e+00+0.j,
 -9.57428581e-14+0.j,  1.74695127e-01+0.j,  0.00000000e+00+0.j,
  0.00000000e+00+0.j,  0.00000000e+00+0.j,  0.00000000e+00+0.j,
  0.00000000e+00+0.j]
e=-1.7037077186606393
x0=[0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,
 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j]
rv=[ 0.00000000e+00+0.j,  0.00000000e+00+0.j,  0.00000000e+00+0.j,
  0.00000000e+00+0.j,  0.00000000e+00+0.j,  0.00000000e+00+0.j,
 -9.91433685e-14+0.j,  0.00000000e+00+0.j,  0.00000000e+00+0.j,
 -9.91289999e-14+0.j,  8.88073154e-02+0.j,  0.00000000e+00+0.j,
  0.00000000e+00+0.j,  0.00000000e+00+0.j,  0.00000000e+00+0.j,
  0.00000000e+00+0.j]


def test_precond():
    spmat = H2.to_matrix()
    rv2 = precond(spmat, dx, e, x0)
    assert np.allclose(rv, rv2)
