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

from pymisc import *

def test_foo():
    assert foo() == 1

def test_axpy():
    x = np.array([1.0 + 2j, 2.0 + 3j, 3.0 + 4j, 4.0 + 5j])
    y = np.ones(4, dtype=complex)
    a = 0.0 + 1.0j
    assert np.allclose(axpy(a,x,y), a * x + y)
