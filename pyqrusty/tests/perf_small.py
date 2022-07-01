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

timer("H2+RowwiseUnsafeChunked/100", lambda: H2.to_matrix_mode("RowwiseUnsafeChunked/100"))
timer("H4+RowwiseUnsafeChunked/100", lambda: H4.to_matrix_mode("RowwiseUnsafeChunked/100"))
timer("H6+RowwiseUnsafeChunked/100", lambda: H6.to_matrix_mode("RowwiseUnsafeChunked/100"))
timer("H8+RayonChunked/100", lambda: H8.to_matrix_mode("RayonChunked/100"))
timer("H8+RowwiseUnsafeChunked/100", lambda: H8.to_matrix_mode("RowwiseUnsafeChunked/100"))
