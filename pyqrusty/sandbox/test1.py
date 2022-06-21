import math
import time
import numpy as np
import scipy.sparse as scisparse
import scipy.sparse.linalg
import primme

from qiskit.quantum_info.operators.symplectic import *
from qiskit.quantum_info.operators.symplectic.base_pauli import *
from qiskit.quantum_info.operators.symplectic.pauli import *
from fixtures import *

from pyqrusty import *

def timer(msg, f):
    t0 = time.time()
    print("START %s" % (msg,))
    rv = f()
    t1 = time.time()
    print('END ELAPSED %s: %.03fms' % (msg, 1000 * (t1-t0)))
    return rv

def build_spop(fixH, m=0, n=None):
    paulis = [Pauli(label) for label in fixH[0]]
    coeffs = fixH[1]
    print("length of sum: %s" % (len(paulis),))
    if n is None: n = len(paulis)
    return SparsePauliOp(paulis[m:n], coeffs[m:n])

primme_methods = [
   'PRIMME_DEFAULT_METHOD',
   'PRIMME_DYNAMIC',
   'PRIMME_DEFAULT_MIN_TIME',
   'PRIMME_DEFAULT_MIN_MATVECS',
   'PRIMME_Arnoldi',
   'PRIMME_GD',
   'PRIMME_GD_plusK',
   'PRIMME_GD_Olsen_plusK',
   'PRIMME_JD_Olsen_plusK',
   'PRIMME_RQI',
   'PRIMME_JDQR',
   'PRIMME_JDQMR',
   'PRIMME_JDQMR_ETol',
   'PRIMME_STEEPEST_DESCENT',
   'PRIMME_LOBPCG_OrthoBasis',
   'PRIMME_LOBPCG_OrthoBasis_Window'
]

def round_to_zero(v):
    re = v.real
    im = v.imag
    if math.isclose(re, 0.0, abs_tol=1e-7) and math.isclose(im, 0.0, abs_tol=1e-7): return 0.0
    if math.isclose(re, 0.0, abs_tol=1e-7): re = 0.0
    if math.isclose(im, 0.0, abs_tol=1e-7): im = 0.0
    return (re + im * 1j)

def postprocess(evs):
    return [round_to_zero(ev) for ev in evs]

def try_matrix(spmat):
    #dense = spmat.todense()
    #if not scipy.linalg.ishermitian(dense):
    #    print("matrix was not hermitian")
    #    return
    #print("numpy.linalg.eigh: ",  timer("np.linalg.eigh", lambda: np.linalg.eigh(dense)[0][0:1]))
    print("scipy.sparse.linalg.eigsh: ",  timer("scisparse.linalg.eigsh", lambda: scisparse.linalg.eigsh(spmat, k=1)[0]))

    for m in primme_methods:
        x = np.zeros(spmat.shape[0])
        x[np.argmin(spmat.diagonal())] = 1
        v0 = np.array([x])
        print("primme.eigsh(method=%s): %s" %
              (m,postprocess(timer("primme.eigsh", lambda: primme.eigsh(spmat, k=1, tol=1e-6, which='SA', return_eigenvectors=False, method=m)))))
        print("primme.eigsh(method=%s with v0): %s" %
              (m,postprocess(timer("primme.eigsh", lambda: primme.eigsh(spmat, k=1, v0=v0, tol=1e-6, which='SA', return_eigenvectors=False, method=m)))))

spop = build_spop(H8)
#spop = build_spop(H2,1, 2)
#spop = SparsePauliOp([Pauli('YY'), Pauli('XX')], [1.0, 1.0])
#print(spop)

spmat = timer("build sparse matrix", lambda: csr_matrix(spop.to_matrix_binary()))

#print(spmat.diagonal())
#print(spmat.todense())
print(repr(spmat))
try_matrix(spmat)
