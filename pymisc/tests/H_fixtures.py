# (C) Copyright IBM 2020, 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

import numpy as np
from pyqrusty import *


def build_spop(labels, coeffs):
    paulis = [Pauli(l) for l in labels]
    return SparsePauliOp(paulis, coeffs)


# ================================================================
H2 = build_spop(
    [
        "IIII",
        "IIIZ",
        "IIYY",
        "IIXX",
        "IIZI",
        "IZII",
        "YYII",
        "XXII",
        "ZIII",
        "IIZZ",
        "IZIZ",
        "IZYY",
        "IZXX",
        "YYIZ",
        "XXIZ",
        "YYYY",
        "XXYY",
        "YYXX",
        "XXXX",
        "ZIIZ",
        "ZIYY",
        "ZIXX",
        "IZZI",
        "YYZI",
        "XXZI",
        "ZIZI",
        "ZZII",
    ],
    coeffs=[
        -0.53063928 + 0.0j,
        0.18388659 + 0.0j,
        -0.00170526 + 0.0j,
        -0.00170526 + 0.0j,
        -0.30789476 + 0.0j,
        0.18388659 + 0.0j,
        -0.00170526 + 0.0j,
        -0.00170526 + 0.0j,
        -0.30789476 + 0.0j,
        0.09918328 + 0.0j,
        0.1459377 + 0.0j,
        0.03292021 + 0.0j,
        0.03292021 + 0.0j,
        0.03292021 + 0.0j,
        0.03292021 + 0.0j,
        0.04367378 + 0.0j,
        0.04367378 + 0.0j,
        0.04367378 + 0.0j,
        0.04367378 + 0.0j,
        0.14285706 + 0.0j,
        0.03462547 + 0.0j,
        0.03462547 + 0.0j,
        0.14285706 + 0.0j,
        0.03462547 + 0.0j,
        0.03462547 + 0.0j,
        0.14863725 + 0.0j,
        0.09918328 + 0.0j,
    ],
)
