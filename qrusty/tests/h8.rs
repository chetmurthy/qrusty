// (C) Copyright IBM 2022
//
// This code is licensed under the Apache License, Version 2.0. You may
// obtain a copy of this license in the LICENSE.txt file in the root directory
// of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
//
// Any modifications or derivative works of this code must retain this
// copyright notice, and modified files need to carry a notice indicating
// that they have been altered from the originals.

use num_complex::Complex64;
use ndarray::array ;
use ndarray::Array2;
use sprs::*;
use ::approx::{ abs_diff_eq, assert_abs_diff_eq, abs_diff_ne, assert_abs_diff_ne, } ;
use qrusty::* ;

#[test]
fn h8() {
    return ;
    let tc = &qrusty::fixtures::H8 ;
    let ll = &tc.labels[..] ;
    let cl = &tc.coeffs[..] ;
    let spop = SparsePauliOp::from_labels(ll, cl).unwrap() ;
    assert_abs_diff_eq!(spop.to_matrix_rayon(),
                        spop.to_matrix(),
                        epsilon = 1e-7
    ) ;
}
