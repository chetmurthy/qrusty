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
