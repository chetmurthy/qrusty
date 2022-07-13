// (C) Copyright IBM 2022
//
// This code is licensed under the Apache License, Version 2.0. You may
// obtain a copy of this license in the LICENSE.txt file in the root directory
// of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
//
// Any modifications or derivative works of this code must retain this
// copyright notice, and modified files need to carry a notice indicating
// that they have been altered from the originals.

#![allow(dead_code)]
use std::cmp::* ;
use num_complex::Complex64;
use qrusty::* ;

use qrusty::util::list::* ;
use qrusty::util::BinaryTreeFold ;

pub fn sparse_pauli_op(tc : &fixtures::TestCase, n : usize) -> SparsePauliOp {
    let labels = &tc.labels ;
    let coeffs = &tc.coeffs ;
    assert_eq!(labels.len(), coeffs.len()) ;
    
    let n = min(n, labels.len()) ;
    let ll = &labels[..n] ;
    let cl = &coeffs[..n] ;

    let spop = SparsePauliOp::from_labels(&ll, &cl) ;
    spop.unwrap()
}

fn main() {
    let spop = sparse_pauli_op(&fixtures::H8, 100000) ;
    let _m = spop.to_matrix_accel() ;
}

fn main2() {
    let p = Pauli::new("Y").unwrap() ;

    let sp_mat = p.to_matrix() ;
    
    println!("Y={:?}\n", sp_mat) ;

    let sp_mat_accel = p.to_matrix_accel() ;

    println!("accel Y={:?}\n", sp_mat_accel) ;

    let usv = p.to_unsafe_vectors() ;
    println!("data={:?}\nindices={:?}\nindptr={:?}\n", usv.data, usv.indices, usv.indptr) ;

}
fn main3() {
    let mut bt = BinaryTreeFold::begin(atom(0), cons) ;
    (1..8).for_each(|n| bt.add(atom(n))) ;
    let rv = bt.end() ;
    println!("{}", (*rv).to_string()) ;
}
fn main4() {
    use sprs::{CsMat};
    let eye : CsMat<f64> = CsMat::eye(2);
    let b = &eye + &eye;
    println!("{}",b.to_dense()) ;
}

fn main5() {
    let p = Pauli::new("I").unwrap() ;
    let sp_mat = p.to_matrix() ;
    let b = &sp_mat + &sp_mat ;
    println!("{}",b.to_dense()) ;
}

fn main6() {
    let p = Pauli::new("I").unwrap() ;
    let q = Pauli::new("X").unwrap() ;
    let p_mat = p.to_matrix() ;
    let q_mat = q.to_matrix() ;
    let a = Complex64::new(2.0, 0.0);
    let b = Complex64::new(4.0, 0.0);
    let b  = sprs::binop::csmat_binop(p_mat.view(),  q_mat.view(), |x,y| a * x + b * y) ;
    println!("{}",b.to_dense()) ;
}

fn main7() {
    let spop = SparsePauliOp::from_labels(
        &["I","X"][..],
        &[Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0)][..]).unwrap() ;
    let sp_mat = spop.to_matrix() ;
    println!("{}",sp_mat.to_dense()) ;
}
