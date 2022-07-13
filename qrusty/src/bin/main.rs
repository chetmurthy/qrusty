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
use num_complex::Complex64;
use qrusty::* ;

use qrusty::util::list::* ;
use qrusty::util::BinaryTreeFold ;

pub fn sparse_pauli_op(n : usize) -> SparsePauliOp {
    let labels = vec![
        "IIIIIIIIIIIIIIIIIIYXXY", "IIIIIIIIIIIIIIIIIIYYYY",
        "IIIIIIIIIIIIIIIIIIXXYY", "IIIIIIIIIIIIIIIIIIYYXX",
        "IIIIIIIIIIIIIIIIIIXXXX", "IIIIIIIIIIIIIIIIIIXYYX",
        "IIIIIIIIIIIIIIIIIYZXXY", "IIIIIIIIIIIIIIIIIYZYYY",
        "IIIIIIIIIIIIIIIIIXZXYY", "IIIIIIIIIIIIIIIIIYZYXX",
        "IIIIIIIIIIIIIIIIIXZXXX", "IIIIIIIIIIIIIIIIIXZYYX",
        "IIIIIIIIIIIIIIIIYZZXXY", "IIIIIIIIIIIIIIIIYZZYYY",
        "IIIIIIIIIIIIIIIIXZZXYY", "IIIIIIIIIIIIIIIIYZZYXX",
        "IIIIIIIIIIIIIIIIXZZXXX", "IIIIIIIIIIIIIIIIXZZYYX",
        "IIIIIIIIIIIIIIIYZZZXXY", "IIIIIIIIIIIIIIIYZZZYYY",
        "IIIIIIIIIIIIIIIXZZZXYY", "IIIIIIIIIIIIIIIYZZZYXX",
        "IIIIIIIIIIIIIIIXZZZXXX", "IIIIIIIIIIIIIIIXZZZYYX",
        "IIIIIIIIIIIIIIYZZZZXXY", "IIIIIIIIIIIIIIYZZZZYYY",
        "IIIIIIIIIIIIIIXZZZZXYY", "IIIIIIIIIIIIIIYZZZZYXX",
        "IIIIIIIIIIIIIIXZZZZXXX", "IIIIIIIIIIIIIIXZZZZYYX",
        "IIIIIIIIIIIIIYZZZZZXXY", "IIIIIIIIIIIIIYZZZZZYYY",
        "IIIIIIIIIIIIIXZZZZZXYY", "IIIIIIIIIIIIIYZZZZZYXX",
        "IIIIIIIIIIIIIXZZZZZXXX", "IIIIIIIIIIIIIXZZZZZYYX",
        "IIIIIIIIIIIIYZZZZZZXXY", "IIIIIIIIIIIIYZZZZZZYYY",
        "IIIIIIIIIIIIXZZZZZZXYY", "IIIIIIIIIIIIYZZZZZZYXX",
        "IIIIIIIIIIIIXZZZZZZXXX", "IIIIIIIIIIIIXZZZZZZYYX",
        "IIIIIIIIIIIYZZZZZZZXXY", "IIIIIIIIIIIYZZZZZZZYYY",
        "IIIIIIIIIIIXZZZZZZZXYY", "IIIIIIIIIIIYZZZZZZZYXX",
        "IIIIIIIIIIIXZZZZZZZXXX", "IIIIIIIIIIIXZZZZZZZYYX",
        "IIIIIIIIIIIIIIIIIYXIXY", "IIIIIIIIIIIIIIIIIYYIYY",
        "IIIIIIIIIIIIIIIIIXXIYY", "IIIIIIIIIIIIIIIIIYYIXX",
        "IIIIIIIIIIIIIIIIIXXIXX", "IIIIIIIIIIIIIIIIIXYIYX",
        "IIIIIIIIIIIIIIIIYZXIXY", "IIIIIIIIIIIIIIIIYZYIYY",
        "IIIIIIIIIIIIIIIIXZXIYY", "IIIIIIIIIIIIIIIIYZYIXX",
        "IIIIIIIIIIIIIIIIXZXIXX", "IIIIIIIIIIIIIIIIXZYIYX",
        "IIIIIIIIIIIIIIIYZZXIXY", "IIIIIIIIIIIIIIIYZZYIYY",
        "IIIIIIIIIIIIIIIXZZXIYY", "IIIIIIIIIIIIIIIYZZYIXX",
        "IIIIIIIIIIIIIIIXZZXIXX", "IIIIIIIIIIIIIIIXZZYIYX",
        "IIIIIIIIIIIIIIYZZZXIXY", "IIIIIIIIIIIIIIYZZZYIYY",
        "IIIIIIIIIIIIIIXZZZXIYY", "IIIIIIIIIIIIIIYZZZYIXX",
        "IIIIIIIIIIIIIIXZZZXIXX", "IIIIIIIIIIIIIIXZZZYIYX",
        "IIIIIIIIIIIIIYZZZZXIXY", "IIIIIIIIIIIIIYZZZZYIYY",
        "IIIIIIIIIIIIIXZZZZXIYY", "IIIIIIIIIIIIIYZZZZYIXX",
        "IIIIIIIIIIIIIXZZZZXIXX", "IIIIIIIIIIIIIXZZZZYIYX",
        "IIIIIIIIIIIIYZZZZZXIXY", "IIIIIIIIIIIIYZZZZZYIYY",
        "IIIIIIIIIIIIXZZZZZXIYY", "IIIIIIIIIIIIYZZZZZYIXX",
        "IIIIIIIIIIIIXZZZZZXIXX", "IIIIIIIIIIIIXZZZZZYIYX",
        "IIIIIIIIIIIYZZZZZZXIXY", "IIIIIIIIIIIYZZZZZZYIYY",
        "IIIIIIIIIIIXZZZZZZXIYY", "IIIIIIIIIIIYZZZZZZYIXX",
        "IIIIIIIIIIIXZZZZZZXIXX", "IIIIIIIIIIIXZZZZZZYIYX",
        "IIIIIIIIIIIIIIIIYXIIXY", "IIIIIIIIIIIIIIIIYYIIYY",
        "IIIIIIIIIIIIIIIIXXIIYY", "IIIIIIIIIIIIIIIIYYIIXX",
        "IIIIIIIIIIIIIIIIXXIIXX", "IIIIIIIIIIIIIIIIXYIIYX",
        "IIIIIIIIIIIIIIIYZXIIXY", "IIIIIIIIIIIIIIIYZYIIYY",
        "IIIIIIIIIIIIIIIXZXIIYY", "IIIIIIIIIIIIIIIYZYIIXX"] ;

    let coeffs = [
        "-2.38476799e-06+0.j", "-2.54069063e-06+0.j",
        "-1.55922634e-07+0.j", "-1.55922634e-07+0.j",
        "-2.54069063e-06+0.j", "-2.38476799e-06+0.j",
        "-3.25786104e-06+0.j", "-7.12962163e-06+0.j",
        "-3.87176059e-06+0.j", "-3.87176059e-06+0.j",
        "-7.12962163e-06+0.j", "-3.25786104e-06+0.j",
        "-1.34019018e-04+0.j", "-1.74138457e-04+0.j",
        "-4.01194385e-05+0.j", "-4.01194385e-05+0.j",
        "-1.74138457e-04+0.j", "-1.34019018e-04+0.j",
        "4.94958014e-05+0.j", "6.41626617e-05+0.j",
        "1.46668603e-05+0.j", "1.46668603e-05+0.j",
        "6.41626617e-05+0.j", "4.94958014e-05+0.j",
        "8.55602904e-05+0.j", "9.18732766e-05+0.j",
        "6.31298618e-06+0.j", "6.31298618e-06+0.j",
        "9.18732766e-05+0.j", "8.55602904e-05+0.j",
        "-7.31341568e-03+0.j", "-7.63751298e-03+0.j",
        "-3.24097301e-04+0.j", "-3.24097301e-04+0.j",
        "-7.63751298e-03+0.j", "-7.31341568e-03+0.j",
        "-5.83847754e-05+0.j", "-6.64310595e-05+0.j",
        "-8.04628410e-06+0.j", "-8.04628410e-06+0.j",
        "-6.64310595e-05+0.j", "-5.83847754e-05+0.j",
        "9.48252613e-05+0.j", "1.17580843e-04+0.j",
        "2.27555813e-05+0.j", "2.27555813e-05+0.j",
        "1.17580843e-04+0.j", "9.48252613e-05+0.j",
        "1.93016760e-05+0.j", "2.27560060e-05+0.j",
        "3.45432999e-06+0.j", "3.45432999e-06+0.j",
        "2.27560060e-05+0.j", "1.93016760e-05+0.j",
        "-2.27329986e-06+0.j", "5.24808543e-06+0.j",
        "7.52138529e-06+0.j", "7.52138529e-06+0.j",
        "5.24808543e-06+0.j", "-2.27329986e-06+0.j",
        "-1.43603748e-05+0.j", "-2.37960262e-05+0.j",
        "-9.43565140e-06+0.j", "-9.43565140e-06+0.j",
        "-2.37960262e-05+0.j", "-1.43603748e-05+0.j",
        "2.89120785e-05+0.j", "3.74584947e-05+0.j",
        "8.54641621e-06+0.j", "8.54641621e-06+0.j",
        "3.74584947e-05+0.j", "2.89120785e-05+0.j",
        "-1.64124583e-03+0.j", "-1.83869099e-03+0.j",
        "-1.97445158e-04+0.j", "-1.97445158e-04+0.j",
        "-1.83869099e-03+0.j", "-1.64124583e-03+0.j",
        "-6.22799134e-06+0.j", "3.67293903e-06+0.j",
        "9.90093037e-06+0.j", "9.90093037e-06+0.j",
        "3.67293903e-06+0.j", "-6.22799134e-06+0.j",
         "1.10777376e-05+0.j", "8.53136587e-06+0.j",
        "-2.54637170e-06+0.j", "-2.54637170e-06+0.j",
        "8.53136587e-06+0.j", "1.10777376e-05+0.j",
        "-1.10711295e-05+0.j", "6.60760260e-05+0.j",
        "7.71471555e-05+0.j", "7.71471555e-05+0.j",
        "6.60760260e-05+0.j", "-1.10711295e-05+0.j",
        "-9.61269282e-06+0.j", "-3.07028288e-05+0.j",
        "-2.10901360e-05+0.j", "-2.10901360e-05+0.j"
    ] ;
    assert_eq!(labels.len(), coeffs.len()) ;
    
    let ll = &labels[..n] ;
    let cl: Result<Vec<Complex64>, _> = coeffs[..n]
        .iter()
        .map(|s| util::complex64_from_string(s))
        .collect::<Result<Vec<Complex64>, _>>();
    let cl = cl.unwrap() ;

    let spop = SparsePauliOp::from_labels(&ll, &cl) ;
    spop.unwrap()
}

fn main() {
        let spop = sparse_pauli_op(100) ;
        spop.to_matrix_accel() ;
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
