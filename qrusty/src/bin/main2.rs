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
use ndarray::Array;
use ndarray::linalg::Dot;
use std::cmp::min;
use qrusty::* ;

use qrusty::util::list::* ;
use qrusty::util::BinaryTreeFold ;
use rayon::prelude::*;
use rayon_subslice::{ concat_slices, par_concat_slices, unsafe_concat_slices, unsafe_par_concat_slices };

pub fn sparse_pauli_op() -> SparsePauliOp {
    let labels = [
        "IIIIIIIIIIIIIIIIIIZZ",
        "IIIIIIIIIIIIZIIZIIII",
        "IIIIIIIZIZIIIIIIIIII",
        "IZIIZIIIIIIIIIIIIIII",
        "IIIZIZIIIIIIIIIIIIII",
        "IZZIIIIIIIIIIIIIIIII",
        "IIIIZIIZIIIIIIIIIIII",
        "IIIIIIIIIIZIZIIIIIII",
        "IIIIIIIIIIIIIIIZIIZI",
        "IIIIIIIIIIIIIIZIZIII",
        "IIIIIIIIIIIIIIIZIIZI",
        "IIIIIIIIIIIIZZIIIIII",
        "IIIIIIZZIIIIIIIIIIII",
        "IZZIIIIIIIIIIIIIIIII",
        "IIZZIIIIIIIIIIIIIIII",
        "ZZIIIIIIIIIIIIIIIIII",
    ] ;
    let coeffs = vec![Complex64::new(1.0, 0.0); 16] ;

    assert_eq!(labels.len(), coeffs.len()) ;
    
    let spop = SparsePauliOp::from_labels(&labels, &coeffs) ;
    spop.unwrap()
}

fn main() {
    let spop = sparse_pauli_op() ;
    let sp_mat = spop.to_matrix_mode(&Some(AccelMode::RayonChunked(100))) ;
    let rows = sp_mat.rows() ;
    let step = 1024 ;
    let chunks : Vec<(u64, u64)> =
	(0..(rows as u64))
	.into_iter()
	.step_by(step)
            .map(|n| (n,min(n + step as u64, rows as u64)))
            .collect();

    let v : Array::<Complex64, _> = Array::zeros(rows) ;

    let wchunks : Vec<Vec<Complex64>> = chunks.par_iter()
	.map(|(lo,hi)| {
	    let v_rc : Vec<Complex64> =
		(*lo..*hi).map(|rowind| {
		    sp_mat.outer_view(rowind as usize).unwrap().dot(&v)
		})
		.collect();
	    v_rc
	})
    .collect() ;

    let w_slices : Vec<&[Complex64]> = wchunks.iter().map(|v| { &v[..] }).collect()  ;
    let w = Array::from_vec(unsafe_par_concat_slices(&w_slices[..])) ;
}
