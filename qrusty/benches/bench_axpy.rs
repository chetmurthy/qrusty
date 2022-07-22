// (C) Copyright IBM 2022
//
// This code is licensed under the Apache License, Version 2.0. You may
// obtain a copy of this license in the LICENSE.txt file in the root directory
// of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
//
// Any modifications or derivative works of this code must retain this
// copyright notice, and modified files need to carry a notice indicating
// that they have been altered from the originals.

use std::mem::MaybeUninit;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use qrusty::Pauli ;
use num_complex::Complex64;
use ndarray::{Array, Dim};
use ndarray_rand::rand_distr::Uniform;
use ndarray_rand::RandomExt;
use qrusty::accel ;

fn prepare_axpy(rows : usize) -> (Complex64, Array<Complex64, Dim<[usize; 1]>>, Array<Complex64, Dim<[usize; 1]>>) {
    let x_re : Array<f64, _> = Array::random(rows, Uniform::new(0., 10.));
    let x_im : Array<f64, _> = Array::random(rows, Uniform::new(0., 10.));
    let x  : Array<Complex64, _> = x_re.iter()
	.zip(x_im.iter())
	.map(|(re,im)| Complex64::new(*re, *im))
	.collect() ;
    let y_re : Array<f64, _> = Array::random(rows, Uniform::new(0., 10.));
    let y_im : Array<f64, _> = Array::random(rows, Uniform::new(0., 10.));
    let y  : Array<Complex64, _> = y_re.iter()
	.zip(y_im.iter())
	.map(|(re,im)| Complex64::new(*re, *im))
	.collect() ;
    let a = Complex64::new(1.0, 2.0) ;
    (a, x, y)
}

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("axpy");
    group.sample_size(20);
    for n in (27..28) {
	let (a,x,y) = prepare_axpy(1<<n) ;
	let label = format!("axpy:2^{}", n) ;
        group.bench_function(label, |b| b.iter(|| {
	    let z = a * &x + &y ;
        }));
	let label = format!("accel::axpy:2^{}", n) ;
        group.bench_function(label, |b| b.iter(|| {
	    let z = accel::axpy(a, &x.view(), &y.view()) ;
        }));
	let label = format!("alloc:2^{}", n) ;
        group.bench_function(label, |b| b.iter(|| {
	    let mut z : Array<Complex64, Dim<[usize; 1]>> = Array::zeros(1<<n) ;
        }));
	let label = format!("uninit:2^{}", n) ;
        group.bench_function(label, |b| b.iter(|| {
	    let mut z : Array<MaybeUninit<Complex64>, Dim<[usize; 1]>> = Array::uninit(1<<n) ;
        }));
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
