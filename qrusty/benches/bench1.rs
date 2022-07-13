// (C) Copyright IBM 2022
//
// This code is licensed under the Apache License, Version 2.0. You may
// obtain a copy of this license in the LICENSE.txt file in the root directory
// of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
//
// Any modifications or derivative works of this code must retain this
// copyright notice, and modified files need to carry a notice indicating
// that they have been altered from the originals.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use qrusty::Pauli ;

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("Pauli::new", |b| b.iter(|| Pauli::new(&String::from("I"))));
    c.bench_function("Pauli::new+to_matrix", |b| b.iter(|| Pauli::new(&String::from("I")).unwrap().to_matrix()));
    c.bench_function("Pauli::new(IIIIIIIIIIIIIIIIIIYXXY)", |b| b.iter(|| Pauli::new(&String::from("IIIIIIIIIIIIIIIIIIYXXY"))));
    c.bench_function("Pauli::new(IIIIIIIIIIIIIIIIIIYXXY)+to_matrix", |b| b.iter(|| Pauli::new(&String::from("IIIIIIIIIIIIIIIIIIYXXY")).unwrap().to_matrix()));
    c.bench_function("Pauli::new(IIIIIIIIIIIIIIIIIIYXXY)+to_matrix_accel", |b| b.iter(|| Pauli::new(&String::from("IIIIIIIIIIIIIIIIIIYXXY")).unwrap().to_matrix_accel()));
    c.bench_function("Pauli::new(IIIIIIIIIIIIIIIIIIYXXY)+to_unsafe_vectors", |b| b.iter(|| Pauli::new(&String::from("IIIIIIIIIIIIIIIIIIYXXY")).unwrap().to_unsafe_vectors()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
