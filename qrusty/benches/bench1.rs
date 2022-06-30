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
