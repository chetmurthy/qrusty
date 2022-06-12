use criterion::{black_box, criterion_group, criterion_main, Criterion};
use qrusty::Pauli ;

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("Pauli::new", |b| b.iter(|| Pauli::new("I")));
    c.bench_function("Pauli::new+to_matrix", |b| b.iter(|| Pauli::new("I").unwrap().to_matrix()));
    c.bench_function("Pauli::new(IIIIIIIIIIIIIIIIIIYXXY)", |b| b.iter(|| Pauli::new("IIIIIIIIIIIIIIIIIIYXXY")));
    c.bench_function("Pauli::new(IIIIIIIIIIIIIIIIIIYXXY)+to_matrix", |b| b.iter(|| Pauli::new("IIIIIIIIIIIIIIIIIIYXXY").unwrap().to_matrix()));
    c.bench_function("Pauli::new(IIIIIIIIIIIIIIIIIIYXXY)+to_matrix_accel", |b| b.iter(|| Pauli::new("IIIIIIIIIIIIIIIIIIYXXY").unwrap().to_matrix_accel()));
    c.bench_function("Pauli::new(IIIIIIIIIIIIIIIIIIYXXY)+to_triplets", |b| b.iter(|| Pauli::new("IIIIIIIIIIIIIIIIIIYXXY").unwrap().to_triplets()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
