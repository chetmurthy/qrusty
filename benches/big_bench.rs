use criterion::{black_box, criterion_group, criterion_main, Criterion};
use qrusty::Pauli ;

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("Pauli::new", |b| b.iter(|| Pauli::new(&String::from("I"))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
