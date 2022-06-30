use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkGroup};
use num_complex::Complex64;

use qrusty::fixtures ;
use qrusty::util ;
use qrusty::Pauli ;
use qrusty::SparsePauliOp ;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("SparsePauliOp");
    group.sample_size(10);
    fixtures::AllTests.iter()
        .for_each(|tc| {
	    let ll = &tc.labels[..] ;
	    let cl = &tc.coeffs[..] ;
	    let spop = SparsePauliOp::from_labels(ll, cl).unwrap() ;
            tc.ops.iter()
                .for_each(|mode| {
                    group.bench_function(format!("{}+{:?}", tc.name, mode), |b| b.iter(|| {
                        spop.to_matrix_mode(&mode) ;
                    }));
                }) ;
        }) ;

    group.finish() ;
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
