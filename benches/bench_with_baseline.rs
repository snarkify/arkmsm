use criterion::{BenchmarkId, black_box, Criterion, criterion_group, criterion_main};
use ark_bls12_381::G1Affine;
use ark_msm::{
    utils::generate_msm_inputs,
    msm::VariableBaseMSM
};
use ark_ec::msm::VariableBaseMSM as BaselineVariableBaseMSM;


fn benchmark_with_baseline(c: &mut Criterion) {
    let mut group = c.benchmark_group("msm");
    for size in 8..13 {
        let (point_vec, scalar_vec) = generate_msm_inputs::<G1Affine>(1 << size);
        let point_vec = black_box(point_vec);
        let scalar_vec = black_box(scalar_vec);
        group.bench_with_input(BenchmarkId::new("Baseline", size), &size, |b, _size| {
            b.iter(|| {
                let _ = BaselineVariableBaseMSM::multi_scalar_mul(&point_vec, &scalar_vec);
            })
        });
        group.bench_with_input(BenchmarkId::new("ArkMSM", size), &size, |b, _size| {
            b.iter(|| {
                let _ = VariableBaseMSM::multi_scalar_mul(&point_vec, &scalar_vec);
            })
        });
    }
    group.finish();
}

criterion_group!(benches, benchmark_with_baseline);
criterion_main!(benches);
