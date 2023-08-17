use ark_bls12_381::G1Affine;
use ark_msm::{msm::VariableBaseMSM, utils::generate_msm_inputs};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

fn bench_window_size(c: &mut Criterion) {
    let mut group = c.benchmark_group("msm");
    group.sample_size(20);
    for size in 10..20 {
        let (point_vec, scalar_vec) = generate_msm_inputs::<G1Affine>(1 << size);
        let point_vec = black_box(point_vec);
        let scalar_vec = black_box(scalar_vec);

        for window_size in (size - 3)..(size + 3) {
            let input = (size, window_size);
            let benchmark_id =
                BenchmarkId::new("ArkMSM", format!("k={}, ws={}", size, window_size));
            group.bench_with_input(benchmark_id, &input, |b, _input| {
                b.iter(|| {
                    let _ = VariableBaseMSM::multi_scalar_mul_custom(
                        &point_vec,
                        &scalar_vec,
                        window_size,
                        2048,
                        256,
                        true,
                    );
                })
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_window_size);
criterion_main!(benches);
