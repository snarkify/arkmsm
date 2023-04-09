#[path = "utils/bench_helper.rs"]
mod bench_helper;

use ark_bls12_381::G1Affine;
use ark_msm::utils::generate_msm_inputs;
use criterion::{BenchmarkId, black_box, Criterion, criterion_group, criterion_main};
use bench_helper::{compute_msm_baseline, compute_msm_opt, get_opt_window_size};


fn benchmark_with_baseline(c: &mut Criterion) {
    let mut group = c.benchmark_group("msm");
    for size in 8..13 {
        let (point_vec, scalar_vec) = generate_msm_inputs::<G1Affine>(1 << size);
        let point_vec = black_box(point_vec);
        let scalar_vec = black_box(scalar_vec);
        group.bench_with_input(BenchmarkId::new("Baseline", size), &size, |b, _size| {
            b.iter(|| {
                let _ = compute_msm_baseline::<G1Affine>(&point_vec, &scalar_vec);
            })
        });
        let window_size = get_opt_window_size(1 << size) as u32;
        group.bench_with_input(BenchmarkId::new("Ark-MSM", size), &size, |b, _size| {
            b.iter(|| {
                let _ = compute_msm_opt(&point_vec, &scalar_vec, window_size);
            })
        });
    }
    group.finish();
}

criterion_group!(benches, benchmark_with_baseline);
criterion_main!(benches);
