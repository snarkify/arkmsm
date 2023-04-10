use ark_bls12_381::G1Affine;
use ark_ec::{AffineCurve, ProjectiveCurve, msm::VariableBaseMSM as BaselineVariableBaseMSM};
use ark_ff::{PrimeField};
use ark_std::vec::Vec;
use ark_msm::{
    types::{G1BigInt, G1Projective},
    msm::VariableBaseMSM
};

// The result of this function is only approximately `ln(a)`
// [`Explanation of usage`]: https://github.com/scipr-lab/zexe/issues/79#issue-556220473
fn ln_without_floats(a: usize) -> usize {
    // log2(a) * ln(2)
    (ark_std::log2(a) * 69 / 100) as usize
}

pub fn get_opt_window_size(num_points: usize) -> usize {
    if num_points < 32 {
        3
    } else {
        let mut ws = ln_without_floats(num_points) + 2;
        if num_points >= (1 << 9) && num_points < (1 << 15) {
            ws += 2;
        }
        if num_points >= (1 << 18) {
            ws -= 1;
        }
        ws
    }
}

pub fn compute_msm_opt(
    point: &Vec<G1Affine>,
    scalar: &Vec<G1BigInt>,
    window_bits: u32,
) -> G1Projective {
    VariableBaseMSM::multi_scalar_mul_custom(
        &point.to_vec(),
        &scalar.to_vec(),
        window_bits,
        2048,
        256
    )
}

pub fn compute_msm_baseline<A>(
    point_vec: &Vec<<A::Projective as ProjectiveCurve>::Affine>,
    scalar_vec: &Vec<<A::ScalarField as PrimeField>::BigInt>,
) -> A::Projective
where
    A: AffineCurve,
{
    BaselineVariableBaseMSM::multi_scalar_mul(point_vec.as_slice(), scalar_vec.as_slice())
}
