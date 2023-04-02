use ark_bls12_381::G1Affine;
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::{PrimeField, prelude::*};
use ark_std::vec::Vec;

use crate::{
    collision_method_pippenger::{quick_msm, MSMRun},
    types::{G1BigInt, G1Projective},
};

// The result of this function is only approximately `ln(a)`
// [`Explanation of usage`]: https://github.com/scipr-lab/zexe/issues/79#issue-556220473
fn ln_without_floats(a: usize) -> usize {
    // log2(a) * ln(2)
    (ark_std::log2(a) * 69 / 100) as usize
}


// WARNING: this function is for test purpose only
// when scalar_bit_length % ws == 0, the sign-bucket-index assert would be
// triggered
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

pub fn multi_scalar_mul<G: AffineCurve>(
    bases: &[G],
    scalars: &[<G::ScalarField as PrimeField>::BigInt],
) -> G::Projective {
    let size = ark_std::cmp::min(bases.len(), scalars.len());
    let scalars = &scalars[..size];
    let bases = &bases[..size];
    let scalars_and_bases_iter = scalars.iter().zip(bases).filter(|(s, _)| !s.is_zero());

    let c = if size < 32 {
        3
    } else {
        ln_without_floats(size) + 2
    };

    let num_bits = <G::ScalarField as PrimeField>::Params::MODULUS_BITS as usize;
    let fr_one = G::ScalarField::one().into_repr();

    let zero = G::Projective::zero();
    let window_starts: Vec<_> = (0..num_bits).step_by(c).collect();

    // Each window is of size `c`.
    // We divide up the bits 0..num_bits into windows of size `c`, and
    // in parallel process each such window.
    let window_sums: Vec<_> = ark_std::cfg_into_iter!(window_starts)
        .map(|w_start| {
            let mut res = zero;
            // We don't need the "zero" bucket, so we only have 2^c - 1 buckets.
            let mut buckets = vec![zero; (1 << c) - 1];
            // This clone is cheap, because the iterator contains just a
            // pointer and an index into the original vectors.
            scalars_and_bases_iter.clone().for_each(|(&scalar, base)| {
                if scalar == fr_one {
                    // We only process unit scalars once in the first window.
                    if w_start == 0 {
                        res.add_assign_mixed(base);
                    }
                } else {
                    let mut scalar = scalar;

                    // We right-shift by w_start, thus getting rid of the
                    // lower bits.
                    scalar.divn(w_start as u32);

                    // We mod the remaining bits by 2^{window size}, thus taking `c` bits.
                    let scalar = scalar.as_ref()[0] % (1 << c);

                    // If the scalar is non-zero, we update the corresponding
                    // bucket.
                    // (Recall that `buckets` doesn't have a zero bucket.)
                    if scalar != 0 {
                        buckets[(scalar - 1) as usize].add_assign_mixed(base);
                    }
                }
            });

            // Compute sum_{i in 0..num_buckets} (sum_{j in i..num_buckets} bucket[j])
            // This is computed below for b buckets, using 2b curve additions.
            //
            // We could first normalize `buckets` and then use mixed-addition
            // here, but that's slower for the kinds of groups we care about
            // (Short Weierstrass curves and Twisted Edwards curves).
            // In the case of Short Weierstrass curves,
            // mixed addition saves ~4 field multiplications per addition.
            // However normalization (with the inversion batched) takes ~6
            // field multiplications per element,
            // hence batch normalization is a slowdown.

            // `running_sum` = sum_{j in i..num_buckets} bucket[j],
            // where we iterate backward from i = num_buckets to 0.
            let mut running_sum = G::Projective::zero();
            buckets.into_iter().rev().for_each(|b| {
                running_sum += &b;
                res += &running_sum;
            });
            res
        })
        .collect();

    // We store the sum for the lowest window.
    let lowest = *window_sums.first().unwrap();

    // We're traversing windows from high to low.
    // sum_i * 2^c
    // ((sum_i * 2^c) + sum_i_1) * 2^c  == sum_i * 2^2c + sum_i_1 * 2^c
    // ...
    // sum_i * 2^(n-1)c  + sum_i_1 * 2^(n-2)c ... sum1 * 2^c + sum0
    lowest
        + &window_sums[1..]
            .iter()
            .rev()
            .fold(zero, |mut total, sum_i| {
                total += sum_i;
                for _ in 0..c {
                    total.double_in_place();
                }
                total
            })
}

pub fn compute_msm_opt(
    point_vec: &Vec<G1Affine>,
    scalar_vec: &Vec<G1BigInt>,
    window_bits: u32,
) -> G1Projective {
    let msm_run = MSMRun {
        points: point_vec.to_vec(),
        scalars: scalar_vec.to_vec(),
        window_bits: window_bits,
        max_batch: 2048,
        max_collisions: 256,
    };
    return quick_msm(&msm_run);
}

pub fn compute_msm_baseline<A>(
    point_vec: &Vec<<A::Projective as ProjectiveCurve>::Affine>,
    scalar_vec: &Vec<<A::ScalarField as PrimeField>::BigInt>,
) -> A::Projective
where
    A: AffineCurve,
{
    multi_scalar_mul(point_vec.as_slice(), scalar_vec.as_slice())
}
