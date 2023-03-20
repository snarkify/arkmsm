use crate::{
    bucket_msm::BucketMSM,
    types::{G1BigInt, G1Projective, G1_SCALAR_SIZE},
};
use ark_bls12_381::G1Affine;
use ark_ff::BigInteger;

pub struct MSMRun {
    pub points: Vec<G1Affine>,
    pub scalars: Vec<G1BigInt>,
    pub window_bits: u32,
    pub max_batch: u32,
    pub max_collisions: u32,
}

fn msm_slice(scalar: G1BigInt, slices: &mut Vec<u32>, window_bits: u32) {
    let mut temp = scalar;
    for i in 0..slices.len() {
        slices[i] = (temp.as_ref()[0] % (1 << window_bits)) as u32;
        temp.divn(window_bits);
    }

    //    for i in 0..8 {
    //        if sliced[i] > 0x4000 {
    //            sliced[i + 1] += 1;
    //            sliced[i] = if sliced[i] == 0x8000 {
    //                0
    //            } else {
    //                sliced[i] ^ 0x80007FFF
    //            };
    //        }
    //    }
}

pub fn quick_msm(run: &MSMRun) -> G1Projective {
    let mut bucket_msm: BucketMSM<G1Affine> = BucketMSM::new(
        G1_SCALAR_SIZE,
        run.window_bits,
        run.max_batch,
        run.max_collisions,
    );
    let num_slices: u32 = (G1_SCALAR_SIZE + run.window_bits - 1) / run.window_bits;
    let mut slices: Vec<u32> = vec![0; num_slices as usize];

    let scalars_and_bases_iter = run
        .scalars
        .iter()
        .zip(&run.points)
        .filter(|(s, _)| !s.is_zero());
    scalars_and_bases_iter.for_each(|(&scalar, point)| {
        msm_slice(scalar, &mut slices, run.window_bits);
        bucket_msm.process_point_and_slices(&point, &slices);
    });

    bucket_msm.process_complete();
    return bucket_msm.msm_reduce();
}

#[cfg(test)]
mod collision_method_pippenger_tests {
    use super::*;

    #[test]
    fn test_msm_slice_window_size_1() {
        let scalar = G1BigInt::from(0b101);
        let mut slices: Vec<u32> = vec![0; 3];
        msm_slice(scalar, &mut slices, 1);
        assert_eq!(slices.iter().eq([1, 0, 1].iter()), true);
    }

    #[test]
    fn test_msm_slice_window_size_2() {
        let scalar = G1BigInt::from(0b000110);
        let mut slices: Vec<u32> = vec![0; 3];
        msm_slice(scalar, &mut slices, 2);
        assert_eq!(slices.iter().eq([2, 1, 0].iter()), true);
    }

    #[test]
    fn test_msm_slice_window_size_3() {
        let scalar = G1BigInt::from(0b111010000);
        let mut slices: Vec<u32> = vec![0; 3];
        msm_slice(scalar, &mut slices, 3);
        assert_eq!(slices.iter().eq([0, 2, 7].iter()), true);
    }

    #[test]
    fn test_msm_slice_window_size_16() {
        let scalar = G1BigInt::from(0x123400007FFF);
        let mut slices: Vec<u32> = vec![0; 3];
        msm_slice(scalar, &mut slices, 16);
        assert_eq!(slices.iter().eq([0x7FFF, 0, 0x1234].iter()), true);
    }
}
