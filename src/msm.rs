use crate::{
    bucket_msm::BucketMSM,
    glv::{decompose},
    types::{G1BigInt, G1Projective, G1_SCALAR_SIZE_GLV},
};
use ark_bls12_381::{G1Affine, Fr};
use ark_ff::BigInteger;

pub struct VariableBaseMSM;

impl VariableBaseMSM {
    fn msm_slice(scalar: G1BigInt, slices: &mut Vec<u32>, window_bits: u32) {
        assert!(window_bits <= 31); // reserve one bit for marking signed slices
        let mut temp = scalar;
        for i in 0..slices.len() {
            slices[i] = (temp.as_ref()[0] % (1 << window_bits)) as u32;
            temp.divn(window_bits);
        }

        let mut carry = 0;
        let total = 1 << window_bits;
        let half = total >> 1;
        for i in 0..slices.len() {
            slices[i] += carry;
            if slices[i] > half {
                // slices[i] == half is okay, since (slice[i]-1) will be used for bucket_id
                slices[i] = total - slices[i];
                carry = 1;
                slices[i] |= 1 << 31; // mark the highest bit for later
            } else {
                carry = 0;
            }
        }
        assert!(
            carry == 0,
            "msm_slice overflows when apply signed-bucket-index"
        );
    }

    pub fn multi_scalar_mul_custom(
        points: &Vec<G1Affine>,
        scalars: &Vec<G1BigInt>,
        window_bits: u32,
        max_batch: u32,
        max_collisions: u32
    ) -> G1Projective {
        let mut bucket_msm = BucketMSM::new(
            G1_SCALAR_SIZE_GLV,
            window_bits,
            max_batch,
            max_collisions,
        );
        let num_slices: u32 = (G1_SCALAR_SIZE_GLV + window_bits - 1) / window_bits;
        // scalar = phi * lambda + normal
        let mut phi_slices: Vec<u32> = vec![0; num_slices as usize];
        let mut normal_slices: Vec<u32> = vec![0; num_slices as usize];

        let scalars_and_bases_iter = scalars
            .iter()
            .zip(points)
            .filter(|(s, _)| !s.is_zero());
        scalars_and_bases_iter.for_each(|(&scalar, point)| {
            let (phi, normal, is_neg_scalar, is_neg_normal) = decompose(&Fr::from(scalar), window_bits);
            Self::msm_slice(phi.into(), &mut phi_slices, window_bits);
            Self::msm_slice(normal.into(), &mut normal_slices, window_bits);
            bucket_msm.process_point_and_slices_glv(&point, &normal_slices, &phi_slices, is_neg_scalar, is_neg_normal);
        });

        bucket_msm.process_complete();
        return bucket_msm.batch_reduce();
    }
}


#[cfg(test)]
mod collision_method_pippenger_tests {
    use super::*;

    #[test]
    fn test_msm_slice_window_size_1() {
        let scalar = G1BigInt::from(0b101);
        let mut slices: Vec<u32> = vec![0; 3];
        VariableBaseMSM::msm_slice(scalar, &mut slices, 1);
        // print!("slices {:?}\n", slices);
        assert_eq!(slices.iter().eq([1, 0, 1].iter()), true);
    }
    #[test]
    fn test_msm_slice_window_size_2() {
        let scalar = G1BigInt::from(0b000110);
        let mut slices: Vec<u32> = vec![0; 3];
        VariableBaseMSM::msm_slice(scalar, &mut slices, 2);
        assert_eq!(slices.iter().eq([2, 1, 0].iter()), true);
    }

    #[test]
    fn test_msm_slice_window_size_3() {
        let scalar = G1BigInt::from(0b010111000);
        let mut slices: Vec<u32> = vec![0; 3];
        VariableBaseMSM::msm_slice(scalar, &mut slices, 3);
        assert_eq!(slices.iter().eq([0, 0x80000001, 3].iter()), true);
    }

    #[test]
    fn test_msm_slice_window_size_16() {
        let scalar = G1BigInt::from(0x123400007FFF);
        let mut slices: Vec<u32> = vec![0; 3];
        VariableBaseMSM::msm_slice(scalar, &mut slices, 16);
        assert_eq!(slices.iter().eq([0x7FFF, 0, 0x1234].iter()), true);
    }
}
