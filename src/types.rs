use ark_bls12_381::G1Affine;
use ark_ec::AffineCurve;
use ark_ff::{FpParameters, PrimeField};

pub const G1_SCALAR_SIZE: u32 =
    <<<G1Affine as AffineCurve>::ScalarField as PrimeField>::Params as FpParameters>::MODULUS_BITS
        as u32;

pub type G1BigInt = <<G1Affine as AffineCurve>::ScalarField as PrimeField>::BigInt;
pub type G1Projective = <G1Affine as AffineCurve>::Projective;
