// use ark_bls12_381::G1Affine;
use ark_ec::{models::ModelParameters, AffineCurve};
use ark_ff::{FpParameters, PrimeField};

// pub const G1_SCALAR_SIZE: u32 =
//     <<<G1Affine as AffineCurve>::ScalarField as PrimeField>::Params as FpParameters>::MODULUS_BITS;
pub const G1_SCALAR_SIZE_GLV: u32 = 128u32;
pub const GROUP_SIZE_IN_BITS: usize = 6;
pub const GROUP_SIZE: usize = 1 << GROUP_SIZE_IN_BITS;

// pub type G1BigInt = <<G1Affine as AffineCurve>::ScalarField as PrimeField>::BigInt;
// pub type G1Projective = <G1Affine as AffineCurve>::Projective;
// pub type G1ScalarField = <G1Affine as AffineCurve>::ScalarField;
// pub type G1BaseField = <G1Affine as AffineCurve>::BaseField;

pub type BigInt<P> = <<P as ModelParameters>::ScalarField as PrimeField>::BigInt;
