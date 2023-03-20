use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::{PrimeField, UniformRand};

pub fn generate_msm_inputs<A>(
    size: usize,
) -> (
    Vec<<A::Projective as ProjectiveCurve>::Affine>,
    Vec<<A::ScalarField as PrimeField>::BigInt>,
)
where
    A: AffineCurve,
{
    let mut rng = ark_std::test_rng();
    let scalar_vec = (0..size)
        .map(|_| A::ScalarField::rand(&mut rng).into_repr())
        .collect();
    let point_vec = (0..size)
        .map(|_| A::Projective::rand(&mut rng))
        .collect::<Vec<_>>();
    (
        <A::Projective as ProjectiveCurve>::batch_normalization_into_affine(&point_vec),
        scalar_vec,
    )
}
