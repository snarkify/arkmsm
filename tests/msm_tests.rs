use ark_bls12_381::G1Affine;
use ark_ec::msm::VariableBaseMSM as BaselineVariableBaseMSM;
use ark_ec::AffineCurve;
use ark_ff::PrimeField;
use ark_msm::{msm::VariableBaseMSM, utils::generate_msm_inputs};

#[cfg(test)]
mod msm_test {
    use super::*;
    use ark_msm::types::G1BigInt;
    use ark_std::UniformRand;

    fn verify_correctness(points: &[G1Affine], scalars: &[G1BigInt], window_size: u32) {
        let baseline = BaselineVariableBaseMSM::multi_scalar_mul(points, scalars);
        let opt = VariableBaseMSM::multi_scalar_mul_custom(points, scalars, window_size, 2048, 256);
        assert_eq!(baseline, opt);
    }

    #[test]
    fn test_msm_correctness_few_points_no_collision_c14() {
        let num_points = 2;
        let mut rng = ark_std::test_rng();
        let mut points: Vec<_> = Vec::new();
        let mut scalars: Vec<_> = Vec::new();
        for _ in 0..num_points {
            points.push(G1Affine::from(<G1Affine as AffineCurve>::Projective::rand(
                &mut rng,
            )));
            scalars.push(<G1Affine as AffineCurve>::ScalarField::rand(&mut rng).into_repr());
        }

        verify_correctness(&points, &scalars, 14);
    }

    #[test]
    fn test_msm_correctness_few_points_no_collision_c15() {
        let num_points = 2;
        let mut rng = ark_std::test_rng();
        let mut points: Vec<_> = Vec::new();
        let mut scalars: Vec<_> = Vec::new();
        for _ in 0..num_points {
            points.push(G1Affine::from(<G1Affine as AffineCurve>::Projective::rand(
                &mut rng,
            )));
            scalars.push(<G1Affine as AffineCurve>::ScalarField::rand(&mut rng).into_repr());
        }

        verify_correctness(&points, &scalars, 15);
    }

    #[test]
    fn test_msm_correctness_few_points_no_collision_c16() {
        let num_points = 2;
        let mut rng = ark_std::test_rng();
        let mut points: Vec<_> = Vec::new();
        let mut scalars: Vec<_> = Vec::new();
        for _ in 0..num_points {
            points.push(G1Affine::from(<G1Affine as AffineCurve>::Projective::rand(
                &mut rng,
            )));
            scalars.push(<G1Affine as AffineCurve>::ScalarField::rand(&mut rng).into_repr());
        }

        verify_correctness(&points, &scalars, 16);
    }

    #[test]
    fn test_msm_correctness_few_points_u32_all_collision() {
        let num_points = 2;
        let mut rng = ark_std::test_rng();
        let mut points: Vec<_> = Vec::new();
        for _ in 0..num_points {
            points.push(G1Affine::from(<G1Affine as AffineCurve>::Projective::rand(
                &mut rng,
            )));
        }
        let scalars: Vec<_> = vec![G1BigInt::from(0xFFFFFFFF); num_points];

        verify_correctness(&points, &scalars, 14);
    }

    #[test]
    fn test_msm_correctness_few_points_bigint_all_collision() {
        let num_points = 2;
        let mut rng = ark_std::test_rng();
        let mut points: Vec<_> = Vec::new();
        for _ in 0..num_points {
            points.push(G1Affine::from(<G1Affine as AffineCurve>::Projective::rand(
                &mut rng,
            )));
        }
        let scalars: Vec<_> =
            vec![<G1Affine as AffineCurve>::ScalarField::rand(&mut rng).into_repr(); num_points];

        verify_correctness(&points, &scalars, 15);
    }

    #[test]
    fn test_msm_correctness_tremendous_points_c15() {
        let size = 1 << 10;
        let (points, scalars) = generate_msm_inputs::<G1Affine>(size);
        verify_correctness(&points, &scalars, 15);
    }

    #[test]
    fn test_msm_correctness_tremendous_points_c16() {
        let size = 1 << 10;
        let (points, scalars) = generate_msm_inputs::<G1Affine>(size);
        verify_correctness(&points, &scalars, 16);
    }
}
