use ark_bls12_381::G1Affine;
use ark_ec::AffineCurve;
use ark_ff::PrimeField;
use ark_msm::{
    msm::{compute_msm_baseline, compute_msm_opt},
    utils::generate_msm_inputs,
};

#[cfg(test)]
mod msm_tests {
    use super::*;
    use ark_msm::types::G1BigInt;
    use ark_std::UniformRand;

    #[test]
    fn test_msm_correctness_few_points_no_collision_c14() {
        let num_points = 2;
        let mut rng = ark_std::test_rng();
        let mut points: Vec<_> = Vec::new();
        let mut scalars: Vec<_> = Vec::new();
        for _ in 0..num_points {
            points.push(G1Affine::from(<G1Affine as AffineCurve>::Projective::rand(&mut rng)));
            scalars.push(<G1Affine as AffineCurve>::ScalarField::rand(&mut rng).into_repr());
        }

        let baseline = compute_msm_baseline::<G1Affine>(&points, &scalars);
        let opt = compute_msm_opt(&points, &scalars, 14);
        assert_eq!(baseline, opt);
    }

/*
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

        let baseline = compute_msm_baseline::<G1Affine>(&points, &scalars);
        let opt = compute_msm_opt(&points, &scalars, 15);
        assert_eq!(baseline, opt);
    }
*/
    #[test]
    fn test_msm_correctness_few_points_no_collision_c16() {
        let num_points = 2;
        let mut rng = ark_std::test_rng();
        let mut points: Vec<_> = Vec::new();
        let mut scalars: Vec<_> = Vec::new();
        for _ in 0..num_points {
            points.push(G1Affine::from(<G1Affine as AffineCurve>::Projective::rand(&mut rng)));
            scalars.push(<G1Affine as AffineCurve>::ScalarField::rand(&mut rng).into_repr());
        }

        let baseline = compute_msm_baseline::<G1Affine>(&points, &scalars);
        let opt = compute_msm_opt(&points, &scalars, 16);
        assert_eq!(baseline, opt);
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

        let baseline = compute_msm_baseline::<G1Affine>(&points, &scalars);
        let opt = compute_msm_opt(&points, &scalars, 16);
        assert_eq!(baseline, opt);
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

        let baseline = compute_msm_baseline::<G1Affine>(&points, &scalars);
        let opt = compute_msm_opt(&points, &scalars, 15);
        assert_eq!(baseline, opt);
    }

/*
    #[test]
    fn test_msm_correctness_tremendous_points_c15() {
        let size = 1 << 10;
        let (point_vec, scalar_vec) = generate_msm_inputs::<G1Affine>(size);
        let baseline = compute_msm_baseline::<G1Affine>(&point_vec, &scalar_vec);
        let opt = compute_msm_opt(&point_vec, &scalar_vec, 15);
        assert_eq!(baseline, opt);
    }
*/

    #[test]
    fn test_msm_correctness_tremendous_points_c16() {
        let size = 1 << 10;
        let (point_vec, scalar_vec) = generate_msm_inputs::<G1Affine>(size);
        let baseline = compute_msm_baseline::<G1Affine>(&point_vec, &scalar_vec);
        let opt = compute_msm_opt(&point_vec, &scalar_vec, 8);
        assert_eq!(baseline, opt);
    }
}
