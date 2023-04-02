#[cfg(test)]
mod glv_tests {
    use ark_msm::{glv::*, types::G1ScalarField};
    use ark_bls12_381::{G1Affine, Fr};
    use ark_ec::AffineCurve;
    use ark_ff::{UniformRand, PrimeField, BigInteger256};
    use num_bigint::BigUint;

    #[test]
    fn test_decompose_1() {
        let s: G1ScalarField = Fr::from(1);
        let (q, r, _, _) = decompose(&s, 16);

        assert_eq!(q, Fr::from(0));
        assert_eq!(r, Fr::from(1));
    }

    #[test]
    fn test_decompose_lambda_plus_1() {
        let s: G1ScalarField = Fr::from(BigInteger256([
            0x0000000100000000,
            0xac45a4010001a402,
            0x0000000000000000,
            0x0000000000000000,
        ]));
        let (q, r, _, _) = decompose(&s, 16);

        assert_eq!(q, Fr::from(1));
        assert_eq!(r, Fr::from(1));
    }

    #[test]
    fn test_decompose_random() {
        let mut rng = ark_std::test_rng();
        let s = Fr::rand(&mut rng);
        let (q, r, _, _) = decompose(&s, 16);

        let lambda: u128 = 0xac45a4010001a40200000000ffffffff;
        let mut ss: BigUint = s.into();
        let modulus: BigUint = BigUint::parse_bytes(b"73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16).unwrap();
        if s.into_repr().as_ref()[3] >= 0x3FFFFFFFF {
            ss = modulus - ss;
        }
        let quotient = ss.clone() / lambda;
        let remainder = ss - quotient.clone() * lambda;
        assert_eq!(q, Fr::from(quotient));
        assert_eq!(r, Fr::from(remainder));
    }

    #[test]
    fn test_endo() {
        let mut rng = ark_std::test_rng();
        let lambda: G1ScalarField = Fr::from(BigInteger256([
            0x00000000ffffffff,
            0xac45a4010001a402,
            0x0000000000000000,
            0x0000000000000000,
        ]));
        let p = G1Affine::from(<G1Affine as AffineCurve>::Projective::rand(&mut rng));

        let mut endo_p = p.clone();
        endomorphism(&mut endo_p);

        let lambda_p = G1Affine::from(p.mul(lambda));
        assert_eq!(lambda_p, endo_p);
    }
}
