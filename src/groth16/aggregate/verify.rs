use digest::Digest;
use ff::{Field, PrimeField};
use groupy::{CurveAffine, CurveProjective};

use super::{
    inner_product, prove::polynomial_evaluation_product_form_from_transcript,
    structured_scalar_power, AggregateProof, GIPAProof, GIPAProofWithSSM,
    MultiExpInnerProductCProof, PairingInnerProductABProof, VerifierSRS,
};
use crate::bls::Engine;
use crate::groth16::VerifyingKey;

pub fn verify_aggregate_proof<E: Engine, D: Digest>(
    ip_verifier_srs: &VerifierSRS<E>,
    vk: &VerifyingKey<E>,
    public_inputs: &Vec<Vec<E::Fr>>,
    proof: &AggregateProof<E, D>,
) -> bool {
    // Random linear combination of proofs
    let mut counter_nonce: usize = 0;
    let r = loop {
        let mut hash_input = Vec::new();
        hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
        hash_input.extend_from_slice(&proof.com_a.as_bytes());
        hash_input.extend_from_slice(&proof.com_b.as_bytes());
        hash_input.extend_from_slice(&proof.com_c.as_bytes());
        if let Some(r) =
            E::Fr::from_bytes(&D::digest(&hash_input).as_slice()[..E::Fr::SERIALIZED_BYTES])
        {
            break r;
        };
        counter_nonce += 1;
    };

    // Check TIPA proofs
    let tipa_proof_ab_valid = verify_with_srs_shift::<E, D>(
        ip_verifier_srs,
        (&proof.com_a, &proof.com_b, &vec![proof.ip_ab]),
        &proof.tipa_proof_ab,
        &r,
    );
    let tipa_proof_c_valid = verify_with_structured_scalar_message::<E, D>(
        ip_verifier_srs,
        (&proof.com_c, &vec![proof.agg_c]),
        &r,
        &proof.tipa_proof_c,
    );
    // assert!(tipa_proof_ab_valid, "TIPP failed.");
    // assert!(tipa_proof_c_valid, "MIPP failed.");

    // Check aggregate pairing product equation
    // {
    // let mut r_sum2 = r.pow(&[public_inputs.len() as u64]);
    // r_sum2.sub_assign(&E::Fr::one());
    // let mut r_div = r.clone();
    // r_div.sub_assign(&E::Fr::one());
    // r_div.inverse();
    // r_sum.mul_assign(&r_div);
    // }
    let mut r_sum = E::Fr::zero();
    for j in 0..public_inputs.len() {
        r_sum.add_assign(&r.clone().pow(&[2 * j as u64]));
    }

    let p1 = {
        let mut alpha_g1_r_sum = vk.alpha_g1.into_projective();
        alpha_g1_r_sum.mul_assign(r_sum);
        E::pairing(alpha_g1_r_sum, vk.beta_g2)
    };

    assert_eq!(vk.ic.len(), public_inputs[0].len() + 1);
    let r_vec = structured_scalar_power(public_inputs.len(), &r);
    let mut g_ic = vk.ic[0].into_projective();
    g_ic.mul_assign(r_sum);

    for (i, b) in vk.ic.iter().skip(1).enumerate() {
        let mut x = b.into_projective();
        x.mul_assign(inner_product::scalar(
            &public_inputs
                .iter()
                .map(|inputs| inputs[i].clone())
                .collect::<Vec<E::Fr>>(),
            &r_vec,
        ));
        g_ic.add_assign(&x);
    }
    let p2 = E::pairing(g_ic, vk.gamma_g2);
    let p3 = E::pairing(proof.agg_c, vk.delta_g2);

    let mut p1_p2_p3 = p1;
    p1_p2_p3.mul_assign(&p2);
    p1_p2_p3.mul_assign(&p3);
    let ppe_valid = proof.ip_ab == p1_p2_p3;
    assert!(ppe_valid, "outer verification failed");

    tipa_proof_ab_valid && tipa_proof_c_valid && ppe_valid
}

fn verify_with_srs_shift<E: Engine, D: Digest>(
    v_srs: &VerifierSRS<E>,
    com: (&E::Fqk, &E::Fqk, &Vec<E::Fqk>),
    proof: &PairingInnerProductABProof<E, D>,
    r_shift: &E::Fr,
) -> bool {
    let (base_com, transcript) = gipa_verify_recursive_challenge_transcript(com, &proof.gipa_proof);
    let transcript_inverse = transcript.iter().map(|x| x.inverse().unwrap()).collect();

    // Verify commitment keys wellformed
    let (ck_a_final, ck_b_final) = &proof.final_ck;
    let (ck_a_proof, ck_b_proof) = &proof.final_ck_proof;

    // KZG challenge point
    let mut counter_nonce: usize = 0;
    let c = loop {
        let mut hash_input = Vec::new();
        hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
        //TODO: Should use CanonicalSerialize instead of ToBytes
        hash_input.extend_from_slice(&transcript.first().unwrap().as_bytes());
        hash_input.extend_from_slice(ck_a_final.into_affine().into_uncompressed().as_ref());
        hash_input.extend_from_slice(&ck_b_final.into_affine().into_uncompressed().as_ref());

        if let Some(c) =
            E::Fr::from_bytes(&D::digest(&hash_input).as_slice()[..E::Fr::SERIALIZED_BYTES])
        {
            break c;
        };
        counter_nonce += 1;
    };

    let ck_a_valid = verify_commitment_key_g2_kzg_opening(
        v_srs,
        &ck_a_final,
        &ck_a_proof,
        &transcript_inverse,
        &r_shift.inverse().unwrap(),
        &c,
    );
    let ck_b_valid = verify_commitment_key_g1_kzg_opening(
        v_srs,
        &ck_b_final,
        &ck_b_proof,
        &transcript,
        &E::Fr::one(),
        &c,
    );

    // Verify base inner product commitment
    let (com_a, com_b, _com_t) = base_com;
    let a_base = vec![proof.gipa_proof.r_base.0.clone()];
    let b_base = vec![proof.gipa_proof.r_base.1.clone()];
    let _t_base = vec![inner_product::pairing::<E>(&a_base, &b_base)];

    let base_valid = inner_product::pairing::<E>(&a_base, &vec![ck_a_final.clone()]) == com_a
        && inner_product::pairing::<E>(&vec![ck_b_final.clone()], &b_base) == com_b;

    ck_a_valid && ck_b_valid && base_valid
}

fn gipa_verify_recursive_challenge_transcript<E: Engine, D: Digest>(
    com: (&E::Fqk, &E::Fqk, &Vec<E::Fqk>),
    proof: &GIPAProof<E, D>,
) -> ((E::Fqk, E::Fqk, Vec<E::Fqk>), Vec<E::Fr>) {
    let (com_0, com_1, com_2) = com.clone();
    let (mut com_a, mut com_b, mut com_t) = (*com_0, *com_1, com_2.clone());
    let mut r_transcript = Vec::new();
    for (com_1, com_2) in proof.r_commitment_steps.iter().rev() {
        // Fiat-Shamir challenge
        let mut counter_nonce: usize = 0;
        let default_transcript = E::Fr::zero();
        let transcript = r_transcript.last().unwrap_or(&default_transcript);
        let (c, c_inv) = 'challenge: loop {
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
            hash_input.extend_from_slice(&transcript.as_bytes());
            hash_input.extend_from_slice(&com_1.0.as_bytes());
            hash_input.extend_from_slice(&com_1.1.as_bytes());
            for c in &com_1.2 {
                hash_input.extend_from_slice(&c.as_bytes());
            }
            hash_input.extend_from_slice(&com_2.0.as_bytes());
            hash_input.extend_from_slice(&com_2.1.as_bytes());
            for c in &com_2.2 {
                hash_input.extend_from_slice(&c.as_bytes());
            }
            let c =
                E::Fr::from_bytes(&D::digest(&hash_input).as_slice()[..E::Fr::SERIALIZED_BYTES]);
            if let Some(c) = c {
                if let Some(c_inv) = c.inverse() {
                    // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                    // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                    break 'challenge (c_inv, c);
                }
            }
            counter_nonce += 1;
        };

        com_a = add!(
            add!(mul_scalar::<E>(com_1.0, &c), &com_a),
            &mul_scalar::<E>(com_2.0, &c_inv)
        );
        com_b = add!(
            add!(mul_scalar::<E>(com_1.1, &c), &com_b),
            &mul_scalar::<E>(com_2.1, &c_inv)
        );
        com_t = {
            let x = com_2
                .2
                .iter()
                .map(|com_2_2| mul_scalar::<E>(*com_2_2, &c_inv));
            let a = com_1.2.iter().map(|com_1_2| mul_scalar::<E>(*com_1_2, &c));
            let b = a.zip(com_t.iter()).map(|(a, com_t)| add!(a, com_t));
            b.zip(x).map(|(b, x)| add!(b, &x)).collect::<Vec<_>>()
        };

        r_transcript.push(c);
    }
    r_transcript.reverse();
    ((com_a, com_b, com_t), r_transcript)
}

pub fn verify_commitment_key_g2_kzg_opening<E: Engine>(
    v_srs: &VerifierSRS<E>,
    ck_final: &E::G2,
    ck_opening: &E::G2,
    transcript: &Vec<E::Fr>,
    r_shift: &E::Fr,
    kzg_challenge: &E::Fr,
) -> bool {
    let ck_polynomial_c_eval =
        polynomial_evaluation_product_form_from_transcript(transcript, kzg_challenge, r_shift);
    E::pairing(
        v_srs.g.clone(),
        sub!(*ck_final, &mul!(v_srs.h, ck_polynomial_c_eval)),
    ) == E::pairing(
        sub!(v_srs.g_beta, &mul!(v_srs.g, kzg_challenge.clone())),
        ck_opening.clone(),
    )
}

pub fn verify_commitment_key_g1_kzg_opening<E: Engine>(
    v_srs: &VerifierSRS<E>,
    ck_final: &E::G1,
    ck_opening: &E::G1,
    transcript: &Vec<E::Fr>,
    r_shift: &E::Fr,
    kzg_challenge: &E::Fr,
) -> bool {
    let ck_polynomial_c_eval =
        polynomial_evaluation_product_form_from_transcript(transcript, kzg_challenge, r_shift);
    E::pairing(
        sub!(*ck_final, &mul!(v_srs.g, ck_polynomial_c_eval)),
        v_srs.h,
    ) == E::pairing(
        *ck_opening,
        sub!(v_srs.h_alpha, &mul!(v_srs.h, *kzg_challenge)),
    )
}

fn mul_scalar<E: Engine>(a: E::Fqk, b: &E::Fr) -> E::Fqk {
    a.pow(b.into_repr())
}

fn verify_with_structured_scalar_message<E: Engine, D: Digest>(
    v_srs: &VerifierSRS<E>,
    com: (&E::Fqk, &Vec<E::G1>),
    scalar_b: &E::Fr,
    proof: &MultiExpInnerProductCProof<E, D>,
) -> bool {
    let (base_com, transcript) = gipa_with_ssm_verify_recursive_challenge_transcript(
        (com.0, scalar_b, com.1),
        &proof.gipa_proof,
    );
    let transcript_inverse = transcript.iter().map(|x| x.inverse().unwrap()).collect();

    let ck_a_final = &proof.final_ck;
    let ck_a_proof = &proof.final_ck_proof;

    // KZG challenge point
    let mut counter_nonce: usize = 0;
    let c = loop {
        let mut hash_input = Vec::new();
        hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
        //TODO: Should use CanonicalSerialize instead of ToBytes
        hash_input.extend_from_slice(&transcript.first().unwrap().as_bytes());
        hash_input.extend_from_slice(ck_a_final.into_affine().into_uncompressed().as_ref());
        if let Some(c) =
            E::Fr::from_bytes(&D::digest(&hash_input).as_slice()[..E::Fr::SERIALIZED_BYTES])
        {
            break c;
        };
        counter_nonce += 1;
    };

    // Check commitment key
    let ck_a_valid = verify_commitment_key_g2_kzg_opening(
        v_srs,
        &ck_a_final,
        &ck_a_proof,
        &transcript_inverse,
        &E::Fr::one(),
        &c,
    );

    // Compute final scalar
    let mut power_2_b = scalar_b.clone();
    let mut b_base = E::Fr::one();
    for x in transcript.iter() {
        b_base.mul_assign(&add!(
            E::Fr::one(),
            &(mul!(x.inverse().unwrap(), &power_2_b))
        ));
        power_2_b.mul_assign(&power_2_b.clone());
    }

    // Verify base inner product commitment
    let (com_a, com_t) = base_com;
    let a_base = vec![proof.gipa_proof.r_base.0.clone()];
    let t_base = vec![inner_product::multiexponentiation(&a_base, &vec![b_base])];
    let base_valid = inner_product::pairing::<E>(&a_base, &vec![ck_a_final.clone()]) == com_a
        && &t_base == &com_t;

    ck_a_valid && base_valid
}

fn gipa_with_ssm_verify_recursive_challenge_transcript<E: Engine, D: Digest>(
    com: (&E::Fqk, &E::Fr, &Vec<E::G1>),
    proof: &GIPAProofWithSSM<E, D>,
) -> ((E::Fqk, Vec<E::G1>), Vec<E::Fr>) {
    let (com_0, com_1, com_2) = com.clone();
    let (mut com_a, _com_b, mut com_t) = (*com_0, *com_1, com_2.clone());
    let mut r_transcript = Vec::new();
    for (com_1, com_2) in proof.r_commitment_steps.iter().rev() {
        // Fiat-Shamir challenge
        let mut counter_nonce: usize = 0;
        let default_transcript = E::Fr::zero();
        let transcript = r_transcript.last().unwrap_or(&default_transcript);
        let (c, c_inv) = 'challenge: loop {
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
            hash_input.extend_from_slice(&transcript.as_bytes());
            hash_input.extend_from_slice(&com_1.0.as_bytes());
            for c in &com_1.1 {
                hash_input.extend_from_slice(c.into_affine().into_uncompressed().as_ref());
            }
            hash_input.extend_from_slice(&com_2.0.as_bytes());
            for c in &com_2.1 {
                hash_input.extend_from_slice(c.into_affine().into_uncompressed().as_ref());
            }
            let c =
                E::Fr::from_bytes(&D::digest(&hash_input).as_slice()[..E::Fr::SERIALIZED_BYTES]);
            if let Some(c) = c {
                if let Some(c_inv) = c.inverse() {
                    // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                    // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                    break 'challenge (c_inv, c);
                }
            }
            counter_nonce += 1;
        };

        com_a = add!(
            add!(mul_scalar::<E>(com_1.0, &c), &com_a),
            &mul_scalar::<E>(com_2.0, &c_inv)
        );

        com_t = {
            let x = com_2.1.iter().map(|com_2_1| mul!(*com_2_1, c_inv));
            let a = com_1.1.iter().map(|com_1_1| mul!(*com_1_1, c));
            let b = a.zip(com_t.iter()).map(|(a, com_t)| add!(a, com_t));
            b.zip(x).map(|(b, x)| add!(b, &x)).collect::<Vec<_>>()
        };

        r_transcript.push(c);
    }
    r_transcript.reverse();
    ((com_a, com_t), r_transcript)
}
