use halo2_curves::bn256::{pairing, Bn256, Fr, G1Affine};
use merlin::Transcript;

type Scalar = Fr;

#[derive(Clone, Debug)]
pub struct MvKzgProof {
    pub w: Vec<Option<G1Affine>>,
}

#[derive(Clone, Debug)]
pub struct GrandProductData {
    pub f_0_r: Scalar,
    pub f_1_r: Scalar,
    pub f_r_0: Scalar,
    pub f_r_1: Scalar,
    pub product: Scalar,
    pub f_0_r_proof: G1Affine,
    pub f_1_r_proof: G1Affine,
    pub f_r_0_proof: G1Affine,
    pub f_r_1_proof: G1Affine,
    pub product_proof: G1Affine,
}

/// Lasso proof message round 1
#[derive(Clone, Debug)]
pub struct Message1 {
    pub a_comm: G1Affine,
    pub logm: u64,
    pub dim_comm: Vec<G1Affine>,
}

/// Lasso proof message round 2
#[derive(Clone, Debug)]
pub struct Message2 {
    pub a_eval: Scalar,
    pub a_eval_proof: MvKzGProof,
    pub e_comm: Vec<G1Affine>,
    pub read_ts_comm: Vec<G1Affine>,
    pub final_cts_comm: Vec<G1Affine>,
}

/// Lasso proof message round 3
#[derive(Clone, Debug)]
pub struct Message3 {
    pub h_sumcheck_proof: Vec<Vec<Scalar>>,
    pub rz: Vec<Scalar>,
    pub e_eval: Vec<Scalar>,
    pub e_eval_proof: Vec<MvKzGProof>,
}

/// Lasso proof message round 4
#[derive(Clone, Debug)]
pub struct Message4 {
    pub s0_comm: Vec<G1Affine>,
    pub s_comm: Vec<G1Affine>,
    pub rs_comm: Vec<G1Affine>,
    pub ws_comm: Vec<G1Affine>,
}

/// Lasso proof message round 5
#[derive(Clone, Debug)]
pub struct Message5 {
    pub s0_sumcheck_proof: Vec<Vec<Vec<Scalar>>>,
    pub s_sumcheck_proof: Vec<Vec<Vec<Scalar>>>,
    pub rs_sumcheck_proof: Vec<Vec<Vec<Scalar>>>,
    pub ws_sumcheck_proof: Vec<Vec<Vec<Scalar>>>,
    pub r_prime: Vec<Vec<Scalar>>,
    pub r_prime2: Vec<Vec<Scalar>>,
    pub r_prime3: Vec<Vec<Scalar>>,
    pub r_prime4: Vec<Vec<Scalar>>,
    pub s0_data: Vec<GrandProductData>,
    pub s_data: Vec<GrandProductData>,
    pub rs_data: Vec<GrandProductData>,
    pub ws_data: Vec<GrandProductData>,
    pub e_eval2: Vec<Scalar>,
    pub dim_eval: Vec<Scalar>,
    pub read_ts_eval: Vec<Scalar>,
    pub final_cts_eval: Vec<Scalar>,
    pub e_eval2_proof: Vec<MvKzGProof>,
    pub dim_eval_proof: Vec<MvKzGProof>,
    pub read_ts_eval_proof: Vec<MvKzGProof>,
    pub final_cts_eval_proof: Vec<MvKzGProof>,
}
