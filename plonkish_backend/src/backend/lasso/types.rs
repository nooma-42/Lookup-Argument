use halo2_curves::bn256::{Fr, G1Affine};
use std::collections::HashMap;

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
    pub logm: usize,
    pub dim_comm: Vec<G1Affine>,
}

/// Lasso proof message round 2
#[derive(Clone, Debug)]
pub struct Message2 {
    pub a_eval: Scalar,
    pub a_eval_proof: MvKzgProof,
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
    pub e_eval_proof: Vec<MvKzgProof>,
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
    pub e_eval2_proof: Vec<MvKzgProof>,
    pub dim_eval_proof: Vec<MvKzgProof>,
    pub read_ts_eval_proof: Vec<MvKzgProof>,
    pub final_cts_eval_proof: Vec<MvKzgProof>,
}

#[derive(Clone, Debug)]
pub struct Proof {
    pub msg_1: Message1,
    pub msg_2: Message2,
    pub msg_3: Message3,
    pub msg_4: Message4,
    pub msg_5: Message5,
}

#[derive(Clone, Debug)]
pub enum FlattenItem {
    Scalar(Scalar),
    Scalars(Vec<Scalar>),
    Point(G1Affine),
    Points(Vec<G1Affine>),
    Proof(MvKzgProof),
    Proofs(Vec<MvKzgProof>),
    SumcheckProof(Vec<Vec<Scalar>>),
    SumcheckProofs(Vec<Vec<Vec<Scalar>>>),
    GrandData(Vec<GrandProductData>),
    Usize(usize),
}

#[rustfmt::skip]
impl Proof {
    pub fn flatten(&self) -> HashMap<String, FlattenItem> {
        let mut proof = HashMap::new();

        // msg_1
        proof.insert("a_comm".to_string(), FlattenItem::Point(self.msg_1.a_comm.clone()));
        proof.insert("logm".to_string(), FlattenItem::Usize(self.msg_1.logm));
        proof.insert("dim_comm".to_string(), FlattenItem::Points(self.msg_1.dim_comm.clone()));

        // msg_2
        proof.insert("a_eval".to_string(), FlattenItem::Scalar(self.msg_2.a_eval.clone()));
        proof.insert("a_eval_proof".to_string(), FlattenItem::Proof(self.msg_2.a_eval_proof.clone()));
        proof.insert("E_comm".to_string(), FlattenItem::Points(self.msg_2.e_comm.clone()));
        proof.insert("read_ts_comm".to_string(), FlattenItem::Points(self.msg_2.read_ts_comm.clone()));
        proof.insert("final_cts_comm".to_string(), FlattenItem::Points(self.msg_2.final_cts_comm.clone()));

        // msg_3
        proof.insert("h_sumcheck_proof".to_string(), FlattenItem::SumcheckProof(self.msg_3.h_sumcheck_proof.clone()));
        proof.insert("rz".to_string(), FlattenItem::Scalars(self.msg_3.rz.clone()));
        proof.insert("E_eval".to_string(), FlattenItem::Scalars(self.msg_3.e_eval.clone()));
        proof.insert("E_eval_proof".to_string(), FlattenItem::Proofs(self.msg_3.e_eval_proof.clone()));

        // msg_4
        proof.insert("S0_comm".to_string(), FlattenItem::Points(self.msg_4.s0_comm.clone()));
        proof.insert("S_comm".to_string(), FlattenItem::Points(self.msg_4.s_comm.clone()));
        proof.insert("RS_comm".to_string(), FlattenItem::Points(self.msg_4.rs_comm.clone()));
        proof.insert("WS_comm".to_string(), FlattenItem::Points(self.msg_4.ws_comm.clone()));

        // msg_5
        proof.insert("S0_sumcheck_proof".to_string(), FlattenItem::SumcheckProofs(self.msg_5.s0_sumcheck_proof.clone()));
        proof.insert("S_sumcheck_proof".to_string(), FlattenItem::SumcheckProofs(self.msg_5.s_sumcheck_proof.clone()));
        proof.insert("RS_sumcheck_proof".to_string(), FlattenItem::SumcheckProofs(self.msg_5.rs_sumcheck_proof.clone()));
        proof.insert("WS_sumcheck_proof".to_string(), FlattenItem::SumcheckProofs(self.msg_5.ws_sumcheck_proof.clone()));
        proof.insert("r_prime".to_string(), FlattenItem::SumcheckProof(self.msg_5.r_prime.clone()));
        proof.insert("r_prime2".to_string(), FlattenItem::SumcheckProof(self.msg_5.r_prime2.clone()));
        proof.insert("r_prime3".to_string(), FlattenItem::SumcheckProof(self.msg_5.r_prime3.clone()));
        proof.insert("r_prime4".to_string(), FlattenItem::SumcheckProof(self.msg_5.r_prime4.clone()));

        proof.insert("S0_data".to_string(), FlattenItem::GrandData(self.msg_5.s0_data.clone()));
        proof.insert("S_data".to_string(), FlattenItem::GrandData(self.msg_5.s_data.clone()));
        proof.insert("RS_data".to_string(), FlattenItem::GrandData(self.msg_5.rs_data.clone()));
        proof.insert("WS_data".to_string(), FlattenItem::GrandData(self.msg_5.ws_data.clone()));

        proof.insert("E_eval2".to_string(), FlattenItem::Scalars(self.msg_5.e_eval2.clone()));
        proof.insert("dim_eval".to_string(), FlattenItem::Scalars(self.msg_5.dim_eval.clone()));
        proof.insert("read_ts_eval".to_string(), FlattenItem::Scalars(self.msg_5.read_ts_eval.clone()));
        proof.insert("final_cts_eval".to_string(), FlattenItem::Scalars(self.msg_5.final_cts_eval.clone()));

        proof.insert("E_eval2_proof".to_string(), FlattenItem::Proofs(self.msg_5.e_eval2_proof.clone()));
        proof.insert("dim_eval_proof".to_string(), FlattenItem::Proofs(self.msg_5.dim_eval_proof.clone()));
        proof.insert("read_ts_eval_proof".to_string(), FlattenItem::Proofs(self.msg_5.read_ts_eval_proof.clone()));
        proof.insert("final_cts_eval_proof".to_string(), FlattenItem::Proofs(self.msg_5.final_cts_eval_proof.clone()));

        return proof
    }
}
