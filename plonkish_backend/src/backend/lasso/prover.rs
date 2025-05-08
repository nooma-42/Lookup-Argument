use crate::backend::lasso::{
    mle_poly::*,
    sostable::SOSTable,
    types::*,
    util::{hash_tuple, log_ceil},
};
use crate::pcs::{
    multilinear::{
        MultilinearKzg, MultilinearKzgCommitment, MultilinearKzgParam, MultilinearKzgProverParam,
    },
    PolynomialCommitmentScheme,
};
use crate::piop::sum_check::{
    classic::{ClassicSumCheck, EvaluationsProver},
    SumCheck, VirtualPolynomial,
};
use crate::poly::multilinear::MultilinearPolynomial;
use crate::util::{
    expression::{rotate::BinaryField, Expression, Query},
    hash::Hash,
    transcript::{FiatShamirTranscript, InMemoryTranscript, Keccak256Transcript},
};
use halo2_curves::bn256::{Bn256, Fr, G1Affine};
use std::collections::BTreeMap;
use std::io;

use super::types::MvKzgProof;

type Scalar = Fr;
type Pcs = MultilinearKzg<Bn256>;

pub struct Prover<'b> {
    table: &'b SOSTable,
    param: &'b MultilinearKzgParam<Bn256>,
    pp: &'b MultilinearKzgProverParam<Bn256>,
    m: usize,
    logm: usize,
    r: Vec<Scalar>,
    tau: Scalar,
    gamma: Scalar,
}

impl Prover<'_> {
    pub fn new<'a>(
        table: &'a SOSTable,
        param: &'a MultilinearKzgParam<Bn256>,
        pp: &'a MultilinearKzgProverParam<Bn256>,
    ) -> Prover<'a> {
        Prover {
            table,
            param,
            pp,
            m: 0,
            logm: 0,
            r: vec![],
            tau: Scalar::zero(),
            gamma: Scalar::zero(),
        }
    }
    pub fn prove(&mut self, witness: &Vec<u64>) -> Proof {
        // [round 0]
        let mut transcript = Keccak256Transcript::new(());

        // [round 1]
        self.m = witness.len();
        self.logm = log_ceil(self.m);
        assert!(self.logm <= self.table.l);

        let mut a = vec![Scalar::zero(); 2 << self.logm];
        for i in 0..self.m {
            a[i] = Scalar::from(witness[i]);
        }

        // all commits in prove() are different to the ones in pylookup;
        let a_poly: MultilinearPolynomial<Scalar> = MultilinearPolynomial::new(a.clone());
        let a_comm: MultilinearKzgCommitment<G1Affine> =
            Pcs::commit_and_write(&self.pp, &a_poly, &mut transcript).unwrap();
        let mut indice = Vec::new();
        for w in witness {
            indice.push(self.table.get_index(*w));
        }
        let mut dim_ploy: Vec<MultilinearPolynomial<Scalar>> = Vec::new();
        let mut dim_values: Vec<Vec<Scalar>> = Vec::new();
        for i in 0..self.table.c {
            let mut values = vec![0; 2 << self.logm];
            for j in 0..self.m {
                values[j] = indice[j][i];
            }
            self.append_poly_and_values(&mut dim_ploy, &mut dim_values, &values, self.logm);
        }
        let dim_comm = Pcs::batch_commit_and_write(self.pp, &dim_ploy, &mut transcript).unwrap();
        let msg_1 = Message1::new(a_comm, self.logm, dim_comm);

        // [Round 2]
        let (a_eval, a_eval_proof) = self.eval_and_prove(&a_poly, &self.r, &mut transcript);
        let mut e_poly = Vec::new();
        let mut read_poly = Vec::new();
        let mut write_poly = Vec::new();
        let mut final_poly = Vec::new();

        let mut e_values = Vec::new();
        let mut read_values = Vec::new();
        let mut write_values = Vec::new();
        let mut final_values = Vec::new();
        for i in 0..self.table.alpha {
            let mut e_val = vec![0; 2 << self.logm];
            let mut read_ts = vec![0; 2 << self.logm];
            let mut write_cts = vec![0; 2 << self.logm];
            let mut final_cts = vec![0; 2 << self.logm];
            for j in 0..self.m {
                let index = indice[j][i / self.table.k];
                e_val[j] = self.table.tables[i][index as usize];
                let ts = final_cts[j];
                read_ts[j] = ts;
                write_cts[j] = ts + 1;
                final_cts[index as usize] = ts + 1;
            }
            self.append_poly_and_values(&mut e_poly, &mut e_values, &e_val, self.logm);
            self.append_poly_and_values(&mut read_poly, &mut read_values, &read_ts, self.logm);
            self.append_poly_and_values(&mut write_poly, &mut write_values, &write_cts, self.logm);
            self.append_poly_and_values(
                &mut final_poly,
                &mut final_values,
                &final_cts,
                self.table.l,
            );
        }
        let e_comm = Pcs::batch_commit_and_write(self.pp, &e_poly, &mut transcript).unwrap();

        let read_comm = Pcs::batch_commit_and_write(self.pp, &read_poly, &mut transcript).unwrap();
        let final_comm =
            Pcs::batch_commit_and_write(self.pp, &final_poly, &mut transcript).unwrap();
        let msg_2 = Message2::new(a_eval, a_eval_proof, e_comm, read_comm, final_comm);

        // [round 3]
        let p1 = eq_mle_poly(&self.r).to_mle();
        let p2 = self.table.g_func(&e_poly);
        let h_poly = mle_mul(&p1, &p2);
        let (_, rz, h_sumcheck_proof) = self.prove_sumcheck(&h_poly, &mut transcript);

        let mut e_eval = Vec::new();
        let mut e_eval_proof = Vec::new();
        for poly in &e_poly {
            self.append_eval_and_proof(&mut e_eval, &mut e_eval_proof, poly, &rz, &mut transcript);
        }
        let msg_3 = Message3::new(h_sumcheck_proof, rz, e_eval, e_eval_proof);

        // [round 4]
        let mut s0_polys = Vec::new();
        let mut s_polys = Vec::new();
        let mut rs_polys = Vec::new();
        let mut ws_polys = Vec::new();
        let mut s0_comms = Vec::new();
        let mut s_comms = Vec::new();
        let mut rs_comms = Vec::new();
        let mut ws_comms = Vec::new();
        for i in 0..self.table.alpha {
            let mut s0 = Vec::new();
            let mut s = Vec::new();
            let mut rs = Vec::new();
            let mut ws = Vec::new();
            let range_l = 2 << self.table.l;
            let range_m = 2 << self.logm;
            for j in 0..range_l {
                s0.push((
                    Scalar::from(j),
                    Scalar::from(self.table.tables[i][j as usize]), // should know if j is not bigger than usize
                    Scalar::zero(),
                ));
                s.push((
                    Scalar::from(j),
                    Scalar::from(self.table.tables[i][j as usize]),
                    final_values[i][j as usize],
                ))
            }
            for j in 0..range_m {
                rs.push((dim_values[i][j], e_values[i][j], read_values[i][j]));
                ws.push((dim_values[i][j], e_values[i][j], write_values[i][j]));
            }
            let s0_poly = self.grand_product_poly(&s0, self.table.l);
            let s_poly = self.grand_product_poly(&s, self.table.l);
            let rs_poly = self.grand_product_poly(&rs, self.logm);
            let ws_poly = self.grand_product_poly(&ws, self.logm);
            s0_polys.push(s0_poly.clone());
            s_polys.push(s_poly.clone());
            rs_polys.push(rs_poly.clone());
            ws_polys.push(ws_poly.clone());
            s0_comms.push(Pcs::commit_and_write(&self.pp, &s0_poly, &mut transcript).unwrap());
            s_comms.push(Pcs::commit_and_write(&self.pp, &s_poly, &mut transcript).unwrap());
            rs_comms.push(Pcs::commit_and_write(&self.pp, &rs_poly, &mut transcript).unwrap());
            ws_comms.push(Pcs::commit_and_write(&self.pp, &ws_poly, &mut transcript).unwrap());
        }
        let msg_4 = Message4::new(s0_comms, s_comms, rs_comms, ws_comms);

        // [round 5]
        let mut s0_sumcheck_proof = Vec::new();
        let mut s_sumcheck_proof = Vec::new();
        let mut rs_sumcheck_proof = Vec::new();
        let mut ws_sumcheck_proof = Vec::new();
        let mut r_prime = Vec::new();
        let mut r_prime2 = Vec::new();
        let mut r_prime3 = Vec::new();
        let mut r_prime4 = Vec::new();
        let mut s0_data = Vec::new();
        let mut s_data = Vec::new();
        let mut rs_data = Vec::new();
        let mut ws_data = Vec::new();
        let mut e_eval2 = Vec::new();
        let mut dim_eval = Vec::new();
        let mut read_eval = Vec::new();
        let mut final_eval = Vec::new();
        let mut e_eval2_proof = Vec::new();
        let mut dim_eval_proof = Vec::new();
        let mut read_eval_proof = Vec::new();
        let mut final_eval_proof = Vec::new();
        for i in 0..self.table.alpha {
            self.handle_grand_product_sumcheck(
                &mut s0_sumcheck_proof,
                &mut r_prime,
                &s0_polys[i],
                self.table.l,
            );
            self.handle_grand_product_sumcheck(
                &mut s_sumcheck_proof,
                &mut r_prime2,
                &s_polys[i],
                self.table.l,
            );
            self.handle_grand_product_sumcheck(
                &mut rs_sumcheck_proof,
                &mut r_prime3,
                &rs_polys[i],
                self.logm,
            );
            self.handle_grand_product_sumcheck(
                &mut ws_sumcheck_proof,
                &mut r_prime4,
                &ws_polys[i],
                self.logm,
            );
            s0_data.push(self.generate_grand_product_data(
                &s0_polys[i],
                &r_prime[i],
                &mut transcript,
            ));
            s_data.push(self.generate_grand_product_data(
                &s_polys[i],
                &r_prime2[i],
                &mut transcript,
            ));
            rs_data.push(self.generate_grand_product_data(
                &rs_polys[i],
                &r_prime3[i],
                &mut transcript,
            ));
            ws_data.push(self.generate_grand_product_data(
                &ws_polys[i],
                &r_prime4[i],
                &mut transcript,
            ));
            self.append_eval_and_proof(
                &mut e_eval2,
                &mut e_eval2_proof,
                &e_poly[i],
                &r_prime3[i],
                &mut transcript,
            );
            self.append_eval_and_proof(
                &mut dim_eval,
                &mut dim_eval_proof,
                &dim_ploy[i / self.table.k],
                &r_prime3[i],
                &mut transcript,
            );
            self.append_eval_and_proof(
                &mut read_eval,
                &mut read_eval_proof,
                &read_poly[i],
                &r_prime3[i],
                &mut transcript,
            );
            self.append_eval_and_proof(
                &mut final_eval,
                &mut final_eval_proof,
                &final_poly[i],
                &r_prime2[i],
                &mut transcript,
            );
        }

        let msg_5 = Message5::new(
            s0_sumcheck_proof,
            s_sumcheck_proof,
            rs_sumcheck_proof,
            ws_sumcheck_proof,
            r_prime,
            r_prime2,
            r_prime3,
            r_prime4,
            s0_data,
            s_data,
            rs_data,
            ws_data,
            e_eval2,
            dim_eval,
            read_eval,
            final_eval,
            e_eval2_proof,
            dim_eval_proof,
            read_eval_proof,
            final_eval_proof,
        );

        Proof {
            msg_1,
            msg_2,
            msg_3,
            msg_4,
            msg_5,
        }
    }

    fn eval_and_prove<H, S>(
        &self,
        poly: &MultilinearPolynomial<Scalar>,
        points: &Vec<Scalar>,
        transcript: &mut FiatShamirTranscript<H, S>,
    ) -> (Scalar, MvKzgProof)
    where
        S: io::Read + io::Write,
        H: Hash,
    {
        unimplemented!();
    }

    fn prove_sumcheck<H, S>(
        &self,
        poly: &MultilinearPolynomial<Scalar>,
        transcript: &mut FiatShamirTranscript<H, S>,
    ) -> (Scalar, Vec<Scalar>, BTreeMap<Query, Scalar>)
    where
        S: io::Read + io::Write,
        H: Hash,
    {
        // not sure what expression is.
        unimplemented!();
        //let expression: Expression<Scalar> = Expression::Constant(Scalar::zero());
        //let virtual_poly = VirtualPolynomial::new(&expression, vec![poly], &[], &[]);
        //let (_, rz, h_sumcheck_proof) =
        //    ClassicSumCheck::<EvaluationsProver<_>, BinaryField>::prove(
        //        &(),
        //        poly.num_vars(),
        //        virtual_poly,
        //        Scalar::from(0), // not sure what the value of sum is.
        //        transcript,
        //    )
        //    .unwrap();
        //(Scalar::zero(), rz, h_sumcheck_proof)
    }

    fn append_poly_and_values(
        &self,
        poly_list: &mut Vec<MultilinearPolynomial<Scalar>>,
        value_list: &mut Vec<Vec<Scalar>>,
        values: &[u64],
        _: usize,
    ) {
        let scalars: Vec<Scalar> = values.iter().map(|&v| Scalar::from(v)).collect();
        value_list.push(scalars.clone());
        poly_list.push(MultilinearPolynomial::new(scalars));
    }

    fn append_eval_and_proof<H, S>(
        &self,
        eval_list: &mut Vec<Scalar>,
        proof_list: &mut Vec<MvKzgProof>,
        poly: &MultilinearPolynomial<Scalar>,
        points: &Vec<Scalar>,
        transcript: &mut FiatShamirTranscript<H, S>,
    ) where
        S: io::Read + io::Write,
        H: Hash,
    {
        let (eval, proof) = self.eval_and_prove(poly, points, transcript);
        eval_list.push(eval);
        proof_list.push(proof);
    }

    fn grand_product_poly(
        &self,
        multiset: &[(Scalar, Scalar, Scalar)],
        _: usize,
    ) -> MultilinearPolynomial<Scalar> {
        let mut f: Vec<Scalar> = multiset
            .iter()
            .map(|s| hash_tuple(*s, self.tau, self.gamma))
            .collect();
        let n = f.len();
        for i in 0..n - 1 {
            f.push(f[2 * i] * f[2 * i + 1]);
        }
        f.push(Scalar::zero());
        MultilinearPolynomial::new(f)
    }

    fn grand_product_sumcheck(
        &self,
        poly: &MultilinearPolynomial<Scalar>,
        _: usize,
    ) -> (Scalar, Vec<Scalar>, BTreeMap<Query, Scalar>) {
        // Here is 100% wrong, need to know how to do the sumcheck.
        unimplemented!()
    }

    fn handle_grand_product_sumcheck(
        &self,
        data_list: &mut Vec<BTreeMap<Query, Scalar>>,
        r_list: &mut Vec<Vec<Scalar>>,
        poly: &MultilinearPolynomial<Scalar>,
        length: usize,
    ) {
        let (_, r, proof) = self.grand_product_sumcheck(poly, length);
        data_list.push(proof);
        r_list.push(r);
    }

    fn generate_grand_product_data<H, S>(
        &self,
        f: &MultilinearPolynomial<Scalar>,
        r: &Vec<Scalar>,
        transcript: &mut FiatShamirTranscript<H, S>,
    ) -> GrandProductData
    where
        S: io::Read + io::Write,
        H: Hash,
    {
        let mut arr_0_r = vec![Scalar::zero()];
        let mut arr_1_r = vec![Scalar::one()];
        let mut arr_r_0 = r.clone();
        let mut arr_r_1 = r.clone();
        let mut arr_prodcut = vec![Scalar::one(); r.len()];
        arr_0_r.extend(r);
        arr_1_r.extend(r);
        arr_r_0.push(Scalar::zero());
        arr_r_1.push(Scalar::one());
        arr_prodcut.push(Scalar::zero());

        let (f_0_r, f_0_r_proof) = self.eval_and_prove(f, &arr_0_r, transcript);
        let (f_1_r, f_1_r_proof) = self.eval_and_prove(f, &arr_1_r, transcript);
        let (f_r_0, f_r_0_proof) = self.eval_and_prove(f, &arr_r_0, transcript);
        let (f_r_1, f_r_1_proof) = self.eval_and_prove(f, &arr_r_1, transcript);
        let (product, product_proof) = self.eval_and_prove(f, &arr_prodcut, transcript);
        GrandProductData {
            f_0_r,
            f_1_r,
            f_r_0,
            f_r_1,
            product,
            f_0_r_proof,
            f_1_r_proof,
            f_r_0_proof,
            f_r_1_proof,
            product_proof,
        }
    }
}
