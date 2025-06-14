use std::{collections::{HashMap, HashSet}, iter::repeat, marker::PhantomData};

use halo2_curves::ff::{Field, PrimeField};
use itertools::Itertools;

#[cfg(feature = "parallel")]
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

use crate::{
    backend::lasso::DecomposableTable,
    pcs::{CommitmentChunk, Evaluation, PolynomialCommitmentScheme},
    piop::sum_check::{
        classic::{ClassicSumCheck, EvaluationsProver},
        SumCheck as _, VirtualPolynomial,
    },
    poly::multilinear::MultilinearPolynomial,
    util::{
        arithmetic::{fe_to_bits_le, usize_from_bits_le},
        expression::{Expression, Query, Rotation},
        transcript::TranscriptWrite,
    },
    Error,
};

use super::Poly;

type SumCheck<F> = ClassicSumCheck<EvaluationsProver<F>>;

pub struct Surge<
    F: Field + PrimeField,
    Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
> {
    lookup_indices: Vec<Vec<usize>>,
    _marker: PhantomData<F>,
    _marker2: PhantomData<Pcs>,
}

impl<
        F: Field + PrimeField,
        Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
    > Surge<F, Pcs>
{
    pub fn new() -> Self {
        Self {
            lookup_indices: vec![vec![]],
            _marker: PhantomData,
            _marker2: PhantomData,
        }
    }

    pub fn indices(&'_ self) -> Vec<&[usize]> {
        self.lookup_indices
            .iter()
            .map(|lookup_indices| lookup_indices.as_slice())
            .collect_vec()
    }

    /// computes dim_1, ..., dim_c where c == DecomposableTable::C
    pub fn commit(
        &mut self,
        table: &Box<dyn DecomposableTable<F>>,
        index_poly: &MultilinearPolynomial<F>,
    ) -> Vec<MultilinearPolynomial<F>> {
        let num_rows: usize = 1 << index_poly.num_vars();
        let num_chunks = table.chunk_bits().len();
        // get indices of non-zero columns of all rows where each index is chunked
        let indices = (0..num_rows)
            .map(|i| {
                let mut index_bits = fe_to_bits_le(index_poly[i]);
                index_bits.truncate(table.chunk_bits().iter().sum());
                assert_eq!(
                    usize_from_bits_le(&fe_to_bits_le(index_poly[i])),
                    usize_from_bits_le(&index_bits)
                );

                let mut chunked_index = repeat(0).take(num_chunks).collect_vec();
                let chunked_index_bits = table.subtable_indices(index_bits);
                chunked_index
                    .iter_mut()
                    .zip(chunked_index_bits)
                    .map(|(chunked_index, index_bits)| {
                        *chunked_index = usize_from_bits_le(&index_bits);
                    })
                    .collect_vec();
                chunked_index
            })
            .collect_vec();
        let mut dims = Vec::with_capacity(num_chunks);
        self.lookup_indices.resize(num_chunks, vec![]);
        self.lookup_indices
            .iter_mut()
            .enumerate()
            .for_each(|(i, lookup_indices)| {
                let indices = indices
                    .iter()
                    .map(|indices| {
                        lookup_indices.push(indices[i]);
                        indices[i]
                    })
                    .collect_vec();
                dims.push(MultilinearPolynomial::from_usize(indices));
            });

        dims
    }

    pub fn counter_polys(
        &self,
        table: &Box<dyn DecomposableTable<F>>,
    ) -> (Vec<MultilinearPolynomial<F>>, Vec<MultilinearPolynomial<F>>) {
        let num_chunks = table.chunk_bits().len();
        let mut read_ts_polys = Vec::with_capacity(num_chunks);
        let mut final_cts_polys = Vec::with_capacity(num_chunks);
        let chunk_bits = table.chunk_bits();
        self.lookup_indices
            .iter()
            .enumerate()
            .for_each(|(i, lookup_indices)| {
                let num_reads = lookup_indices.len();
                let memory_size = 1 << chunk_bits[i];
                let mut final_timestamps = vec![0usize; memory_size];
                let mut read_timestamps = vec![0usize; num_reads];
                (0..num_reads).for_each(|i| {
                    let memory_address = lookup_indices[i];
                    let ts = final_timestamps[memory_address];
                    read_timestamps[i] = ts;
                    let write_timestamp = ts + 1;
                    final_timestamps[memory_address] = write_timestamp;
                });
                read_ts_polys.push(MultilinearPolynomial::from_usize(read_timestamps));
                final_cts_polys.push(MultilinearPolynomial::from_usize(final_timestamps));
            });

        (read_ts_polys, final_cts_polys)
    }

    pub fn prove_sum_check(
        table: &Box<dyn DecomposableTable<F>>,
        lookup_output_poly: &Poly<F>,
        e_polys: &[&Poly<F>],
        r: &[F],
        num_vars: usize,
        points_offset: usize,
        lookup_opening_points: &mut Vec<Vec<F>>,
        lookup_opening_evals: &mut Vec<Evaluation<F>>,
        transcript: &mut impl TranscriptWrite<CommitmentChunk<F, Pcs>, F>,
    ) -> Result<(), Error> {
        let claimed_sum = Self::sum_check_claim(&r, &table, &e_polys);
        assert_eq!(claimed_sum, lookup_output_poly.evaluate(r));

        transcript.write_field_element(&claimed_sum)?;

        let expression = Self::sum_check_expression(&table);

        // proceed sumcheck
        let (final_eval, x, evals) = SumCheck::prove(
            &(),
            num_vars,
            VirtualPolynomial::new(
                &expression,
                e_polys.iter().map(|e_poly| &e_poly.poly),
                &[],
                &[r.to_vec()],
            ),
            claimed_sum,
            transcript,
        )?;

        lookup_opening_points.extend_from_slice(&[r.to_vec(), x]);
        let evals = expression
            .used_query()
            .into_iter()
            .map(|query| {
                let evaluation_value = evals.get(&query).expect("Evaluation for this query should exist in the map");
                transcript.write_field_element(evaluation_value).unwrap();
                Evaluation::new(
                    e_polys[query.poly()].offset,
                    points_offset + 1,
                    evaluation_value.clone(),
                )
            })
            .chain([Evaluation::new(
                lookup_output_poly.offset,
                points_offset,
                claimed_sum,
            )])
            .collect_vec();
        lookup_opening_evals.extend_from_slice(&evals);

        Ok(())
    }

    pub fn sum_check_claim(
        r: &[F],
        table: &Box<dyn DecomposableTable<F>>,
        e_polys: &[&Poly<F>],
    ) -> F {
        let num_memories = table.num_memories();
        assert_eq!(e_polys.len(), num_memories);
        let num_vars = e_polys[0].num_vars();
        let bh_size = 1 << num_vars;
        let eq = MultilinearPolynomial::eq_xy(r);
        // \sum_{k \in \{0, 1\}^{\log m}} (\tilde{eq}(r, k) * g(E_1(k), ..., E_{\alpha}(k)))
        
        // 根据 parallel feature 使用不同的迭代器
        #[cfg(feature = "parallel")]
        let claim = (0..bh_size)
            .into_par_iter()
            .map(|k| {
                let operands = e_polys.iter().map(|e_poly| e_poly[k]).collect_vec();
                eq[k] * table.combine_lookups(&operands)
            })
            .sum();

        #[cfg(not(feature = "parallel"))]
        let claim = (0..bh_size)
            .into_iter()
            .map(|k| {
                let operands = e_polys.iter().map(|e_poly| e_poly[k]).collect_vec();
                eq[k] * table.combine_lookups(&operands)
            })
            .sum();

        claim
    }

    // (\tilde{eq}(r, k) * g(E_1(k), ..., E_{\alpha}(k)))
    pub fn sum_check_expression(table: &Box<dyn DecomposableTable<F>>) -> Expression<F> {
        let num_memories = table.num_memories();
        let exprs = table.combine_lookup_expressions(
            (0..num_memories)
                .map(|idx| Expression::Polynomial(Query::new(idx, Rotation::cur())))
                .collect_vec(),
        );
        let eq_xy = Expression::<F>::eq_xy(0);
        eq_xy * exprs
    }
}
