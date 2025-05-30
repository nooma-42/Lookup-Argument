use std::{
    collections::{HashMap, HashSet},
    marker::PhantomData,
};

use halo2_curves::ff::{Field, PrimeField};
use itertools::{chain, Itertools};

use crate::{
    pcs::{CommitmentChunk, Evaluation, PolynomialCommitmentScheme},
    poly::multilinear::MultilinearPolynomial,
    util::{
        expression::{rotate::{Rotatable, BinaryField}, CommonPolynomial, Expression},
        impl_index,
        parallel::parallelize,
        transcript::{FieldTranscriptWrite, TranscriptWrite},
    },
    Error,
};

use super::{memory_checking::MemoryCheckingProver, DecomposableTable};

mod surge;

pub use surge::Surge;

#[derive(Clone, Debug)]
pub struct Poly<F> {
    /// polynomial offset in batch opening
    pub(crate) offset: usize,
    pub(crate) poly: MultilinearPolynomial<F>,
}

impl_index!(Poly, poly);

impl<F: PrimeField> Poly<F> {
    pub fn num_vars(&self) -> usize {
        self.poly.num_vars()
    }

    pub fn evaluate(&self, x: &[F]) -> F {
        self.poly.evaluate(x)
    }
}

#[derive(Clone, Debug)]
pub struct Chunk<'a, F: PrimeField> {
    pub(super) chunk_index: usize,
    pub(super) chunk_bits: usize,
    pub(super) dim: &'a Poly<F>,
    pub(super) read_ts_poly: &'a Poly<F>,
    pub(super) final_cts_poly: &'a Poly<F>,
    pub(super) memories: Vec<Memory<'a, F>>,
}

impl<'a, F: PrimeField> Chunk<'a, F> {
    fn new(
        chunk_index: usize,
        chunk_bits: usize,
        dim: &'a Poly<F>,
        read_ts_poly: &'a Poly<F>,
        final_cts_poly: &'a Poly<F>,
        memory: Memory<'a, F>,
    ) -> Self {
        // sanity check
        assert_eq!(dim.num_vars(), read_ts_poly.num_vars());

        Self {
            chunk_index,
            chunk_bits,
            dim,
            read_ts_poly,
            final_cts_poly,
            memories: vec![memory],
        }
    }

    pub fn chunk_index(&self) -> usize {
        self.chunk_index
    }

    pub fn chunk_bits(&self) -> usize {
        self.chunk_bits
    }

    pub fn num_reads(&self) -> usize {
        1 << self.dim.num_vars()
    }

    pub fn chunk_polys(&self) -> impl Iterator<Item = &'a MultilinearPolynomial<F>> {
        chain!([
            &self.dim.poly,
            &self.read_ts_poly.poly,
            &self.final_cts_poly.poly
        ])
    }

    pub fn chunk_poly_evals(&self, x: &[F], y: &[F]) -> Vec<F> {
        // Pad y to match the extended polynomial's variable count if needed
        let final_cts_num_vars = self.final_cts_poly.num_vars();
        let y_padded = if y.len() < final_cts_num_vars {
            let mut padded = y.to_vec();
            padded.resize(final_cts_num_vars, F::ZERO);
            padded
        } else {
            y.to_vec()
        };
        
        vec![
            self.dim.evaluate(x),
            self.read_ts_poly.evaluate(x),
            self.final_cts_poly.evaluate(&y_padded),
        ]
    }

    pub fn e_poly_evals(&self, x: &[F]) -> Vec<F> {
        self.memories
            .iter()
            .map(|memory| memory.e_poly.evaluate(x))
            .collect_vec()
    }

    pub(super) fn memories(&self) -> impl Iterator<Item = &Memory<F>> {
        self.memories.iter()
    }

    pub(super) fn add_memory(&mut self, memory: Memory<'a, F>) {
        // sanity check
        let chunk_bits = self.chunk_bits();
        let num_reads = self.num_reads();
        assert_eq!(chunk_bits, memory.subtable_poly.num_vars());
        assert_eq!(num_reads, memory.e_poly.num_vars());

        self.memories.push(memory);
    }
}

#[derive(Clone, Debug)]
pub(super) struct Memory<'a, F: PrimeField> {
    subtable_poly: &'a MultilinearPolynomial<F>,
    pub(crate) e_poly: &'a Poly<F>,
}

impl<'a, F: PrimeField> Memory<'a, F> {
    fn new(subtable_poly: &'a MultilinearPolynomial<F>, e_poly: &'a Poly<F>) -> Self {
        Self {
            subtable_poly,
            e_poly,
        }
    }

    pub fn polys(&'a self) -> impl Iterator<Item = &'a MultilinearPolynomial<F>> {
        chain!([&self.subtable_poly, &self.e_poly.poly])
    }
}

pub struct LassoProver<
    F: Field + PrimeField,
    Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
>(PhantomData<F>, PhantomData<Pcs>);

impl<
        F: Field + PrimeField,
        Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
    > LassoProver<F, Pcs>
{
    pub fn lookup_poly(
        lookup: &(&Expression<F>, &Expression<F>),
        polys: &[&MultilinearPolynomial<F>],
    ) -> (MultilinearPolynomial<F>, MultilinearPolynomial<F>) {
        let num_vars = polys[0].num_vars();
        let expression = lookup.0 + lookup.1;
        let lagranges = {
            let bh = BinaryField::new(num_vars).iter().collect_vec();
            expression
                .used_langrange()
                .into_iter()
                .map(|i| (i, bh[i.rem_euclid(1 << num_vars) as usize]))
                .collect::<HashSet<_>>()
        };
        let bh = BinaryField::new(num_vars);

        let evaluate = |expression: &Expression<F>| {
            let mut evals = vec![F::ZERO; 1 << num_vars];
            parallelize(&mut evals, |(evals, start)| {
                for (b, eval) in (start..).zip(evals) {
                    *eval = expression.evaluate(
                        &|constant| constant,
                        &|common_poly| match common_poly {
                            CommonPolynomial::Identity => F::from(b as u64),
                            CommonPolynomial::Lagrange(i) => {
                                if lagranges.contains(&(i, b)) {
                                    F::ONE
                                } else {
                                    F::ZERO
                                }
                            }
                            CommonPolynomial::EqXY(_) => unreachable!(),
                        },
                        &|query| polys[query.poly()][bh.rotate(b, query.rotation())],
                        &|_| unreachable!(),
                        &|value| -value,
                        &|lhs, rhs| lhs + &rhs,
                        &|lhs, rhs| lhs * &rhs,
                        &|value, scalar| value * &scalar,
                    );
                }
            });
            MultilinearPolynomial::new(evals)
        };

        let (input, index) = lookup;
        (evaluate(input), evaluate(index))
    }
}

impl<
        F: Field + PrimeField,
        Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
    > LassoProver<F, Pcs>
{
    fn e_polys(
        table: &Box<dyn DecomposableTable<F>>,
        subtable_polys: &[&MultilinearPolynomial<F>],
        indices: &Vec<&[usize]>,
    ) -> Vec<MultilinearPolynomial<F>> {
        let num_chunks = table.chunk_bits().len();
        let num_memories = table.num_memories();
        assert_eq!(indices.len(), num_chunks);
        let num_reads = indices[0].len();
        (0..num_memories)
            .map(|i| {
                let mut e_poly = Vec::with_capacity(num_reads);
                let subtable_poly = subtable_polys[table.memory_to_subtable_index(i)];
                let index = indices[table.memory_to_chunk_index(i)];
                (0..num_reads).for_each(|j| {
                    e_poly.push(subtable_poly[index[j]]);
                });
                MultilinearPolynomial::new(e_poly)
            })
            .collect_vec()
    }

    fn chunks<'a>(
        table: &Box<dyn DecomposableTable<F>>,
        subtable_polys: &'a [&MultilinearPolynomial<F>],
        dims: &'a [Poly<F>],
        read_ts_polys: &'a [Poly<F>],
        final_cts_polys: &'a [Poly<F>],
        e_polys: &'a [Poly<F>],
    ) -> Vec<Chunk<'a, F>> {
        // key: chunk index, value: chunk
        let mut chunk_map: HashMap<usize, Chunk<F>> = HashMap::new();

        let num_memories = table.num_memories();
        let memories = (0..num_memories).map(|memory_index| {
            let subtable_poly = subtable_polys[table.memory_to_subtable_index(memory_index)];
            Memory::new(subtable_poly, &e_polys[memory_index])
        });
        memories.enumerate().for_each(|(memory_index, memory)| {
            let chunk_index = table.memory_to_chunk_index(memory_index);
            if chunk_map.get(&chunk_index).is_some() {
                chunk_map.entry(chunk_index).and_modify(|chunk| {
                    chunk.add_memory(memory);
                });
            } else {
                let dim = &dims[chunk_index];
                let read_ts_poly = &read_ts_polys[chunk_index];
                let final_cts_poly = &final_cts_polys[chunk_index];
                let chunk_bits = table.chunk_bits()[chunk_index];
                chunk_map.insert(
                    chunk_index,
                    Chunk::new(chunk_index, chunk_bits, dim, read_ts_poly, final_cts_poly, memory),
                );
            }
        });

        // sanity check
        {
            let num_chunks = table.chunk_bits().len();
            assert_eq!(chunk_map.len(), num_chunks);
        }

        let mut chunks = chunk_map.into_iter().collect_vec();
        chunks.sort_by_key(|(chunk_index, _)| *chunk_index);
        chunks.into_iter().map(|(_, chunk)| chunk).collect_vec()
    }

    pub fn prove_sum_check(
        points_offset: usize,
        lookup_opening_points: &mut Vec<Vec<F>>,
        lookup_opening_evals: &mut Vec<Evaluation<F>>,
        table: &Box<dyn DecomposableTable<F>>,
        lookup_output_poly: &Poly<F>,
        e_polys: &[&Poly<F>],
        r: &[F],
        num_vars: usize,
        transcript: &mut impl TranscriptWrite<CommitmentChunk<F, Pcs>, F>,
    ) -> Result<(), Error> {
        Surge::<F, Pcs>::prove_sum_check(
            table,
            lookup_output_poly,
            &e_polys,
            r,
            num_vars,
            points_offset,
            lookup_opening_points,
            lookup_opening_evals,
            transcript,
        )
    }

    fn prepare_memory_checking<'a>(
        points_offset: usize,
        table: &Box<dyn DecomposableTable<F>>,
        subtable_polys: &'a [&MultilinearPolynomial<F>],
        dims: &'a [Poly<F>],
        read_ts_polys: &'a [Poly<F>],
        final_cts_polys: &'a [Poly<F>],
        e_polys: &'a [Poly<F>],
        gamma: &F,
        tau: &F,
        unified_num_vars: usize,
    ) -> Vec<MemoryCheckingProver<'a, F>> {
        let chunks = Self::chunks(
            table,
            subtable_polys,
            dims,
            read_ts_polys,
            final_cts_polys,
            e_polys,
        );
        let chunk_bits = table.chunk_bits();
        // key: chunk bits, value: chunks
        let mut chunk_map: HashMap<usize, Vec<Chunk<F>>> = HashMap::new();

        chunks.iter().enumerate().for_each(|(chunk_index, chunk)| {
            let chunk_bits = chunk_bits[chunk_index];
            if let Some(_) = chunk_map.get(&chunk_bits) {
                chunk_map.entry(chunk_bits).and_modify(|chunks| {
                    chunks.push(chunk.clone());
                });
            } else {
                chunk_map.insert(chunk_bits, vec![chunk.clone()]);
            }
        });

        chunk_map
            .into_iter()
            .enumerate()
            .map(|(index, (_, chunks))| {
                let points_offset = points_offset + 2 + 2 * index;
                MemoryCheckingProver::new(points_offset, chunks, tau, gamma, unified_num_vars)
            })
            .collect_vec()
    }

    pub fn memory_checking<'a>(
        points_offset: usize,
        lookup_opening_points: &mut Vec<Vec<F>>,
        lookup_opening_evals: &mut Vec<Evaluation<F>>,
        table: &Box<dyn DecomposableTable<F>>,
        subtable_polys: &'a [&MultilinearPolynomial<F>],
        dims: &'a [Poly<F>],
        read_ts_polys: &'a [Poly<F>],
        final_cts_polys: &'a [Poly<F>],
        e_polys: &'a [Poly<F>],
        gamma: &F,
        tau: &F,
        unified_num_vars: usize,
        transcript: &mut impl FieldTranscriptWrite<F>,
    ) -> Result<(), Error> {
        let mut memory_checking = LassoProver::<F, Pcs>::prepare_memory_checking(
            points_offset,
            &table,
            &subtable_polys,
            &dims,
            &read_ts_polys,
            &final_cts_polys,
            &e_polys,
            &gamma,
            &tau,
            unified_num_vars,
        );

        memory_checking
            .iter_mut()
            .map(|memory_checking| {
                memory_checking.prove(
                    points_offset,
                    lookup_opening_points,
                    lookup_opening_evals,
                    transcript,
                )
            })
            .collect::<Result<Vec<()>, Error>>()?;
        Ok(())
    }

    pub fn commit(
        pp: &Pcs::ProverParam,
        lookup_polys_offset: usize,
        table: &Box<dyn DecomposableTable<F>>,
        subtable_polys: &[&MultilinearPolynomial<F>],
        lookup_output_poly: MultilinearPolynomial<F>,
        lookup_index_poly: &MultilinearPolynomial<F>,
        transcript: &mut impl TranscriptWrite<Pcs::CommitmentChunk, F>,
    ) -> Result<(Vec<Vec<Poly<F>>>, Vec<Vec<Pcs::Commitment>>), Error> {
        let num_chunks = table.chunk_bits().len();

        // Calculate unified_num_vars for memory checking consistency first
        let chunk_bits = table.chunk_bits();
        let original_num_vars_for_lookups = lookup_index_poly.num_vars();
        let unified_num_vars = *chunk_bits.iter().max().unwrap().max(&original_num_vars_for_lookups);

        // commit to lookup_output_poly
        let lookup_output_comm = Pcs::commit_and_write(&pp, &lookup_output_poly, transcript)?;

        // get surge and dims
        let mut surge = Surge::<F, Pcs>::new();

        // get dims first, then extend and commit
        let dims = surge.commit(&table, lookup_index_poly);

        // Helper function to extend a polynomial to unified_num_vars
        let extend_poly_to_unified_vars = |poly: MultilinearPolynomial<F>| -> MultilinearPolynomial<F> {
            let current_vars = poly.num_vars();
            if current_vars >= unified_num_vars {
                poly
            } else {
                let current_size = 1 << current_vars;
                let unified_size = 1 << unified_num_vars;
                let mut extended_evals = vec![F::ZERO; unified_size];
                
                // Copy the original evaluations and pad with zeros
                for i in 0..current_size {
                    extended_evals[i] = poly[i];
                }
                
                MultilinearPolynomial::new(extended_evals)
            }
        };

        // get e_polys & read_ts_polys & final_cts_polys
        let e_polys = {
            let indices = surge.indices();
            LassoProver::<F, Pcs>::e_polys(&table, subtable_polys, &indices)
        };
        let (read_ts_polys, final_cts_polys) = surge.counter_polys(&table);

        // Extend dims to unified_num_vars and commit
        let extended_dims = dims.into_iter().map(extend_poly_to_unified_vars).collect_vec();
        let dim_comms = Pcs::batch_commit_and_write(pp, &extended_dims, transcript)?;
        
        // Extend read_ts_polys, final_cts_polys, and e_polys to unified_num_vars before committing
        let extended_read_ts_polys = read_ts_polys.into_iter().map(extend_poly_to_unified_vars).collect_vec();
        let extended_final_cts_polys = final_cts_polys.into_iter().map(extend_poly_to_unified_vars).collect_vec();
        let extended_e_polys = e_polys.into_iter().map(extend_poly_to_unified_vars).collect_vec();

        // commit to read_ts_polys & final_cts_polys & e_polys (now extended)
        let read_ts_comms = Pcs::batch_commit_and_write(&pp, &extended_read_ts_polys, transcript)?;
        let final_cts_comms = Pcs::batch_commit_and_write(&pp, &extended_final_cts_polys, transcript)?;
        let e_comms = Pcs::batch_commit_and_write(&pp, extended_e_polys.as_slice(), transcript)?;

        let lookup_output_poly = Poly {
            offset: lookup_polys_offset,
            poly: lookup_output_poly,
        };

        let dims = extended_dims
            .into_iter()
            .enumerate()
            .map(|(chunk_index, dim)| Poly {
                offset: lookup_polys_offset + 1 + chunk_index,
                poly: dim,
            })
            .collect_vec();

        let read_ts_polys = extended_read_ts_polys
            .into_iter()
            .enumerate()
            .map(|(chunk_index, read_ts_poly)| Poly {
                offset: lookup_polys_offset + 1 + num_chunks + chunk_index,
                poly: read_ts_poly,
            })
            .collect_vec();

        let final_cts_polys = extended_final_cts_polys
            .into_iter()
            .enumerate()
            .map(|(chunk_index, final_cts_poly)| Poly {
                offset: lookup_polys_offset + 1 + 2 * num_chunks + chunk_index,
                poly: final_cts_poly,
            })
            .collect_vec();

        let e_polys = extended_e_polys
            .into_iter()
            .enumerate()
            .map(|(memory_index, e_poly)| Poly {
                offset: lookup_polys_offset + 1 + 3 * num_chunks + memory_index,
                poly: e_poly,
            })
            .collect_vec();

        Ok((
            vec![
                vec![lookup_output_poly],
                dims,
                read_ts_polys,
                final_cts_polys,
                e_polys,
            ],
            vec![
                vec![lookup_output_comm],
                dim_comms,
                read_ts_comms,
                final_cts_comms,
                e_comms,
            ],
        ))
    }
}
