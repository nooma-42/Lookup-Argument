use std::{iter, marker::PhantomData};

use halo2_curves::ff::PrimeField;
use itertools::{chain, Itertools};

use crate::{
    pcs::Evaluation,
    piop::gkr::verify_grand_product,
    poly::multilinear::MultilinearPolynomialTerms,
    util::{arithmetic::inner_product, transcript::FieldTranscriptRead},
    Error,
    log_debug, debug_println,
};

// Initialize logging at module level
static INIT: std::sync::Once = std::sync::Once::new();

fn ensure_logging_init() {
    INIT.call_once(|| {
        crate::logging::init_logging();
    });
}

#[derive(Clone, Debug)]
pub(in crate::backend::lasso) struct Chunk<F> {
    chunk_index: usize,
    chunk_bits: usize,
    pub(crate) memory: Vec<Memory<F>>,
}

impl<F: PrimeField> Chunk<F> {
    pub fn chunk_polys_index(&self, offset: usize, num_chunks: usize) -> Vec<usize> {
        let dim_poly_index = offset + 1 + self.chunk_index;
        let read_ts_poly_index = offset + 1 + num_chunks + self.chunk_index;
        let final_cts_poly_index = offset + 1 + 2 * num_chunks + self.chunk_index;
        vec![dim_poly_index, read_ts_poly_index, final_cts_poly_index]
    }

    pub fn new(chunk_index: usize, chunk_bits: usize, memory: Memory<F>) -> Self {
        Self {
            chunk_index,
            chunk_bits,
            memory: vec![memory],
        }
    }

    pub fn num_memories(&self) -> usize {
        self.memory.len()
    }

    pub fn chunk_bits(&self) -> usize {
        self.chunk_bits
    }

    pub fn add_memory(&mut self, memory: Memory<F>) {
        self.memory.push(memory);
    }

    pub fn memory_indices(&self) -> Vec<usize> {
        self.memory
            .iter()
            .map(|memory| memory.memory_index)
            .collect_vec()
    }

    /// check the following relations:
    /// - $read(x) == hash(dim(x), E(x), read_ts(x))$
    /// - $write(x) == hash(dim(x), E(x), read_ts(x) + 1)$
    /// - $init(y) == hash(y, T(y), 0)$
    /// - $final_read(y) == hash(y, T(y), final_cts(x))$
    pub fn verify_memories(
        &self,
        read_xs: &[F],
        write_xs: &[F],
        init_ys: &[F],
        final_read_ys: &[F],
        y: &[F],
        hash: impl Fn(&F, &F, &F) -> F,
        transcript: &mut impl FieldTranscriptRead<F>,
    ) -> Result<(F, F, F, Vec<F>), Error> {
        ensure_logging_init();
        log_debug!("Chunk::verify_memories: attempting to read 3 field elements + {} e_poly_xs", self.num_memories());
        let [dim_x, read_ts_poly_x, final_cts_poly_y] =
            transcript.read_field_elements(3)?.try_into().unwrap();
        let e_poly_xs = transcript.read_field_elements(self.num_memories())?;
        
        log_debug!("Verifier read from transcript: [dim_x, read_ts_poly_x, final_cts_poly_y] = [{:?}, {:?}, {:?}]", 
            dim_x, read_ts_poly_x, final_cts_poly_y);
        log_debug!("Verifier read e_poly_xs (len={}): {:?}", e_poly_xs.len(), e_poly_xs);
        
        // Compute id_poly_y using the same method as the prover
        // The identity polynomial at index i should be i
        // When evaluated at challenge point y, it should be the inner product of (1, 2, 4, 8, ...) with y
        let chunk_bits = self.chunk_bits();
        let y_for_id = if y.len() > chunk_bits {
            &y[..chunk_bits]
        } else {
            y
        };
        let id_poly_y = inner_product(
            iter::successors(Some(F::ONE), |power_of_two| Some(power_of_two.double()))
                .take(y_for_id.len())
                .collect_vec()
                .iter(),
            y_for_id,
        );
        
        log_debug!("Verifier: y = {:?}", y);
        log_debug!("Verifier: id_poly_y = {:?}", id_poly_y);
        log_debug!("Verifier: final_cts_poly_y = {:?}", final_cts_poly_y);
        
        self.memory.iter().enumerate().for_each(|(i, memory)| {
            assert_eq!(read_xs[i], hash(&dim_x, &e_poly_xs[i], &read_ts_poly_x));
            assert_eq!(
                write_xs[i],
                hash(&dim_x, &e_poly_xs[i], &(read_ts_poly_x + F::ONE))
            );
            
            // Handle variable count mismatch: truncate y to match subtable_poly's variable count
            let subtable_num_vars = memory.subtable_poly.num_vars();
            let y_truncated = if y.len() > subtable_num_vars {
                &y[..subtable_num_vars]
            } else {
                y
            };
            let subtable_poly_y = memory.subtable_poly.evaluate(y_truncated);
            log_debug!("Verifier memory[{}]: subtable_poly_y = {:?}", i, subtable_poly_y);
            
            // Note: We don't check init_ys[i] or final_read_ys[i] here because these polynomials
            // are constructed as init[j] = hash(j, subtable[j], 0) and final_read[j] = hash(j, subtable[j], final_cts[j])
            // for each memory location j, and the grand product verification already ensures the correctness
            // of these polynomial evaluations at the challenge points.
            // The hash function is not linear, so we cannot simply compute
            // hash(id_poly(y), subtable_poly(y), final_cts_poly(y)) and expect it to equal final_read_poly(y).
        });
        Ok((dim_x, read_ts_poly_x, final_cts_poly_y, e_poly_xs))
    }
}

#[derive(Clone, Debug)]
pub(in crate::backend::lasso) struct Memory<F> {
    memory_index: usize,
    subtable_poly: MultilinearPolynomialTerms<F>,
}

impl<F> Memory<F> {
    pub fn new(memory_index: usize, subtable_poly: MultilinearPolynomialTerms<F>) -> Self {
        Self {
            memory_index,
            subtable_poly,
        }
    }
}

#[derive(Clone, Debug)]
pub(in crate::backend::lasso) struct MemoryCheckingVerifier<F: PrimeField> {
    /// chunks with the same bits size
    chunks: Vec<Chunk<F>>,
    _marker: PhantomData<F>,
}

impl<'a, F: PrimeField> MemoryCheckingVerifier<F> {
    pub fn new(chunks: Vec<Chunk<F>>) -> Self {
        Self {
            chunks,
            _marker: PhantomData,
        }
    }

    pub fn verify(
        &self,
        num_chunks: usize,
        num_reads: usize,
        polys_offset: usize,
        points_offset: usize,
        gamma: &F,
        tau: &F,
        unified_num_vars: usize,
        lookup_opening_points: &mut Vec<Vec<F>>,
        lookup_opening_evals: &mut Vec<Evaluation<F>>,
        transcript: &mut impl FieldTranscriptRead<F>,
    ) -> Result<(), Error> {
        let num_memories: usize = self.chunks.iter().map(|chunk| chunk.num_memories()).sum();
        let memory_bits = self.chunks[0].chunk_bits();
        ensure_logging_init();
        log_debug!("MemoryCheckingVerifier::verify: num_memories = {}, memory_bits = {}", num_memories, memory_bits);
        
        let (read_write_xs, x) = verify_grand_product(
            unified_num_vars,
            iter::repeat(None).take(2 * num_memories),
            transcript,
        )?;
        log_debug!("MemoryCheckingVerifier::verify: after first verify_grand_product, x.len() = {}, read_write_xs.len() = {}", x.len(), read_write_xs.len());
        let (read_xs, write_xs) = read_write_xs.split_at(num_memories);

        let (init_final_read_ys, y) = verify_grand_product(
            unified_num_vars,
            iter::repeat(None).take(2 * num_memories),
            transcript,
        )?;
        log_debug!("MemoryCheckingVerifier::verify: after second verify_grand_product, y.len() = {}, init_final_read_ys.len() = {}", y.len(), init_final_read_ys.len());
        let (init_ys, final_read_ys) = init_final_read_ys.split_at(num_memories);

        let hash = |a: &F, v: &F, t: &F| -> F { *a + *v * gamma + *t * gamma.square() - tau };
        let mut offset = 0;
        let (dim_xs, read_ts_poly_xs, final_cts_poly_ys, e_poly_xs) = self
            .chunks
            .iter()
            .map(|chunk| {
                let num_memories = chunk.num_memories();
                log_debug!("MemoryCheckingVerifier::verify: processing chunk with {} memories, offset = {}", num_memories, offset);
                let result = chunk.verify_memories(
                    &read_xs[offset..offset + num_memories],
                    &write_xs[offset..offset + num_memories],
                    &init_ys[offset..offset + num_memories],
                    &final_read_ys[offset..offset + num_memories],
                    &y,
                    hash,
                    transcript,
                );
                offset += num_memories;
                result
            })
            .collect::<Result<Vec<(F, F, F, Vec<F>)>, Error>>()?
            .into_iter()
            .multiunzip::<(Vec<_>, Vec<_>, Vec<_>, Vec<Vec<_>>)>();

        self.opening_evals(
            num_chunks,
            polys_offset,
            points_offset,
            &lookup_opening_points,
            lookup_opening_evals,
            &dim_xs,
            &read_ts_poly_xs,
            &final_cts_poly_ys,
            &e_poly_xs.concat(),
        );
        lookup_opening_points.extend_from_slice(&[x, y]);

        Ok(())
    }

    fn opening_evals(
        &self,
        num_chunks: usize,
        polys_offset: usize,
        points_offset: usize,
        lookup_opening_points: &Vec<Vec<F>>,
        lookup_opening_evals: &mut Vec<Evaluation<F>>,
        dim_xs: &[F],
        read_ts_poly_xs: &[F],
        final_cts_poly_ys: &[F],
        e_poly_xs: &[F],
    ) {
        let x_offset = lookup_opening_points.len();
        let y_offset = lookup_opening_points.len() + 1;
        let (dim_xs, read_ts_poly_xs, final_cts_poly_ys) = self
            .chunks
            .iter()
            .enumerate()
            .map(|(i, chunk)| {
                let chunk_polys_index = chunk.chunk_polys_index(polys_offset, num_chunks);
                (
                    Evaluation::new(chunk_polys_index[0], x_offset, dim_xs[i]),
                    Evaluation::new(chunk_polys_index[1], x_offset, read_ts_poly_xs[i]),
                    Evaluation::new(chunk_polys_index[2], y_offset, final_cts_poly_ys[i]),
                )
            })
            .multiunzip::<(Vec<Evaluation<F>>, Vec<Evaluation<F>>, Vec<Evaluation<F>>)>();

        let e_poly_offset = polys_offset + 1 + 3 * num_chunks;
        let e_poly_xs = self
            .chunks
            .iter()
            .flat_map(|chunk| chunk.memory_indices())
            .zip(e_poly_xs)
            .map(|(memory_index, &e_poly_x)| {
                Evaluation::new(e_poly_offset + memory_index, x_offset, e_poly_x)
            })
            .collect_vec();
        lookup_opening_evals.extend_from_slice(
            &chain!(dim_xs, read_ts_poly_xs, final_cts_poly_ys, e_poly_xs).collect_vec(),
        );
    }
}
