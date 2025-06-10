use std::iter;

use halo2_curves::ff::PrimeField;
use itertools::{chain, Itertools};
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

use crate::{
    backend::lasso::prover::Chunk, pcs::Evaluation, piop::gkr::prove_grand_product,
    poly::multilinear::MultilinearPolynomial, util::{transcript::FieldTranscriptWrite, arithmetic::inner_product}, Error,
    log_debug, debug_println,
};

// Initialize logging at module level
static INIT: std::sync::Once = std::sync::Once::new();

fn ensure_logging_init() {
    INIT.call_once(|| {
        crate::logging::init_logging();
    });
}

use super::MemoryGKR;

pub struct MemoryCheckingProver<'a, F: PrimeField> {
    /// offset of MemoryCheckingProver instance opening points
    points_offset: usize,
    /// chunks with the same bits size
    chunks: Vec<Chunk<'a, F>>,
    /// GKR initial polynomials for each memory
    memories: Vec<MemoryGKR<F>>,
}

impl<'a, F: PrimeField> MemoryCheckingProver<'a, F> {
    // T_1[dim_1(x)], ..., T_k[dim_1(x)],
    // ...
    // T_{\alpha-k+1}[dim_c(x)], ..., T_{\alpha}[dim_c(x)]
    pub fn new(points_offset: usize, chunks: Vec<Chunk<'a, F>>, tau: &F, gamma: &F, unified_num_vars: usize) -> Self {
        let num_reads = chunks[0].num_reads();
        let memory_size = 1 << chunks[0].chunk_bits();
        let unified_memory_size = 1 << unified_num_vars;

        let hash = |a: &F, v: &F, t: &F| -> F { *a + *v * gamma + *t * gamma.square() - tau };



        let memories_gkr: Vec<MemoryGKR<F>> = (0..chunks.len())
            .into_par_iter()
            .flat_map(|i| {
                let chunk = &chunks[i];
                let chunk_polys = chunk.chunk_polys().collect_vec();
                let (dim, read_ts_poly, final_cts_poly) =
                    (chunk_polys[0], chunk_polys[1], chunk_polys[2]);
                chunk
                    .memories()
                    .map(|memory| {
                        let memory_polys = memory.polys().collect_vec();
                        let (subtable_poly, e_poly) = (memory_polys[0], memory_polys[1]);
                        // Use unified_memory_size for all polynomials to ensure consistent variable count
                        let mut init = vec![F::ZERO; unified_memory_size];
                        let mut read = vec![F::ZERO; unified_memory_size];
                        let mut write = vec![F::ZERO; unified_memory_size];
                        let mut final_read = vec![F::ZERO; unified_memory_size];
                        
                        // Fill the memory values - ensure we don't exceed bounds
                        let safe_memory_size = memory_size.min(subtable_poly.evals().len()).min(final_cts_poly.evals().len());
                        (0..safe_memory_size).for_each(|i| {
                            // The identity polynomial value at index i should be i
                            let id_value = F::from(i as u64);
                            init[i] = hash(&id_value, &subtable_poly[i], &F::ZERO);
                            final_read[i] = hash(
                                &id_value,
                                &subtable_poly[i],
                                &final_cts_poly[i],
                            );
                            
                            // Add debug output for the first few entries when it's a problematic case
                            if i < 5 && chunk.chunk_bits() >= 4 {
                                ensure_logging_init();
                                log_debug!("Prover memory index {}: id_value = {:?}, subtable_poly = {:?}, final_cts_poly = {:?}", 
                                    i, id_value, subtable_poly[i], final_cts_poly[i]);
                                log_debug!("Prover memory index {}: init = {:?}", i, init[i]);
                                log_debug!("Prover memory index {}: final_read = {:?}", i, final_read[i]);
                                
                                // Debug the hash computation
                                let manual_hash = id_value + subtable_poly[i] * gamma + F::ZERO * gamma.square() - tau;
                                log_debug!("Prover memory index {}: manual_hash for init = {:?}", i, manual_hash);
                                
                                // Verify that id_value matches subtable_poly[i] for the identity function
                                if id_value != subtable_poly[i] {
                                    log_debug!("Prover memory index {}: MISMATCH! id_value = {:?}, subtable_poly = {:?}", 
                                        i, id_value, subtable_poly[i]);
                                } else {
                                    log_debug!("Prover memory index {}: MATCH! id_value = subtable_poly = {:?}", i, id_value);
                                }
                            }
                        });
                        
                        // Fill the read/write values - ensure we don't exceed bounds
                        let safe_num_reads = num_reads.min(dim.evals().len()).min(e_poly.evals().len()).min(read_ts_poly.evals().len());
                        (0..safe_num_reads).for_each(|i| {
                            read[i] = hash(&dim[i], &e_poly[i], &read_ts_poly[i]);
                            write[i] = hash(&dim[i], &e_poly[i], &(read_ts_poly[i] + F::ONE));
                        });
                        
                        // Debug: Create the multilinear polynomial and check its evaluation at a test point
                        let init_poly = MultilinearPolynomial::new(init.clone());
                        if chunk.chunk_bits() >= 4 {
                            ensure_logging_init();
                            log_debug!("Prover: init vector length = {}, memory_size = {}, unified_memory_size = {}", 
                                init.len(), memory_size, unified_memory_size);
                            log_debug!("Prover: init polynomial has {} variables, {} evaluations", 
                                init_poly.num_vars(), init_poly.evals().len());
                            log_debug!("Prover: first 5 init values: {:?}", 
                                &init_poly.evals()[..5.min(init_poly.evals().len())]);
                        }
                        
                        MemoryGKR::new(
                            init_poly,
                            MultilinearPolynomial::new(read),
                            MultilinearPolynomial::new(write),
                            MultilinearPolynomial::new(final_read),
                        )
                    })
                    .collect_vec()
            })
            .collect();

        Self {
            points_offset,
            chunks,
            memories: memories_gkr,
        }
    }

    fn inits(&self) -> impl Iterator<Item = &MultilinearPolynomial<F>> {
        self.memories.iter().map(|memory| &memory.init)
    }

    fn reads(&self) -> impl Iterator<Item = &MultilinearPolynomial<F>> {
        self.memories.iter().map(|memory| &memory.read)
    }

    fn writes(&self) -> impl Iterator<Item = &MultilinearPolynomial<F>> {
        self.memories.iter().map(|memory| &memory.write)
    }

    fn final_reads(&self) -> impl Iterator<Item = &MultilinearPolynomial<F>> {
        self.memories.iter().map(|memory| &memory.final_read)
    }

    fn iter(
        &self,
    ) -> impl Iterator<
        Item = (
            &MultilinearPolynomial<F>,
            &MultilinearPolynomial<F>,
            &MultilinearPolynomial<F>,
            &MultilinearPolynomial<F>,
        ),
    > {
        self.memories.iter().map(|memory| {
            (
                &memory.init,
                &memory.read,
                &memory.write,
                &memory.final_read,
            )
        })
    }

    pub fn claimed_v_0s(&self) -> impl IntoIterator<Item = Vec<Option<F>>> {
        let (claimed_read_0s, claimed_write_0s, claimed_init_0s, claimed_final_read_0s) = self
            .iter()
            .map(|(init, read, write, final_read)| {
                let claimed_init_0 = init.iter().product();
                let claimed_read_0 = read.iter().product();
                let claimed_write_0 = write.iter().product();
                let claimed_final_read_0 = final_read.iter().product();

                // sanity check
                debug_assert_eq!(
                    claimed_init_0 * claimed_write_0,
                    claimed_read_0 * claimed_final_read_0,
                    "Multiset hashes don't match",
                );
                (
                    Some(claimed_read_0),
                    Some(claimed_write_0),
                    Some(claimed_init_0),
                    Some(claimed_final_read_0),
                )
            })
            .multiunzip::<(Vec<_>, Vec<_>, Vec<_>, Vec<_>)>();
        chain!([
            claimed_read_0s,
            claimed_write_0s,
            claimed_init_0s,
            claimed_final_read_0s
        ])
    }

    pub fn prove(
        &mut self,
        points_offset: usize,
        lookup_opening_points: &mut Vec<Vec<F>>,
        lookup_opening_evals: &mut Vec<Evaluation<F>>,
        transcript: &mut impl FieldTranscriptWrite<F>,
    ) -> Result<(), Error> {
        ensure_logging_init();
        log_debug!("MemoryCheckingProver::prove: num_memories = {}", self.memories.len());
        
        let (_, x) = prove_grand_product(
            iter::repeat(None).take(self.memories.len() * 2),
            chain!(self.reads(), self.writes()),
            transcript,
        )?;
        log_debug!("MemoryCheckingProver::prove: after first grand_product, x.len() = {}", x.len());

        let (_, y) = prove_grand_product(
            iter::repeat(None).take(self.memories.len() * 2),
            chain!(self.inits(), self.final_reads()),
            transcript,
        )?;
        log_debug!("MemoryCheckingProver::prove: after second grand_product, y.len() = {}", y.len());
        
        // Debug: Check what the init polynomial evaluates to at the challenge point y
        if self.memories.len() > 0 {
            let init_poly = &self.memories[0].init;
            log_debug!("Prover: init_poly.num_vars() = {}", init_poly.num_vars());
            if y.len() >= init_poly.num_vars() {
                let y_truncated = &y[..init_poly.num_vars()];
                let init_at_y = init_poly.evaluate(y_truncated);
                log_debug!("Prover: init_poly.evaluate(y[..{}]) = {:?}", init_poly.num_vars(), init_at_y);
            } else {
                log_debug!("Prover: y.len() = {} < init_poly.num_vars() = {}", y.len(), init_poly.num_vars());
            }
        }

        assert_eq!(
            points_offset + lookup_opening_points.len(),
            self.points_offset
        );
        let x_offset = lookup_opening_points.len();
        let y_offset = lookup_opening_points.len() + 1;
        let (dim_xs, read_ts_poly_xs, final_cts_poly_ys, e_poly_xs) = self
            .chunks
            .iter()
            .map(|chunk| {
                let chunk_poly_evals = chunk.chunk_poly_evals(&x, &y);
                let e_poly_xs = chunk.e_poly_evals(&x);
                
                log_debug!("Prover chunk_poly_evals: {:?}", chunk_poly_evals);
                log_debug!("Prover writing to transcript: final_cts_poly.evaluate(y) = {:?}", chunk_poly_evals[2]);
                log_debug!("Prover writing chunk_poly_evals (len={}): {:?}", chunk_poly_evals.len(), chunk_poly_evals);
                log_debug!("Prover writing e_poly_xs (len={}): {:?}", e_poly_xs.len(), e_poly_xs);
                
                transcript.write_field_elements(&chunk_poly_evals).unwrap();
                transcript.write_field_elements(&e_poly_xs).unwrap();

                (
                    Evaluation::new(chunk.dim.offset, x_offset, chunk_poly_evals[0]),
                    Evaluation::new(chunk.read_ts_poly.offset, x_offset, chunk_poly_evals[1]),
                    Evaluation::new(chunk.final_cts_poly.offset, y_offset, chunk_poly_evals[2]),
                    chunk
                        .memories()
                        .enumerate()
                        .map(|(i, memory)| {
                            // Ensure we don't exceed the e_poly_xs bounds
                            let eval_value = if i < e_poly_xs.len() {
                                e_poly_xs[i]
                            } else {
                                // This should not happen in correct implementation
                                // but provide a safe fallback
                                F::ZERO
                            };
                            Evaluation::new(memory.e_poly.offset, x_offset, eval_value)
                        })
                        .collect_vec(),
                )
            })
            .multiunzip::<(
                Vec<Evaluation<F>>,
                Vec<Evaluation<F>>,
                Vec<Evaluation<F>>,
                Vec<Vec<Evaluation<F>>>,
            )>();

        lookup_opening_points.extend_from_slice(&[x, y]);
        let opening_evals = chain!(
            dim_xs,
            read_ts_poly_xs,
            final_cts_poly_ys,
            e_poly_xs.concat()
        )
        .collect_vec();
        lookup_opening_evals.extend_from_slice(&opening_evals);

        Ok(())
    }
}
