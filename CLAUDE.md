# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Rust workspace for benchmarking lookup argument implementations in zero-knowledge proof systems. The project contains comprehensive benchmarks for various proving systems including CQ, Baloo, LogupGKR, Plookup, Caulk, and Lasso with parallel execution capabilities.

## Build Commands

### Basic Build
```bash
cargo build
```

### Running Tests
```bash
cargo test
```

### Running Benchmarks
```bash
# Basic benchmark - single system
cargo bench --bench proof_system -- --system lasso --k 8 --format table

# Comprehensive benchmark with multiple systems and parameters
cargo bench --bench proof_system -- --system all --k 8..12 --ratio 2,4,8 --format table

# Non-parallel benchmark (for comparison)
cargo bench --bench proof_system --no-default-features --features timer -- --system lasso --k 8..10
```

### Running with Debug Logging
```bash
# Set log level via environment variable
LOOKUP_LOG_LEVEL=debug cargo bench --bench proof_system -- --system lasso --k 8
LOOKUP_LOG_LEVEL=info cargo test  # Default level
LOOKUP_LOG_LEVEL=silent cargo test  # No debug output
```

### Lint and Format
No specific lint/format commands are configured. Use standard Rust tooling:
```bash
cargo fmt
cargo clippy
```

## Workspace Structure

This is a Rust workspace with two main crates:
- **`plonkish_backend/`**: Core proving system implementations and utilities
- **`benchmark/`**: Benchmarking framework and test harnesses

### Core Architecture

The codebase implements multiple lookup argument protocols:

#### Proving Systems (in `plonkish_backend/src/backend/`)
- **CQ**: Constant-time lookup arguments
- **Baloo**: Sublinear proving time for large tables  
- **LogupGKR**: GKR-based lookup arguments with parallel/non-parallel modes
- **Plookup**: Original Plookup protocol
- **Caulk**: Vector commitment-based lookups
- **Lasso**: Memory-checking based lookups with parallel/non-parallel modes

#### Key Components
- **PCS Layer** (`pcs/`): Polynomial commitment schemes (KZG, IPA, Gemini, etc.)
- **PIOP Layer** (`piop/`): Interactive Oracle Proofs (sum-check, GKR protocols)
- **Frontend** (`frontend/`): Circuit interfaces (Halo2 integration)
- **Utilities** (`util/`): Arithmetic, parallel execution, logging, transcripts

### Parallelization Architecture

The project supports two levels of parallelization:
1. **Benchmark-level**: Multiple benchmark tasks run concurrently using `rayon`
2. **Algorithm-level**: Individual proving systems can run in parallel/non-parallel mode based on feature flags

Systems with parallel/non-parallel support:
- LogupGKR: Full parallel/non-parallel support
- Lasso: Full parallel/non-parallel support  
- HyperPlonk: Always parallel

## Feature Flags

### Important Features
- `parallel`: Enables parallel execution within algorithms (default)
- `timer`: Enables performance timing with `ark-std` 
- `frontend-halo2`: Enables Halo2 circuit integration (default)
- `benchmark`: Required for running benchmarks

### Feature Combinations
```bash
# Default build (parallel + frontend-halo2)
cargo build

# Non-parallel build for comparison
cargo build --no-default-features --features timer

# Benchmark-specific build
cargo build --features benchmark
```

## Benchmark Parameters

### Key Parameters
- `--system`: Choose systems (all, cq, baloo, logupgkr, plookup, caulk, lasso)
- `--k`: Log table size, supports ranges like `8..12` or single values `10`
- `--ratio`: N:n ratio(s), supports comma-separated values like `2,4,8,16`
- `--format`: Output format (table, compact, csv, json)
- `--verbose/-v`: Detailed progress output
- `--debug/-d`: Debug-level output

### Example Configurations
```bash
# Quick test
cargo bench --bench proof_system -- --system lasso --k 8 --ratio 2

# Performance comparison
cargo bench --bench proof_system -- --system logupgkr,lasso --k 8..10 --ratio 2,4,8

# Comprehensive benchmark (may take hours)
./run_benchmark_with_csv.sh
```

## Debug Logging System

The project includes a configurable logging system controlled by the `LOOKUP_LOG_LEVEL` environment variable:
- `debug` or `2`: Show all debug and info messages
- `info` or `1`: Show only info messages (default)
- `silent` or `0`: No debug output

Use `log_info!()` and `log_debug!()` macros in code for structured logging.

## Development Notes

### Git Configuration
- Main branch: `main`
- Current branch: `soundness&completeness`

### Dependencies
- Uses forked versions of `halo2` and related libraries for benchmark compatibility
- Requires `rayon` for parallel execution
- Built on `halo2_curves` with custom patches

### Testing Parallel vs Non-Parallel
Use `test_parallel_comparison.sh` to compare performance between parallel and non-parallel implementations for supported systems.

## Known Issues

### Lasso Parallel Execution Race Conditions
**Issue**: Lasso benchmarks fail with `InvalidSumcheck` errors when run in parallel with other Lasso instances.

**Symptoms**: 
- Individual Lasso tests work perfectly
- Multiple concurrent Lasso benchmarks cause sumcheck validation failures
- Error occurs at `lasso.rs:459:7` in verifier with message "Unmatched between Lasso sum_check output and query evaluation"

**Root Cause**: Race conditions in Lasso's internal state when multiple instances run concurrently. This affects:
- Memory checking prover/verifier state
- Transcript generation and challenge derivation
- Polynomial commitment operations
- Shared static state in range table creation

**Workaround**: Run Lasso benchmarks with sequential execution:
```bash
# For individual testing (clean output without verbose timing logs)
cargo bench --bench proof_system --no-default-features -- --system lasso --k 8..10

# For comprehensive benchmarking
./run_benchmark_with_csv.sh  # Now automatically handles this
```

**Status**: This is a known limitation. Lasso works correctly but cannot run multiple instances in parallel safely. The benchmark script has been updated to handle this automatically.