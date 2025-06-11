# Lookup Argument Benchmarks

Benchmark different lookup argument implementations with configurable debug logging and **parallel execution**.

## Table of Contents
- [Benchmarking](#benchmarking)
- [Parallel Execution](#parallel-execution)
- [Debug Logging System](#debug-logging-system)
- [Basic Usage Examples](#basic-usage-examples)
- [Advanced Usage Examples](#advanced-usage-examples)
- [New Features](#new-features)
- [Parameter Reference](#parameter-reference)

## Benchmarking

```sh
Usage: cargo bench --bench proof_system -- [OPTIONS]

Options:
  --system <s>    Lookup system(s) to run. [possible values: all, cq, baloo, logupgkr, plookup, caulk, lasso]
  --k <K>              (Range of) log size of the lookup table.
  --ratio, --n-to-n-ratio <RATIO>  Ratio between table size (N) and lookup size (n). Default is 2. Supports comma-separated multiple values.
  --verbose, -v        Enable verbose output.
  --debug, -d          Enable debug mode with detailed output.
  --format, -f <FORMAT> Output format. [possible values: table, compact, csv, json]
```

## Parallel Execution

**🚀 NEW: Automatic Parallel Benchmark Execution + Individual Algorithm Parallelization**

The benchmark system now supports two levels of parallelization:

### 1. Benchmark-Level Parallelization
Multiple benchmark tasks run in parallel using Rust's `rayon` crate, providing significant performance improvements:

#### Key Features:
- **Automatic parallelization**: All benchmark combinations (systems × k-values × ratios) run in parallel
- **Thread-safe execution**: Safe concurrent access to shared resources and file outputs
- **Progress tracking**: Real-time progress updates with thread IDs for debugging
- **Consistent results**: Results are automatically sorted for consistent display order
- **Load balancing**: Work is automatically distributed across available CPU cores

#### Performance Benefits:
- **Faster execution**: Benchmark suites complete significantly faster than serial execution
- **CPU utilization**: Better utilization of multi-core systems
- **Time savings**: Particularly beneficial for comprehensive benchmarks with multiple systems and parameter combinations

### 2. Algorithm-Level Parallelization
Individual proof systems can now run in parallel or non-parallel mode based on the `parallel` feature flag:

#### Supported Systems:
- **LogupGKR**: ✅ Full parallel/non-parallel support
- **Lasso**: ✅ Full parallel/non-parallel support (NEW!)
- **HyperPlonk**: ✅ Always parallel (uses `util::parallel` module)
- **Plookup**: ❌ Always serial
- **CQ**: ❌ Always serial  
- **Baloo**: ❌ Always serial
- **Caulk**: ❌ Always serial

#### Performance Comparison:
You can now compare the performance difference between parallel and non-parallel implementations:

```bash
# Run with parallel algorithms enabled (default)
cargo bench --bench proof_system -- --system lasso --k 8..12 --format table

# Run with parallel algorithms disabled
cargo bench --bench proof_system --no-default-features --features timer -- --system lasso --k 8..12 --format table
```

### Example Parallel Output:
```bash
→ Running 24 benchmark tasks in parallel...
→ Systems: ["CQ", "Baloo", "LogupGKR", "Plookup", "Caulk", "Lasso"]
→ K values: [8, 9]
→ N:n ratios: [2, 4]
→ [Thread ThreadId(6)] Running CQ benchmark with k = 8, N:n ratio = 2
→ [Thread ThreadId(5)] Running Plookup benchmark with k = 8, N:n ratio = 2
→ [Thread ThreadId(9)] Running Baloo benchmark with k = 8, N:n ratio = 2
...
✓ [Thread ThreadId(7)] Completed LogupGKR benchmark with k = 8, N:n ratio = 2
✓ [Thread ThreadId(2)] Completed Lasso benchmark with k = 8, N:n ratio = 2
→ All 24 benchmark tasks completed!
```

### Technical Implementation:
- Uses `rayon::par_iter()` for parallel task execution
- Thread-safe result collection with `Arc<Mutex<Vec<BenchmarkResult>>>`
- Maintains file output integrity across concurrent writes
- Preserves all original functionality while adding parallelization

## Debug Logging System

This project now supports controllable debug message output. You can control the displayed message level through environment variables.

### Log Levels

The system supports three log levels:

1. **Silent** (`silent` or `0`) - No debug messages
2. **Info** (`info` or `1`) - Only important information (default level)
3. **Debug** (`debug` or `2`) - Detailed debug information

### Setting Environment Variables

```bash
# Set to debug mode (show all messages)
export LOOKUP_LOG_LEVEL=debug

# Set to info mode (only important messages, default)
export LOOKUP_LOG_LEVEL=info

# Set to silent mode (no debug messages)
export LOOKUP_LOG_LEVEL=silent
```

### Logging Usage Examples

```bash
# Run lasso test in debug mode
LOOKUP_LOG_LEVEL=debug cargo test test_lasso

# Run tests in silent mode (no debug output)
LOOKUP_LOG_LEVEL=silent cargo test

# Run with default info mode
cargo test

# Run benchmarks without debug messages
LOOKUP_LOG_LEVEL=silent cargo bench

# Run benchmarks with detailed debug information
LOOKUP_LOG_LEVEL=debug cargo bench
```

### Code Usage

If you're extending the code, you can use the following macros:

```rust
use plonkish_backend::{log_info, log_debug};

// Initialize logging (usually called once at program start)
plonkish_backend::logging::init_logging();

// Info level messages (important information)
log_info!("Processing {} elements", count);

// Debug level messages (detailed debugging information)
log_debug!("Internal state: {:?}", state);

// Conditional debug println (compatible with existing code)
debug_println!("DEBUG: Value = {}", value);
```

### Logging Output Examples

**Before (all debug messages shown)**:
```
DEBUG Prover memory index 0: id_value = 0x0000...
DEBUG Prover memory index 0: init = 0x22fc...
DEBUG Prover memory index 1: id_value = 0x0001...
DEBUG MemoryCheckingProver::prove: num_memories = 2
...massive debug output...
```

**After (controllable)**:

**Info mode (default)**:
```
[INFO] Processing Lasso proof with k=8
[INFO] Verification completed successfully
```

**Debug mode**:
```
[INFO] Processing Lasso proof with k=8
[DEBUG] Prover memory index 0: id_value = 0x0000...
[DEBUG] Prover memory index 0: init = 0x22fc...
[DEBUG] MemoryCheckingProver::prove: num_memories = 2
...detailed debug output...
[INFO] Verification completed successfully
```

**Silent mode**:
```
(only test results, no debug messages)
```

### Logging Benefits

1. **Clean default output**: Only shows important messages by default, not overwhelmed by debug messages
2. **Flexible debugging**: Easy to enable detailed output when debugging is needed
3. **Performance friendly**: Debug formatting and output are skipped in silent/info modes
4. **Backward compatible**: Existing code continues to work normally
5. **Environment variable control**: Change log level without recompilation

## Basic Usage Examples

### Single System with Basic Parameters
```sh
cargo bench --bench proof_system -- \
    --system all \
    --k 8 \
    --format table
```

### Range of Table Sizes
```sh
cargo bench --bench proof_system -- \
    --system logupgkr \
    --k 8..12 \
    --format compact
```

### Custom N:n Ratio
```sh
cargo bench --bench proof_system -- \
    --system all \
    --k 10 \
    --ratio 8 \
    --format table
```

## Advanced Usage Examples

### Multiple N:n Ratios Testing
```sh
# Test single ratio
cargo bench --bench proof_system -- \
    --system all \
    --k 8..12 \
    --ratio 4 \
    --format compact

# Test multiple ratios (comma-separated)
cargo bench --bench proof_system -- \
    --system all \
    --k 8..12 \
    --ratio 2,4,8,16 \
    --format compact
```

### Proof Size Analysis
```sh
# Table format with detailed information including proof size
cargo bench --bench proof_system -- \
    --system all \
    --k 8..12 \
    --ratio 2,4,8 \
    --format table

# Compact format also displays proof size
cargo bench --bench proof_system -- \
    --system all \
    --k 8..12 \
    --ratio 2,4,8 \
    --format compact
```

### Comprehensive Benchmarking
```sh
# Full benchmark test (runs in parallel automatically)
cargo bench --bench proof_system -- \
    --system logupgkr,plookup,lasso \
    --k 8..15 \
    --ratio 1,2,4,8,16,32 \
    --format table \
    --verbose

# Specific system deep analysis (parallel execution)
cargo bench --bench proof_system -- \
    --system lasso \
    --k 10..14 \
    --ratio 2,4,8,16 \
    --format csv \
    --debug
```

### Parallel Performance Examples
```sh
# Small parallel test (6 tasks: 6 systems × 1 k-value × 1 ratio)
cargo bench --bench proof_system -- \
    --system all \
    --k 8 \
    --ratio 2 \
    --verbose

# Medium parallel test (24 tasks: 6 systems × 2 k-values × 2 ratios)
cargo bench --bench proof_system -- \
    --system all \
    --k 8..10 \
    --ratio 2,4 \
    --verbose

# Large parallel test (120 tasks: 5 systems × 4 k-values × 6 ratios)
cargo bench --bench proof_system -- \
    --system cq,baloo,logupgkr,plookup,lasso \
    --k 8..12 \
    --ratio 1,2,4,8,16,32 \
    --format compact
```

## New Features

### 1. **Parallel Benchmark Execution** 🆕
- **Automatic parallelization**: All benchmark tasks run concurrently using `rayon`
- **Significant speedup**: Better utilization of multi-core systems
- **Thread-safe operation**: Safe concurrent file writing and result collection
- **Real-time progress**: Thread-level progress tracking with IDs
- **Consistent output**: Results automatically sorted for predictable display order

### 2. **Debug Logging System**
- **Three log levels**: Silent, Info (default), Debug
- **Environment variable control**: `LOOKUP_LOG_LEVEL=debug/info/silent`
- **Performance optimized**: Debug formatting skipped when disabled
- **Thread-safe**: Global configuration with safe concurrent access

### 3. **Multiple N:n Ratio Testing**
- **Before**: Only single N:n ratio testing
- **Now**: Support comma-separated multiple ratios like `--ratio 2,4,8`
- **Purpose**: Comprehensive analysis of how different table-to-lookup size ratios affect performance across systems

### 4. **Proof Size Measurement**
- **New Field**: `proof_size` (in bytes)
- **Display Format**: 
  - Less than 1KB: "XXXnB" 
  - Greater than 1KB: "X.XKB"
  - No data available: "N/A"
- **Purpose**: Compare proof size efficiency across different proving systems

### 5. **Enhanced Output Formats**
- **Table Format**: Added "Proof Size" column
- **Compact Format**: Includes proof size information
- **Compatibility**: All existing formats (CSV, JSON) support the new field

## Table Size (N) vs Lookup Size (n)

By default, the benchmark uses K to define table_size = 2^K (N) and lookup_size = table_size / ratio (n). The default ratio is 2, giving N:n = 2:1.

This default 2:1 ratio might not fully represent scenarios where sublinear-N protocols like Baloo and CQ are expected to perform best. These protocols shine most when N is much larger than n (e.g., N = 8n or more).

Using the `--ratio` parameter, you can:
- **Single ratio**: `--ratio 8` for N = 8n
- **Multiple ratios**: `--ratio 2,4,8,16` to test various scenarios simultaneously

This allows better showcasing of performance characteristics under different workload ratios.

## Parameter Reference

### Benchmark Parameters
- `--system`: Choose systems to benchmark (all, cq, baloo, logupgkr, plookup, caulk, lasso)
- `--k`: Set k value range, e.g., `8..12` or single value `10`
- `--ratio` or `--n-to-n-ratio`: Set N:n ratio(s), supports single value or comma-separated multiple values
- `--format`: Output format (table, compact, csv, json)
- `--verbose` or `-v`: Display detailed output
- `--debug` or `-d`: Display debug information

### Logging Environment Variables
- `LOOKUP_LOG_LEVEL=debug` or `LOOKUP_LOG_LEVEL=2`: Show all debug and info messages
- `LOOKUP_LOG_LEVEL=info` or `LOOKUP_LOG_LEVEL=1`: Show only info messages (default)
- `LOOKUP_LOG_LEVEL=silent` or `LOOKUP_LOG_LEVEL=0`: No debug output

Results will be formatted according to the specified output format and displayed in the terminal, now including timing information and proof sizes for comprehensive analysis.

## Acknowledgements

This work is forked from Plonkish by Han
