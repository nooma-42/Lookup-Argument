# Lookup Argument Benchmarks

Benchmark different lookup argument implementations with configurable debug logging.

## Table of Contents
- [Benchmarking](#benchmarking)
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
# Full benchmark test
cargo bench --bench proof_system -- \
    --system logupgkr,plookup,lasso \
    --k 8..15 \
    --ratio 1,2,4,8,16,32 \
    --format table \
    --verbose

# Specific system deep analysis
cargo bench --bench proof_system -- \
    --system lasso \
    --k 10..14 \
    --ratio 2,4,8,16 \
    --format csv \
    --debug
```

## New Features

### 1. **Debug Logging System**
- **Three log levels**: Silent, Info (default), Debug
- **Environment variable control**: `LOOKUP_LOG_LEVEL=debug/info/silent`
- **Performance optimized**: Debug formatting skipped when disabled
- **Thread-safe**: Global configuration with safe concurrent access

### 2. **Multiple N:n Ratio Testing**
- **Before**: Only single N:n ratio testing
- **Now**: Support comma-separated multiple ratios like `--ratio 2,4,8`
- **Purpose**: Comprehensive analysis of how different table-to-lookup size ratios affect performance across systems

### 3. **Proof Size Measurement**
- **New Field**: `proof_size` (in bytes)
- **Display Format**: 
  - Less than 1KB: "XXXnB" 
  - Greater than 1KB: "X.XKB"
  - No data available: "N/A"
- **Purpose**: Compare proof size efficiency across different proving systems

### 4. **Enhanced Output Formats**
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
