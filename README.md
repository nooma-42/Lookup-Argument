# Lookup Argument Benchmarks

Benchmark different lookup argument implementations.

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

### 1. **Multiple N:n Ratio Testing**
- **Before**: Only single N:n ratio testing
- **Now**: Support comma-separated multiple ratios like `--ratio 2,4,8`
- **Purpose**: Comprehensive analysis of how different table-to-lookup size ratios affect performance across systems

### 2. **Proof Size Measurement**
- **New Field**: `proof_size` (in bytes)
- **Display Format**: 
  - Less than 1KB: "XXXnB" 
  - Greater than 1KB: "X.XKB"
  - No data available: "N/A"
- **Purpose**: Compare proof size efficiency across different proving systems

### 3. **Enhanced Output Formats**
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

- `--system`: Choose systems to benchmark (all, cq, baloo, logupgkr, plookup, caulk, lasso)
- `--k`: Set k value range, e.g., `8..12` or single value `10`
- `--ratio` or `--n-to-n-ratio`: Set N:n ratio(s), supports single value or comma-separated multiple values
- `--format`: Output format (table, compact, csv, json)
- `--verbose` or `-v`: Display detailed output
- `--debug` or `-d`: Display debug information

Results will be formatted according to the specified output format and displayed in the terminal, now including timing information and proof sizes for comprehensive analysis.

## Acknowledgements

This work is forked from Plonkish by Han
