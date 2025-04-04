# Lookup Argument Benchmarks

Benchmark different lookup argument implementations.

```sh
Usage: cargo bench --bench proof_system -- [OPTIONS]

Options:
  --system <s>    Lookup system(s) to run. [possible values: all, cq, baloo, logupgkr]
  --k <K>              (Range of) log size of the lookup table.
  --ratio, --n-to-n-ratio <RATIO>  Ratio between table size (N) and lookup size (n). Default is 2.
  --verbose, -v        Enable verbose output.
  --debug, -d          Enable debug mode with detailed output.
  --format, -f <FORMAT> Output format. [possible values: table, compact, csv, json]
```

For example, to compare all lookup argument systems with a table size of 2^8:

```sh
cargo bench --bench proof_system -- \
    --system all \
    --k 8 \
    --format table
```

To run benchmarks across a range of table sizes for LogupGKR only:

```sh
cargo bench --bench proof_system -- \
    --system logupgkr \
    --k 8..12 \
    --format compact
```

To compare systems with a larger table-to-lookup ratio (e.g., N = 8n, better showcasing sublinear-N protocols):

```sh
cargo bench --bench proof_system -- \
    --system all \
    --k 10 \
    --ratio 8 \
    --format table
```

## Table Size (N) vs Lookup Size (n)

By default, the benchmark uses K to define table_size = 2^K (N) and lookup_size = table_size / 2 (n), giving a ratio of N:n = 2:1. 

This default 2:1 ratio might not fully represent scenarios where sublinear-N protocols like Baloo and CQ are expected to perform best. These protocols shine most when N is much larger than n (e.g., N = 8n or more).

Using the `--ratio` parameter, you can adjust this relationship to better showcase the performance characteristics of different protocols under various workload ratios.

Results will be formatted according to the specified output format and displayed in the terminal.

## Acknowledgements

This work is forked from Plonkish by Han
