# Lookup Argument Benchmarks

Benchmark different lookup argument implementations.

```sh
Usage: cargo bench --bench proof_system -- [OPTIONS]

Options:
  --system <SYSTEM>    Lookup system(s) to run. [possible values: all, cq, baloo, logupgkr]
  --k <K>              (Range of) log size of the lookup table.
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

Results will be formatted according to the specified output format and displayed in the terminal.

## Acknowledgements

This work is forked from Plonkish by Han
