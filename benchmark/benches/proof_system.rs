use itertools::Itertools;
use plonkish_backend::backend::{self, baloo, cq};
use regex::Regex;
use std::collections::HashMap;
use std::{
    env::args,
    fmt::{write, Display},
    fs::{create_dir, File, OpenOptions},
    io::{stdout, Write},
    iter,
    ops::Range,
    path::Path,
    time::{Duration, Instant},
};

const OUTPUT_DIR: &str = "../target/bench";

fn main() {
    // Parse arguments
    let (systems, k_range, verbose, output_format, debug) = parse_args();
    create_output(&systems);

    // Store all benchmark results for final summary
    let mut all_results: Vec<BenchmarkResult> = Vec::new();

    // Run benchmarks for each system and k value
    k_range.for_each(|k| {
        systems.iter().for_each(|system| {
            if verbose {
                println!("→ Running {} benchmark with k = {}", system.to_string(), k);
            }

            let result = system.bench(k, verbose, debug);
            all_results.push(result);
        })
    });

    // Display results based on selected format
    match output_format {
        OutputFormat::Table => display_table_results(&all_results),
        OutputFormat::Compact => display_compact_results(&all_results),
        OutputFormat::CSV => display_csv_results(&all_results),
        OutputFormat::JSON => display_json_results(&all_results),
    }

    // Display summary statistics
    display_summary(&all_results);
}

// Struct to store benchmark results
#[derive(Debug, Clone)]
struct BenchmarkResult {
    system: System,
    k_value: usize,
    setup_time: u64,  // in milliseconds
    prove_time: u64,  // in milliseconds
    verify_time: u64, // in milliseconds
}

// Enum for different output formats
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputFormat {
    Table,   // Fancy table format
    Compact, // Compact single-line format
    CSV,     // Comma-separated values
    JSON,    // JSON format
}

fn bench_baloo(k: usize, verbose: bool, debug: bool) -> BenchmarkResult {
    // Capture and redirect detailed output if not verbose
    let timings = if !verbose && !debug {
        with_suppressed_output(|| plonkish_backend::backend::baloo::Baloo::test_baloo_by_k(k))
    } else {
        plonkish_backend::backend::baloo::Baloo::test_baloo_by_k(k)
    };

    // Write results to file
    for timing in &timings {
        writeln!(&mut System::Baloo.output(), "{}", timing).unwrap();
    }

    let all_timings = timings.join("\n");

    if debug {
        println!("\nDEBUG: Baloo raw timing output:");
        println!("{}", all_timings);
    }

    // Extract performance metrics from combined timings
    let setup_time = extract_setup_time(&all_timings, debug);
    let prove_time = extract_prove_time(&all_timings, debug);
    let verify_time = extract_verify_time(&all_timings, debug);

    if verbose || debug {
        println!(
            "Baloo times extracted - Setup: {}ms, Prove: {}ms, Verify: {}ms",
            setup_time, prove_time, verify_time
        );
    }

    // Return structured benchmark result
    BenchmarkResult {
        system: System::Baloo,
        k_value: k,
        setup_time,
        prove_time,
        verify_time,
    }
}

fn bench_CQ(k: usize, verbose: bool, debug: bool) -> BenchmarkResult {
    // Capture and redirect detailed output if not verbose
    let timings = if !verbose && !debug {
        with_suppressed_output(|| plonkish_backend::backend::cq::test_cq_by_k(k))
    } else {
        plonkish_backend::backend::cq::test_cq_by_k(k)
    };

    // Write results to file
    for timing in &timings {
        writeln!(&mut System::CQ.output(), "{}", timing).unwrap();
    }

    let all_timings = timings.join("\n");

    if debug {
        println!("\nDEBUG: CQ raw timing output:");
        println!("{}", all_timings);
    }

    // Extract performance metrics from combined timings
    let setup_time = extract_setup_time(&all_timings, debug);
    let prove_time = extract_prove_time(&all_timings, debug);
    let verify_time = extract_verify_time(&all_timings, debug);

    if verbose || debug {
        println!(
            "CQ times extracted - Setup: {}ms, Prove: {}ms, Verify: {}ms",
            setup_time, prove_time, verify_time
        );
    }

    // Return structured benchmark result
    BenchmarkResult {
        system: System::CQ,
        k_value: k,
        setup_time,
        prove_time,
        verify_time,
    }
}

// Helper function to extract timing information from output
fn extract_setup_time(timing: &str, debug: bool) -> u64 {
    if debug {
        println!("DEBUG: Searching for setup time in:\n{}", timing);
    }

    let patterns = [
        r"------------\?Setup and preprocess: (\d+)ms-----------", // 最匹配示例輸出的模式
        r"Setup and preprocess: (\d+)ms",
        r"-+\?Setup and preprocess: (\d+)ms-+",
        r"-+Setup and preprocess: (\d+)ms-+",
        r"Setup and preprocess: (\d+) ms",
        r"Setup: (\d+)ms",
    ];

    for pattern in patterns {
        let re = Regex::new(pattern).unwrap();
        for cap in re.captures_iter(timing) {
            if let Ok(value) = cap[1].parse::<u64>() {
                if debug {
                    println!("DEBUG: Found setup time {} with pattern {}", value, pattern);
                }
                return value;
            }
        }
    }

    let general_re = Regex::new(r"[Ss]etup.*?(\d+)\s*ms").unwrap();
    for cap in general_re.captures_iter(timing) {
        if let Ok(value) = cap[1].parse::<u64>() {
            if debug {
                println!("DEBUG: Found setup time {} with general pattern", value);
            }
            return value;
        }
    }

    if debug {
        println!("DEBUG: Could not find setup time");
    }

    0
}

fn extract_prove_time(timing: &str, debug: bool) -> u64 {
    if debug {
        println!("DEBUG: Searching for prove time in:\n{}", timing);
    }

    let patterns = [
        r"------------prove: (\d+)ms------------",
        r"prove: (\d+)ms",
        r"-+prove: (\d+)ms-+",
        r"prove: (\d+) ms",
        r"Prove: (\d+)ms",
    ];

    for pattern in patterns {
        let re = Regex::new(pattern).unwrap();
        for cap in re.captures_iter(timing) {
            if let Ok(value) = cap[1].parse::<u64>() {
                if debug {
                    println!("DEBUG: Found prove time {} with pattern {}", value, pattern);
                }
                return value;
            }
        }
    }

    let general_re = Regex::new(r"[Pp]rove.*?(\d+)\s*ms").unwrap();
    for cap in general_re.captures_iter(timing) {
        if let Ok(value) = cap[1].parse::<u64>() {
            if debug {
                println!("DEBUG: Found prove time {} with general pattern", value);
            }
            return value;
        }
    }

    if debug {
        println!("DEBUG: Could not find prove time");
    }

    0
}

fn extract_verify_time(timing: &str, debug: bool) -> u64 {
    if debug {
        println!("DEBUG: Searching for verify time in:\n{}", timing);
    }

    let patterns = [
        r"------------verify: (\d+)ms------------",
        r"verify: (\d+)ms",
        r"-+verify: (\d+)ms-+",
        r"verify: (\d+) ms",
        r"Verify: (\d+)ms",
    ];

    for pattern in patterns {
        let re = Regex::new(pattern).unwrap();
        for cap in re.captures_iter(timing) {
            if let Ok(value) = cap[1].parse::<u64>() {
                if debug {
                    println!(
                        "DEBUG: Found verify time {} with pattern {}",
                        value, pattern
                    );
                }
                return value;
            }
        }
    }

    let general_re = Regex::new(r"[Vv]erify.*?(\d+)\s*ms").unwrap();
    for cap in general_re.captures_iter(timing) {
        if let Ok(value) = cap[1].parse::<u64>() {
            if debug {
                println!("DEBUG: Found verify time {} with general pattern", value);
            }
            return value;
        }
    }

    if debug {
        println!("DEBUG: Could not find verify time");
    }

    0
}

// Helper function to suppress stdout output
fn with_suppressed_output<F, T>(f: F) -> T
where
    F: FnOnce() -> T,
{
    // TODO: Implement proper output suppression if needed
    // For now, just call the function directly
    f()
}

// Display results in fancy table format
fn display_table_results(results: &[BenchmarkResult]) {
    println!("\n✓ Results:");
    println!(
        "┌──────────┬─────────┬────────────────────┬──────────────┬───────────────┬────────────┐"
    );
    println!(
        "│ System   │ K-value │ Setup+Preprocess   │ Prove        │ Verify        │ Total      │"
    );
    println!(
        "├──────────┼─────────┼────────────────────┼──────────────┼───────────────┼────────────┤"
    );

    // Add data rows
    for result in results {
        let total = result.setup_time + result.prove_time + result.verify_time;
        println!(
            "│ {:<8} │ {:<7} │ {:<18} │ {:<12} │ {:<13} │ {:<10} │",
            result.system.to_string(),
            result.k_value,
            format!("{}ms", result.setup_time),
            format!("{}ms", result.prove_time),
            format!("{}ms", result.verify_time),
            format!("{}ms", total)
        );
    }

    println!(
        "└──────────┴─────────┴────────────────────┴──────────────┴───────────────┴────────────┘"
    );
}

// Display results in compact single-line format
fn display_compact_results(results: &[BenchmarkResult]) {
    println!("\n✓ Results:");
    println!(
        "{:10} {:7} {:12} {:10} {:10} {:10}",
        "System", "K", "Setup (ms)", "Prove (ms)", "Verify (ms)", "Total (ms)"
    );
    println!("{}", "-".repeat(60));

    for result in results {
        let total = result.setup_time + result.prove_time + result.verify_time;
        println!(
            "{:10} {:7} {:12} {:10} {:10} {:10}",
            result.system.to_string(),
            result.k_value,
            result.setup_time,
            result.prove_time,
            result.verify_time,
            total
        );
    }
}

// Display results in CSV format
fn display_csv_results(results: &[BenchmarkResult]) {
    println!("System,K,SetupTime,ProveTime,VerifyTime,TotalTime");
    for result in results {
        let total = result.setup_time + result.prove_time + result.verify_time;
        println!(
            "{},{},{},{},{},{}",
            result.system,
            result.k_value,
            result.setup_time,
            result.prove_time,
            result.verify_time,
            total
        );
    }
}

// Display results in JSON format
fn display_json_results(results: &[BenchmarkResult]) {
    println!("[");
    for (i, result) in results.iter().enumerate() {
        let total = result.setup_time + result.prove_time + result.verify_time;
        println!("  {{");
        println!("    \"system\": \"{}\",", result.system);
        println!("    \"k_value\": {},", result.k_value);
        println!("    \"setup_time\": {},", result.setup_time);
        println!("    \"prove_time\": {},", result.prove_time);
        println!("    \"verify_time\": {},", result.verify_time);
        println!("    \"total_time\": {}", total);
        if i < results.len() - 1 {
            println!("  }},");
        } else {
            println!("  }}");
        }
    }
    println!("]");
}

// Display summary statistics
fn display_summary(results: &[BenchmarkResult]) {
    // Create system-specific summaries
    let mut system_summaries: HashMap<System, Vec<&BenchmarkResult>> = HashMap::new();
    for result in results {
        system_summaries
            .entry(result.system)
            .or_default()
            .push(result);
    }

    println!("\n📊 Benchmark Summary:");
    println!("┌──────────┬─────────────┬─────────────┬────────────┬─────────────┐");
    println!("│ System   │ K Range     │ Avg Setup   │ Avg Prove  │ Avg Verify  │");
    println!("├──────────┼─────────────┼─────────────┼────────────┼─────────────┤");

    for (system, sys_results) in system_summaries {
        // Calculate averages
        let avg_setup =
            sys_results.iter().map(|r| r.setup_time).sum::<u64>() / sys_results.len() as u64;
        let avg_prove =
            sys_results.iter().map(|r| r.prove_time).sum::<u64>() / sys_results.len() as u64;
        let avg_verify =
            sys_results.iter().map(|r| r.verify_time).sum::<u64>() / sys_results.len() as u64;

        // Get k range
        let min_k = sys_results.iter().map(|r| r.k_value).min().unwrap_or(0);
        let max_k = sys_results.iter().map(|r| r.k_value).max().unwrap_or(0);
        let k_range = if min_k == max_k {
            format!("{}", min_k)
        } else {
            format!("{}..{}", min_k, max_k)
        };

        println!(
            "│ {:<8} │ {:<11} │ {:<11} │ {:<10} │ {:<11} │",
            system.to_string(),
            k_range,
            format!("{}ms", avg_setup),
            format!("{}ms", avg_prove),
            format!("{}ms", avg_verify)
        );
    }

    println!("└──────────┴─────────────┴─────────────┴────────────┴─────────────┘");
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum System {
    CQ,
    Baloo,
}

impl System {
    fn all() -> Vec<System> {
        vec![System::CQ, System::Baloo]
    }

    fn output_path(&self) -> String {
        format!("{OUTPUT_DIR}/{self}")
    }

    fn output(&self) -> File {
        OpenOptions::new()
            .append(true)
            .open(self.output_path())
            .unwrap()
    }

    // Benchmark a system with k value
    fn bench(&self, k: usize, verbose: bool, debug: bool) -> BenchmarkResult {
        match self {
            System::Baloo => bench_baloo(k, verbose, debug),
            System::CQ => bench_CQ(k, verbose, debug),
        }
    }
}

impl Display for System {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            System::Baloo => write!(f, "Baloo"),
            System::CQ => write!(f, "CQ"),
        }
    }
}

fn parse_args() -> (Vec<System>, Range<usize>, bool, OutputFormat, bool) {
    let (systems, k_range, verbose, output_format, debug) =
        args().chain(Some("".to_string())).tuple_windows().fold(
            (Vec::new(), 20..26, false, OutputFormat::Table, false),
            |(mut systems, mut k_range, mut verbose, mut output_format, mut debug),
             (key, value)| {
                match key.as_str() {
                    "--system" => match value.as_str() {
                        "all" => systems = System::all(),
                        "cq" => systems.push(System::CQ),
                        "CQ" => systems.push(System::CQ),
                        "Baloo" => systems.push(System::Baloo),
                        "baloo" => systems.push(System::Baloo),
                        _ => panic!("system should be one of {{all, cq, baloo}}"),
                    },

                    "--k" => {
                        if let Some((start, end)) = value.split_once("..") {
                            k_range = start.parse().expect("k range start to be usize")
                                ..end.parse().expect("k range end to be usize");
                        } else {
                            k_range.start = value.parse().expect("k to be usize");
                            k_range.end = k_range.start + 1;
                        }
                    }
                    "--verbose" | "-v" => {
                        verbose = true;
                    }
                    "--debug" | "-d" => {
                        debug = true;
                    }
                    "--format" | "-f" => match value.to_lowercase().as_str() {
                        "table" => output_format = OutputFormat::Table,
                        "compact" => output_format = OutputFormat::Compact,
                        "csv" => output_format = OutputFormat::CSV,
                        "json" => output_format = OutputFormat::JSON,
                        _ => panic!("format should be one of {{table, compact, csv, json}}"),
                    },
                    _ => {}
                }
                (systems, k_range, verbose, output_format, debug)
            },
        );

    let mut systems = systems.into_iter().sorted().dedup().collect_vec();
    if systems.is_empty() {
        systems = System::all();
    };

    (systems, k_range, verbose, output_format, debug)
}

fn create_output(systems: &[System]) {
    if !Path::new(OUTPUT_DIR).exists() {
        create_dir(OUTPUT_DIR).unwrap();
    }
    for system in systems {
        File::create(system.output_path()).unwrap();
    }
}
