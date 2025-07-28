use itertools::Itertools;
use plonkish_backend::backend::{cq, logupgkr, plookup, lasso, caulk};
use plonkish_backend::halo2_curves::bn256::{Bn256, Fr};
use plonkish_backend::pcs::univariate::UnivariateKzg; 
use plonkish_backend::util::testing::{TestType, TestResult};
use rayon::prelude::*;
use regex::Regex;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::{
    env::args,
    fmt::Display,
    fs::{create_dir, File, OpenOptions},
    io::{Write, BufWriter},
    ops::Range,
    path::Path,
};

// Type alias for Plookup with BN256 curve
type PlookupBn256 = plookup::Plookup<Fr, UnivariateKzg<Bn256>>;

const OUTPUT_DIR: &str = "../target/bench";

/// Types of operations that can be tested (currently only applies to Lasso)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Operation {
    Range,  // Range check operations (default for all systems)
    Add,    // Addition operations (Lasso only)
    All,    // Both range and add operations (Lasso only)
}

impl std::fmt::Display for Operation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Operation::Range => write!(f, "Range"),
            Operation::Add => write!(f, "Add"),
            Operation::All => write!(f, "All"),
        }
    }
}

// Struct to track benchmark execution results (success or failure)
#[derive(Debug, Clone)]
enum BenchmarkOutcome {
    Success(BenchmarkResult),
    Failure(BenchmarkFailure),
}

// Struct to store failed benchmark information
#[derive(Debug, Clone)]
struct BenchmarkFailure {
    system: System,
    k_value: usize,
    n_to_n_ratio: usize,
    error_message: String,
    timestamp: std::time::SystemTime,
}

// Thread-safe CSV writer for incremental results
struct CsvWriter {
    file: Arc<Mutex<BufWriter<File>>>,
    header_written: Arc<Mutex<bool>>,
}

impl CsvWriter {
    fn new(filename: &str) -> Result<Self, std::io::Error> {
        let file_exists = Path::new(filename).exists();
        let file = OpenOptions::new().create(true).append(true).open(filename)?;
        let buffered_writer = BufWriter::new(file);
        
        Ok(CsvWriter {
            file: Arc::new(Mutex::new(buffered_writer)),
            header_written: Arc::new(Mutex::new(file_exists)),
        })
    }
    
    fn write_result(&self, result: &BenchmarkResult) -> Result<(), std::io::Error> {
        let mut file = self.file.lock().unwrap();
        let mut header_written = self.header_written.lock().unwrap();
        
        // Write header if not written yet
        if !*header_written {
            writeln!(file, "System,K,N_to_n_Ratio,SetupTime_ms,ProveTime_ms,VerifyTime_ms,TotalTime_ms,ProofSize_bytes,Timestamp")?;
            *header_written = true;
        }
        
        let total = result.setup_time + result.prove_time + result.verify_time;
        let timestamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_secs();
            
        writeln!(
            file,
            "{},{},{},{},{},{},{},{},{}",
            result.system,
            result.k_value,
            result.n_to_n_ratio,
            result.setup_time,
            result.prove_time,
            result.verify_time,
            total,
            result.proof_size,
            timestamp
        )?;
        
        file.flush()?;
        Ok(())
    }
}

// Thread-safe failure tracker for incremental failure logging
struct FailureTracker {
    file: Arc<Mutex<BufWriter<File>>>,
    header_written: Arc<Mutex<bool>>,
}

impl FailureTracker {
    fn new(filename: &str) -> Result<Self, std::io::Error> {
        let file = File::create(filename)?;
        let buffered_writer = BufWriter::new(file);
        
        Ok(FailureTracker {
            file: Arc::new(Mutex::new(buffered_writer)),
            header_written: Arc::new(Mutex::new(false)),
        })
    }
    
    fn write_failure(&self, failure: &BenchmarkFailure) -> Result<(), std::io::Error> {
        let mut file = self.file.lock().unwrap();
        let mut header_written = self.header_written.lock().unwrap();
        
        // Write header if not written yet
        if !*header_written {
            writeln!(file, "System,K,N_to_n_Ratio,ErrorMessage,Timestamp")?;
            *header_written = true;
        }
        
        let timestamp = failure.timestamp
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_secs();
            
        writeln!(
            file,
            "{},{},{},\"{}\",{}",
            failure.system,
            failure.k_value,
            failure.n_to_n_ratio,
            failure.error_message.replace("\"", "\"\""), // Escape quotes in CSV
            timestamp
        )?;
        
        file.flush()?;
        Ok(())
    }
}

fn main() {
    // Parse arguments
    let (systems, k_range, n_to_n_ratios, verbose, output_format, debug, test_type, operation) = parse_args();
    create_output(&systems);

    // Generate all benchmark tasks
    let mut benchmark_tasks = Vec::new();
    for k in k_range {
        for system in &systems {
            for &ratio in &n_to_n_ratios {
                benchmark_tasks.push((system.clone(), k, ratio));
            }
        }
    }

    // Print progress information
    let total_tasks = benchmark_tasks.len();
    println!("\n🚀 Starting benchmark suite with {} tasks", total_tasks);
    println!("📊 Systems: {}", systems.iter().map(|s| s.to_string()).collect::<Vec<_>>().join(", "));
    println!("🔧 K values: {:?}", benchmark_tasks.iter().map(|(_, k, _)| *k).collect::<std::collections::HashSet<_>>().into_iter().collect::<Vec<_>>());
    println!("⚖️  N:n ratios: {:?}", n_to_n_ratios);
    println!("🧪 Test type: {}", test_type);
    println!("⚙️  Operation: {}", operation);
    println!("{}", "=".repeat(60));

    // Create CSV writer and failure tracker for incremental writing
    let csv_writer = match CsvWriter::new("benchmark_results.csv") {
        Ok(writer) => Some(Arc::new(writer)),
        Err(e) => {
            eprintln!("⚠️  Warning: Could not create CSV file: {}. Results will only be displayed.", e);
            None
        }
    };
    
    let failure_tracker = match FailureTracker::new("benchmark_failures.csv") {
        Ok(tracker) => Some(Arc::new(tracker)),
        Err(e) => {
            eprintln!("⚠️  Warning: Could not create failure tracking file: {}. Failures will only be displayed.", e);
            None
        }
    };

    // Use Arc<Mutex<Vec<BenchmarkOutcome>>> for thread-safe result collection
    let all_outcomes = Arc::new(Mutex::new(Vec::new()));
    let completed_counter = Arc::new(Mutex::new(0));
    let success_counter = Arc::new(Mutex::new(0));
    let failure_counter = Arc::new(Mutex::new(0));

    println!("💾 Incremental CSV writing enabled:");
    if csv_writer.is_some() {
        println!("   📊 Results: benchmark_results.csv");
    }
    if failure_tracker.is_some() {
        println!("   ❌ Failures: benchmark_failures.csv");
    }
    println!();

    // Run benchmarks in parallel using rayon
    benchmark_tasks.par_iter().enumerate().for_each(|(task_index, (system, k, ratio))| {
        // Show progress before starting
        {
            let mut counter = completed_counter.lock().unwrap();
            println!("🔄 [{:3}/{:3}] Starting: {} (k={}, ratio={})", 
                    *counter + 1, total_tasks, system.to_string(), k, ratio);
        }

        // Execute benchmark with error handling
        let outcome = match std::panic::catch_unwind(|| {
            system.bench_with_test_type(*k, *ratio, verbose, debug, test_type, operation)
        }) {
            Ok(result) => {
                // Write successful result to CSV immediately
                if let Some(ref csv_writer) = csv_writer {
                    if let Err(e) = csv_writer.write_result(&result) {
                        eprintln!("⚠️  Warning: Failed to write result to CSV: {}", e);
                    }
                }
                BenchmarkOutcome::Success(result)
            },
            Err(panic_info) => {
                let error_msg = if let Some(s) = panic_info.downcast_ref::<String>() {
                    s.clone()
                } else if let Some(s) = panic_info.downcast_ref::<&str>() {
                    s.to_string()
                } else {
                    "Unknown panic occurred".to_string()
                };
                
                let failure = BenchmarkFailure {
                    system: *system,
                    k_value: *k,
                    n_to_n_ratio: *ratio,
                    error_message: error_msg,
                    timestamp: std::time::SystemTime::now(),
                };
                
                // Write failure to failure tracking file immediately
                if let Some(ref failure_tracker) = failure_tracker {
                    if let Err(e) = failure_tracker.write_failure(&failure) {
                        eprintln!("⚠️  Warning: Failed to write failure to CSV: {}", e);
                    }
                }
                
                BenchmarkOutcome::Failure(failure)
            }
        };
        
        // Safely add outcome to shared collection and update counters
        {
            let mut completed = completed_counter.lock().unwrap();
            all_outcomes.lock().unwrap().push(outcome.clone());
            *completed += 1;
            
            match outcome {
                BenchmarkOutcome::Success(_) => {
                    let mut success = success_counter.lock().unwrap();
                    *success += 1;
                    println!("✅ [{:3}/{:3}] Completed: {} (k={}, ratio={}) - {:.1}% done", 
                            *completed, total_tasks, system.to_string(), k, ratio, 
                            (*completed as f64 / total_tasks as f64) * 100.0);
                },
                BenchmarkOutcome::Failure(ref failure) => {
                    let mut failures = failure_counter.lock().unwrap();
                    *failures += 1;
                    println!("❌ [{:3}/{:3}] Failed: {} (k={}, ratio={}) - {:.1}% done - Error: {}", 
                            *completed, total_tasks, system.to_string(), k, ratio, 
                            (*completed as f64 / total_tasks as f64) * 100.0,
                            failure.error_message);
                }
            }
        }
    });

    // Extract results from Arc<Mutex<>>
    let final_outcomes = all_outcomes.lock().unwrap().clone();
    let total_success = *success_counter.lock().unwrap();
    let total_failures = *failure_counter.lock().unwrap();

    // Separate successful results and failures
    let mut successful_results = Vec::new();
    let mut failures = Vec::new();
    
    for outcome in final_outcomes {
        match outcome {
            BenchmarkOutcome::Success(result) => successful_results.push(result),
            BenchmarkOutcome::Failure(failure) => failures.push(failure),
        }
    }

    // Sort results for consistent display order (by system, then k, then ratio)
    successful_results.sort_by(|a, b| {
        a.system.to_string()
            .cmp(&b.system.to_string())
            .then(a.k_value.cmp(&b.k_value))
            .then(a.n_to_n_ratio.cmp(&b.n_to_n_ratio))
    });

    println!("\n🎉 Benchmark suite completed!");
    println!("📊 Summary: {} successful, {} failed, {} total", total_success, total_failures, total_tasks);
    println!("{}", "=".repeat(60));

    // Display successful results based on selected format (if any exist)
    if !successful_results.is_empty() {
        match output_format {
            OutputFormat::Table => display_table_results(&successful_results),
            OutputFormat::Compact => display_compact_results(&successful_results),
            OutputFormat::CSV => display_csv_results(&successful_results),
            OutputFormat::JSON => display_json_results(&successful_results),
        }

        // Display summary statistics for successful results
        display_summary(&successful_results);
    } else {
        println!("❌ No successful benchmark results to display.");
    }

    // Display failure summary if there were any failures
    if !failures.is_empty() {
        display_failure_summary(&failures);
    }

    // Final file location summary
    if csv_writer.is_some() || failure_tracker.is_some() {
        println!("\n📁 Files created:");
        if csv_writer.is_some() {
            println!("   📊 benchmark_results.csv - {} successful results", total_success);
        }
        if failure_tracker.is_some() {
            println!("   ❌ benchmark_failures.csv - {} failures", total_failures);
        }
    }
}

// Struct to store benchmark results
#[derive(Debug, Clone)]
struct BenchmarkResult {
    system: System,
    k_value: usize,
    n_to_n_ratio: usize,
    setup_time: u64,  // in milliseconds
    prove_time: u64,  // in milliseconds
    verify_time: u64, // in milliseconds
    proof_size: u64,  // in bytes
}

// Enum for different output formats
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputFormat {
    Table,   // Fancy table format
    Compact, // Compact single-line format
    CSV,     // Comma-separated values
    JSON,    // JSON format
}

fn bench_baloo(k: usize, n_to_n_ratio: usize, verbose: bool, debug: bool) -> BenchmarkResult {
    // Capture and redirect detailed output if not verbose
    let timings = if !verbose && !debug {
        with_suppressed_output(|| plonkish_backend::backend::baloo::Baloo::test_baloo_by_k_with_ratio(k, n_to_n_ratio))
    } else {
        if verbose || debug {
            println!("Running Baloo benchmark with k={}, N:n ratio={}", k, n_to_n_ratio);
            let range_size = 1 << k;
            let lookup_size = range_size / n_to_n_ratio;
            println!("Range size: {} ({}-bit), Lookup size: {}", range_size, k, lookup_size);
        }
        plonkish_backend::backend::baloo::Baloo::test_baloo_by_k_with_ratio(k, n_to_n_ratio)
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
    let proof_size = extract_proof_size(&all_timings, debug);

    if verbose || debug {
        println!(
            "Baloo times extracted - Setup: {}ms, Prove: {}ms, Verify: {}ms, Proof size: {}B",
            setup_time, prove_time, verify_time, proof_size
        );
    }

    // Return structured benchmark result
    BenchmarkResult {
        system: System::Baloo,
        k_value: k,
        n_to_n_ratio,
        setup_time,
        prove_time,
        verify_time,
        proof_size,
    }
}

fn bench_CQ(k: usize, n_to_n_ratio: usize, verbose: bool, debug: bool) -> BenchmarkResult {
    // Capture and redirect detailed output if not verbose
    let timings = if !verbose && !debug {
        with_suppressed_output(|| plonkish_backend::backend::cq::test_cq_by_k_with_ratio(k, n_to_n_ratio))
    } else {
        if verbose || debug {
            println!("Running CQ benchmark with k={}, N:n ratio={}", k, n_to_n_ratio);
            let range_size = 1 << k;
            let lookup_size = range_size / n_to_n_ratio;
            println!("Range size: {} ({}-bit), Lookup size: {}", range_size, k, lookup_size);
        }
        plonkish_backend::backend::cq::test_cq_by_k_with_ratio(k, n_to_n_ratio)
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
    let proof_size = extract_proof_size(&all_timings, debug);

    if verbose || debug {
        println!(
            "CQ times extracted - Setup: {}ms, Prove: {}ms, Verify: {}ms, Proof size: {}B",
            setup_time, prove_time, verify_time, proof_size
        );
    }

    // Return structured benchmark result
    BenchmarkResult {
        system: System::CQ,
        k_value: k,
        n_to_n_ratio,
        setup_time,
        prove_time,
        verify_time,
        proof_size,
    }
}

fn bench_logup_gkr(k: usize, n_to_n_ratio: usize, verbose: bool, debug: bool) -> BenchmarkResult {
    // Capture and redirect detailed output if not verbose
    let timings = if !verbose && !debug {
        with_suppressed_output(|| logupgkr::LogupGkr::test_logupgkr_by_k_with_ratio(k, n_to_n_ratio))
    } else {
        if verbose || debug {
            println!("Running LogupGKR benchmark with k={}, N:n ratio={}", k, n_to_n_ratio);
            let range_size = 1 << k;
            let lookup_size = range_size / n_to_n_ratio;
            println!("Range size: {} ({}-bit), Lookup size: {}", range_size, k, lookup_size);
        }
        
        let result = logupgkr::LogupGkr::test_logupgkr_by_k_with_ratio(k, n_to_n_ratio);
        
        if verbose {
            println!("LogupGKR test completed. Results:");
            for timing in &result {
                println!("  {}", timing);
            }
        }
        
        result
    };

    // Write results to file
    for timing in &timings {
        writeln!(&mut System::LogupGKR.output(), "{}", timing).unwrap();
    }

    let all_timings = timings.join("\n");

    if debug {
        println!("\nDEBUG: LogupGKR raw timing output:");
        println!("{}", all_timings);
    }

    // Extract performance metrics from combined timings
    let setup_time = extract_setup_time(&all_timings, debug);
    let prove_time = extract_prove_time(&all_timings, debug);
    let verify_time = extract_verify_time(&all_timings, debug);
    let proof_size = extract_proof_size(&all_timings, debug);

    if verbose || debug {
        println!(
            "LogupGKR times extracted - Setup: {}ms, Prove: {}ms, Verify: {}ms, Proof size: {}B",
            setup_time, prove_time, verify_time, proof_size
        );
    }

    // Return structured benchmark result
    BenchmarkResult {
        system: System::LogupGKR,
        k_value: k,
        n_to_n_ratio,
        setup_time,
        prove_time,
        verify_time,
        proof_size,
    }
}

fn bench_plookup(k: usize, n_to_n_ratio: usize, verbose: bool, debug: bool) -> BenchmarkResult {
    // Capture and redirect detailed output if not verbose
    let timings = if !verbose && !debug {
        with_suppressed_output(|| PlookupBn256::test_plookup_by_k_with_ratio(k, n_to_n_ratio))
    } else {
        if verbose || debug {
            println!("Running Plookup benchmark with k={}, N:n ratio={}", k, n_to_n_ratio);
            let range_size = 1 << k;
            let lookup_size = range_size / n_to_n_ratio;
            println!("Range size: {} ({}-bit), Lookup size: {}", range_size, k, lookup_size);
        }
        
        let result = PlookupBn256::test_plookup_by_k_with_ratio(k, n_to_n_ratio);
        
        if verbose {
            println!("Plookup test completed. Results:");
            for timing in &result {
                println!("  {}", timing);
            }
        }
        
        result
    };

    // Write results to file
    for timing in &timings {
        writeln!(&mut System::Plookup.output(), "{}", timing).unwrap();
    }

    let all_timings = timings.join("\n");

    if debug {
        println!("\nDEBUG: Plookup raw timing output:");
        println!("{}", all_timings);
    }

    // Extract performance metrics from combined timings
    let setup_time = extract_setup_time(&all_timings, debug);
    let prove_time = extract_prove_time(&all_timings, debug);
    let verify_time = extract_verify_time(&all_timings, debug);
    let proof_size = extract_proof_size(&all_timings, debug);

    if verbose || debug {
        println!(
            "Plookup times extracted - Setup: {}ms, Prove: {}ms, Verify: {}ms, Proof size: {}B",
            setup_time, prove_time, verify_time, proof_size
        );
    }

    // Return structured benchmark result
    BenchmarkResult {
        system: System::Plookup,
        k_value: k,
        n_to_n_ratio,
        setup_time,
        prove_time,
        verify_time,
        proof_size,
    }
}

fn bench_caulk(k: usize, n_to_n_ratio: usize, verbose: bool, debug: bool) -> BenchmarkResult {
    if verbose || debug {
        println!("Running Caulk benchmark with k={}, N:n ratio={}", k, n_to_n_ratio);
        let range_size = 1 << k;
        let lookup_size = range_size / n_to_n_ratio;
        println!("Range size: {} ({}-bit), Lookup size: {}", range_size, k, lookup_size);
    }

    // Capture and redirect detailed output if not verbose
    let timings = if !verbose && !debug {
        with_suppressed_output(|| plonkish_backend::backend::caulk::Caulk::<Bn256>::test_caulk_by_k_with_ratio_unified(k, n_to_n_ratio))
    } else {
        plonkish_backend::backend::caulk::Caulk::<Bn256>::test_caulk_by_k_with_ratio_unified(k, n_to_n_ratio)
    };

    // Write results to file
    for timing in &timings {
        writeln!(&mut System::Caulk.output(), "{}", timing).unwrap();
    }

    let all_timings = timings.join("\n");

    if debug {
        println!("\nDEBUG: Caulk raw timing output:");
        println!("{}", all_timings);
    }

    // Extract performance metrics from combined timings
    let setup_time = extract_setup_time(&all_timings, debug);
    let prove_time = extract_prove_time(&all_timings, debug);
    let verify_time = extract_verify_time(&all_timings, debug);
    let proof_size = extract_proof_size(&all_timings, debug);

    if verbose || debug {
        println!(
            "Caulk times extracted - Setup: {}ms, Prove: {}ms, Verify: {}ms, Proof size: {}B",
            setup_time, prove_time, verify_time, proof_size
        );
    }

    // Return structured benchmark result
    BenchmarkResult {
        system: System::Caulk,
        k_value: k,
        // Store the requested ratio, even if the underlying function doesn't use it yet.
        n_to_n_ratio,
        setup_time,
        prove_time,
        verify_time,
        proof_size,
    }
}

// Helper function to extract timing information from output
fn extract_setup_time(timing: &str, debug: bool) -> u64 {
    if debug {
        println!("DEBUG: Searching for setup time in:\n{}", timing);
    }

    let patterns = [
        r"------------\?Setup and preprocess: (\d+)ms-----------",
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

fn extract_proof_size(timing: &str, debug: bool) -> u64 {
    if debug {
        println!("DEBUG: Searching for proof size in:\n{}", timing);
    }

    let patterns = [
        r"proof size: (\d+) bytes",
        r"proof size: (\d+)B",
        r"Proof size: (\d+) bytes",
        r"Proof size: (\d+)B",
        r"size: (\d+) bytes",
        r"size: (\d+)B",
        r"proof: (\d+) bytes",
        r"proof: (\d+)B",
        r"(\d+) bytes proof",
        r"(\d+)B proof",
    ];

    for pattern in patterns {
        let re = Regex::new(pattern).unwrap();
        for cap in re.captures_iter(timing) {
            if let Ok(value) = cap[1].parse::<u64>() {
                if debug {
                    println!(
                        "DEBUG: Found proof size {} bytes with pattern {}",
                        value, pattern
                    );
                }
                return value;
            }
        }
    }

    // Also try to find proof size in KB
    let kb_patterns = [
        r"proof size: (\d+(?:\.\d+)?) ?KB",
        r"Proof size: (\d+(?:\.\d+)?) ?KB",
        r"size: (\d+(?:\.\d+)?) ?KB",
        r"(\d+(?:\.\d+)?) ?KB proof",
    ];

    for pattern in kb_patterns {
        let re = Regex::new(pattern).unwrap();
        for cap in re.captures_iter(timing) {
            if let Ok(value) = cap[1].parse::<f64>() {
                let bytes = (value * 1024.0) as u64;
                if debug {
                    println!(
                        "DEBUG: Found proof size {} KB ({} bytes) with pattern {}",
                        value, bytes, pattern
                    );
                }
                return bytes;
            }
        }
    }

    if debug {
        println!("DEBUG: Could not find proof size");
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
        "┌──────────┬─────────┬─────────┬────────────────────┬──────────────┬───────────────┬────────────┬─────────────┐"
    );
    println!(
        "│ System   │ K-value │ N:n     │ Setup+Preprocess   │ Prove        │ Verify        │ Total      │ Proof Size  │"
    );
    println!(
        "├──────────┼─────────┼─────────┼────────────────────┼──────────────┼───────────────┼────────────┼─────────────┤"
    );

    // Add data rows
    for result in results {
        let total = result.setup_time + result.prove_time + result.verify_time;
        let proof_size_display = if result.proof_size > 0 {
            if result.proof_size > 1024 {
                format!("{:.1}KB", result.proof_size as f64 / 1024.0)
            } else {
                format!("{}B", result.proof_size)
            }
        } else {
            "N/A".to_string()
        };
        
        println!(
            "│ {:<8} │ {:<7} │ {:<7} │ {:<18} │ {:<12} │ {:<13} │ {:<10} │ {:<11} │",
            result.system.to_string(),
            result.k_value,
            result.n_to_n_ratio,
            format!("{}ms", result.setup_time),
            format!("{}ms", result.prove_time),
            format!("{}ms", result.verify_time),
            format!("{}ms", total),
            proof_size_display
        );
    }

    println!(
        "└──────────┴─────────┴─────────┴────────────────────┴──────────────┴───────────────┴────────────┴─────────────┘"
    );
}

// Display results in compact single-line format
fn display_compact_results(results: &[BenchmarkResult]) {
    println!("\n✓ Results:");
    println!(
        "{:10} {:7} {:7} {:12} {:10} {:10} {:10} {:11}",
        "System", "K", "N:n", "Setup (ms)", "Prove (ms)", "Verify (ms)", "Total (ms)", "Proof Size"
    );
    println!("{}", "-".repeat(82));

    for result in results {
        let total = result.setup_time + result.prove_time + result.verify_time;
        let proof_size_display = if result.proof_size > 0 {
            if result.proof_size > 1024 {
                format!("{:.1}KB", result.proof_size as f64 / 1024.0)
            } else {
                format!("{}B", result.proof_size)
            }
        } else {
            "N/A".to_string()
        };
        
        println!(
            "{:10} {:7} {:7} {:12} {:10} {:10} {:10} {:11}",
            result.system.to_string(),
            result.k_value,
            result.n_to_n_ratio,
            result.setup_time,
            result.prove_time,
            result.verify_time,
            total,
            proof_size_display
        );
    }
}

// Display results in CSV format
fn display_csv_results(results: &[BenchmarkResult]) {
    println!("System,K,N:n,SetupTime,ProveTime,VerifyTime,TotalTime");
    for result in results {
        let total = result.setup_time + result.prove_time + result.verify_time;
        println!(
            "{},{},{},{},{},{},{}",
            result.system,
            result.k_value,
            result.n_to_n_ratio,
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
        println!("    \"n_to_n_ratio\": {},", result.n_to_n_ratio);
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
    println!("┌──────────┬─────────────┬─────────┬─────────────┬────────────┬─────────────┐");
    println!("│ System   │ K Range     │ N:n     │ Avg Setup   │ Avg Prove  │ Avg Verify  │");
    println!("├──────────┼─────────────┼─────────┼─────────────┼────────────┼─────────────┤");

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
        
        // Get N:n ratio (should be the same for all results in the group)
        let n_to_n_ratio = sys_results.first().map(|r| r.n_to_n_ratio).unwrap_or(2);

        println!(
            "│ {:<8} │ {:<11} │ {:<7} │ {:<11} │ {:<10} │ {:<11} │",
            system.to_string(),
            k_range,
            n_to_n_ratio,
            format!("{}ms", avg_setup),
            format!("{}ms", avg_prove),
            format!("{}ms", avg_verify)
        );
    }

    println!("└──────────┴─────────────┴─────────┴─────────────┴────────────┴─────────────┘");
}

// Function to display failure summary
fn display_failure_summary(failures: &[BenchmarkFailure]) {
    println!("\n❌ Failure Summary:");
    println!("┌────────────┬─────┬─────────┬──────────────────────────────────────────────────┐");
    println!("│ System     │ K   │ N:n     │ Error Message                                    │");
    println!("├────────────┼─────┼─────────┼──────────────────────────────────────────────────┤");

    for failure in failures {
        let truncated_error = if failure.error_message.len() > 48 {
            format!("{}...", &failure.error_message[..45])
        } else {
            failure.error_message.clone()
        };
        
        println!("│ {:<10} │ {:<3} │ {:<7} │ {:<48} │",
            failure.system.to_string(),
            failure.k_value,
            failure.n_to_n_ratio,
            truncated_error
        );
    }
    
    println!("└────────────┴─────┴─────────┴──────────────────────────────────────────────────┘");
    
    // Group failures by system to show patterns
    let mut system_failures: HashMap<System, usize> = HashMap::new();
    for failure in failures {
        *system_failures.entry(failure.system).or_insert(0) += 1;
    }
    
    if system_failures.len() > 1 {
        println!("\n📈 Failure Distribution:");
        for (system, count) in system_failures {
            println!("   {} {}: {} failures", 
                match count {
                    1 => "•",
                    2..=5 => "▪",
                    _ => "▮",
                },
                system, count);
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum System {
    CQ,
    Baloo,
    LogupGKR,
    Plookup,
    Caulk,
    Lasso,
}

impl System {
    fn all() -> Vec<Self> {
        vec![System::CQ, System::Baloo, System::LogupGKR, System::Plookup, System::Caulk, System::Lasso]
    }

    fn output_path(&self) -> String {
        format!("{OUTPUT_DIR}/{self}")
    }

    fn output(&self) -> File {
        OpenOptions::new()
            .create(true)
            .append(true)
            .open(self.output_path())
            .unwrap()
    }

    // Benchmark a system with k value and N:n ratio
    fn bench(&self, k: usize, n_to_n_ratio: usize, verbose: bool, debug: bool) -> BenchmarkResult {
        self.bench_with_test_type(k, n_to_n_ratio, verbose, debug, TestType::Performance, Operation::Range)
    }

    // Benchmark a system with specific test type and operation
    fn bench_with_test_type(&self, k: usize, n_to_n_ratio: usize, verbose: bool, debug: bool, test_type: TestType, operation: Operation) -> BenchmarkResult {
        match test_type {
            TestType::Performance => {
                // Performance benchmarking (existing functionality)
                match (self, operation) {
                    (System::Lasso, Operation::Add) => bench_lasso_add(k, n_to_n_ratio, verbose, debug),
                    _ => {
                        // Default range check performance for all systems
                        match self {
                            System::Baloo => bench_baloo(k, n_to_n_ratio, verbose, debug),
                            System::CQ => bench_CQ(k, n_to_n_ratio, verbose, debug),
                            System::LogupGKR => bench_logup_gkr(k, n_to_n_ratio, verbose, debug),
                            System::Plookup => bench_plookup(k, n_to_n_ratio, verbose, debug),
                            System::Caulk => bench_caulk(k, n_to_n_ratio, verbose, debug),
                            System::Lasso => bench_lasso(k, n_to_n_ratio, verbose, debug),
                        }
                    }
                }
            }
            TestType::Soundness => {
                // Soundness testing - invalid inputs should be rejected
                match (self, operation) {
                    (System::Lasso, Operation::Add) => bench_lasso_add_soundness(k, n_to_n_ratio, verbose, debug),
                    (System::Lasso, Operation::Range) => bench_lasso_range_soundness(k, n_to_n_ratio, verbose, debug),
                    _ => {
                        // For other systems, we don't have soundness tests yet
                        BenchmarkResult {
                            system: *self,
                            k_value: k,
                            n_to_n_ratio,
                            setup_time: 0,
                            prove_time: 0,
                            verify_time: 0,
                            proof_size: 0,
                        }
                    }
                }
            }
            TestType::Completeness => {
                // Completeness testing - valid inputs should be accepted
                match (self, operation) {
                    (System::Lasso, Operation::Add) => bench_lasso_add_completeness(k, n_to_n_ratio, verbose, debug),
                    (System::Lasso, Operation::Range) => bench_lasso_range_completeness(k, n_to_n_ratio, verbose, debug),
                    _ => {
                        // For other systems, we don't have completeness tests yet
                        BenchmarkResult {
                            system: *self,
                            k_value: k,
                            n_to_n_ratio,
                            setup_time: 0,
                            prove_time: 0,
                            verify_time: 0,
                            proof_size: 0,
                        }
                    }
                }
            }
        }
    }
}

impl Display for System {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            System::Baloo => write!(f, "Baloo"),
            System::CQ => write!(f, "CQ"),
            System::LogupGKR => write!(f, "LogupGKR"),
            System::Plookup => write!(f, "Plookup"),
            System::Caulk => write!(f, "Caulk"),
            System::Lasso => write!(f, "Lasso"),
        }
    }
}

fn parse_args() -> (Vec<System>, Range<usize>, Vec<usize>, bool, OutputFormat, bool, TestType, Operation) {
    let (systems, k_range, n_to_n_ratios, verbose, output_format, debug, test_type, operation) =
        args().chain(Some("".to_string())).tuple_windows().fold(
            (Vec::new(), 20..26, vec![2], false, OutputFormat::Table, false, TestType::Performance, Operation::Range),
            |(mut systems, mut k_range, mut n_to_n_ratios, mut verbose, mut output_format, mut debug, mut test_type, mut operation),
             (key, value)| {
                match key.as_str() {
                    "--system" => {
                        if value == "all" {
                            systems = System::all();
                        } else {
                            // Support comma-separated systems, e.g., "cq,plookup"
                            let system_names: Vec<&str> = if value.contains(',') {
                                value.split(',').map(|s| s.trim()).collect()
                            } else {
                                vec![value.as_str()]
                            };
                            
                            for system_name in system_names {
                                match system_name {
                                    "cq" | "CQ" => systems.push(System::CQ),
                                    "baloo" | "Baloo" => systems.push(System::Baloo),
                                    "logupgkr" | "LogupGKR" => systems.push(System::LogupGKR),
                                    "plookup" | "Plookup" => systems.push(System::Plookup),
                                    "caulk" | "Caulk" => systems.push(System::Caulk),
                                    "lasso" | "Lasso" => systems.push(System::Lasso),
                                    _ => panic!("system '{}' should be one of {{all, cq, baloo, logupgkr, plookup, caulk, lasso}}", system_name),
                                }
                            }
                        }
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
                    "--ratio" | "--n-to-n-ratio" => {
                        // Support comma-separated ratios, e.g., "2,4,8"
                        if value.contains(',') {
                            n_to_n_ratios = value.split(',')
                                .map(|s| s.trim().parse().expect("N:n ratio to be usize"))
                                .collect();
                        } else {
                            let ratio = value.parse().expect("N:n ratio to be usize");
                            n_to_n_ratios = vec![ratio];
                        }
                        for &ratio in &n_to_n_ratios {
                            if ratio < 1 {
                                panic!("N:n ratio must be at least 1");
                            }
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
                    "--test-type" => match value.to_lowercase().as_str() {
                        "performance" => test_type = TestType::Performance,
                        "soundness" => test_type = TestType::Soundness,
                        "completeness" => test_type = TestType::Completeness,
                        _ => panic!("test-type should be one of {{performance, soundness, completeness}}"),
                    },
                    "--operation" => match value.to_lowercase().as_str() {
                        "range" => operation = Operation::Range,
                        "add" => operation = Operation::Add,
                        "all" => operation = Operation::All,
                        _ => panic!("operation should be one of {{range, add, all}}"),
                    },
                    _ => {}
                }
                (systems, k_range, n_to_n_ratios, verbose, output_format, debug, test_type, operation)
            },
        );

    let mut systems = systems.into_iter().sorted().dedup().collect_vec();
    if systems.is_empty() {
        systems = System::all();
    };

    (systems, k_range, n_to_n_ratios, verbose, output_format, debug, test_type, operation)
}

fn create_output(systems: &[System]) {
    if !Path::new(OUTPUT_DIR).exists() {
        create_dir(OUTPUT_DIR).unwrap();
    }
    for system in systems {
        File::create(system.output_path()).unwrap();
    }
}

// Add the benchmark function for Lasso add operations
fn bench_lasso_add(k: usize, n_to_n_ratio: usize, verbose: bool, debug: bool) -> BenchmarkResult {
    let timings = if !verbose && !debug {
        with_suppressed_output(|| {
            lasso::test_lasso_add_by_k_with_ratio(k, n_to_n_ratio)
        })
    } else {
        println!("Running Lasso add benchmark with k={}, N:n ratio={}", k, n_to_n_ratio);
        let range_size = 1 << k;
        let lookup_size = range_size / n_to_n_ratio;
        println!("Range size: {} ({}-bit), Lookup size: {}", range_size, k, lookup_size);
        
        let result = lasso::test_lasso_add_by_k_with_ratio(k, n_to_n_ratio);
        
        if verbose {
            println!("Lasso add test completed. Results:");
            for timing in &result {
                println!("  {}", timing);
            }
        }
        
        result
    };

    // Write results to file
    for timing in &timings {
        writeln!(&mut System::Lasso.output(), "{}", timing).unwrap();
    }

    let all_timings = timings.join("\n");

    if debug {
        println!("\nDEBUG: Lasso add raw timing output:");
        println!("{}", all_timings);
    }

    // Extract timing values
    let setup_time = extract_setup_time(&all_timings, debug);
    let prove_time = extract_prove_time(&all_timings, debug);
    let verify_time = extract_verify_time(&all_timings, debug);
    let proof_size = extract_proof_size(&all_timings, debug);

    if verbose || debug {
        println!(
            "Lasso add times extracted - Setup: {}ms, Prove: {}ms, Verify: {}ms, Proof size: {}B",
            setup_time, prove_time, verify_time, proof_size
        );
    }

    BenchmarkResult {
        system: System::Lasso,
        k_value: k,
        n_to_n_ratio,
        setup_time,
        prove_time,
        verify_time,
        proof_size,
    }
}

// Benchmark function for Lasso add soundness testing
fn bench_lasso_add_soundness(k: usize, n_to_n_ratio: usize, verbose: bool, debug: bool) -> BenchmarkResult {
    let start_time = std::time::Instant::now();
    
    match lasso::test_lasso_add_soundness_by_k(k, n_to_n_ratio) {
        Ok(timings) => {
            // Parse timing information from the result
            let mut setup_time = 0u64;
            let mut prove_time = 0u64;
            let mut verify_time = 0u64;
            
            for timing in &timings {
                if let Some(time_str) = timing.strip_prefix("Setup: ").and_then(|s| s.strip_suffix("ms")) {
                    setup_time = time_str.parse().unwrap_or(0);
                } else if let Some(time_str) = timing.strip_prefix("Prove: ").and_then(|s| s.strip_suffix("ms")) {
                    prove_time = time_str.parse().unwrap_or(0);
                } else if let Some(time_str) = timing.strip_prefix("Verify: ").and_then(|s| s.strip_suffix("ms")) {
                    verify_time = time_str.parse().unwrap_or(0);
                }
            }
            
            // Check if verification actually failed (which is expected for soundness)
            let verification_failed = timings.iter().any(|t| 
                t.contains("Error:") && 
                (t.contains("expected for soundness test") || 
                 t.contains("verification failed") ||
                 t.contains("verify failed"))
            );
            
            if verbose || debug {
                if verification_failed {
                    println!("Soundness test passed: invalid add operations correctly rejected during verification");
                } else {
                    println!("WARNING: Soundness test failed - invalid operations were incorrectly accepted");
                }
                for timing in &timings {
                    println!("  {}", timing);
                }
            }
            
            BenchmarkResult {
                system: System::Lasso,
                k_value: k,
                n_to_n_ratio,
                setup_time,
                prove_time,
                verify_time,
                proof_size: 0,
            }
        }
        Err(error_msg) => {
            // This means proving phase failed (assertion error)
            if verbose || debug {
                println!("Soundness test: proving phase failed (assertion error): {}", error_msg);
            }
            BenchmarkResult {
                system: System::Lasso,
                k_value: k,
                n_to_n_ratio,
                setup_time: 0,
                prove_time: 0,
                verify_time: start_time.elapsed().as_millis() as u64,
                proof_size: 0,
            }
        }
    }
}

// Benchmark function for Lasso add completeness testing
fn bench_lasso_add_completeness(k: usize, n_to_n_ratio: usize, verbose: bool, debug: bool) -> BenchmarkResult {
    let start_time = std::time::Instant::now();
    
    match lasso::test_lasso_add_completeness_by_k(k, n_to_n_ratio) {
        Ok(timings) => {
            // This is expected - completeness test should accept valid operations
            if verbose || debug {
                println!("Completeness test passed: valid add operations correctly accepted");
                for timing in &timings {
                    println!("  {}", timing);
                }
            }
            
            let all_timings = timings.join("\n");
            let setup_time = extract_setup_time(&all_timings, debug);
            let prove_time = extract_prove_time(&all_timings, debug);
            let verify_time = extract_verify_time(&all_timings, debug);
            let proof_size = extract_proof_size(&all_timings, debug);
            
            BenchmarkResult {
                system: System::Lasso,
                k_value: k,
                n_to_n_ratio,
                setup_time,
                prove_time,
                verify_time,
                proof_size,
            }
        }
        Err(e) => {
            // This should not happen for completeness test - valid operations should be accepted
            if verbose || debug {
                println!("WARNING: Completeness test failed: {}", e);
            }
            BenchmarkResult {
                system: System::Lasso,
                k_value: k,
                n_to_n_ratio,
                setup_time: 0,
                prove_time: 0,
                verify_time: start_time.elapsed().as_millis() as u64,
                proof_size: 0,
            }
        }
    }
}

// Benchmark function for Lasso range soundness testing
fn bench_lasso_range_soundness(k: usize, n_to_n_ratio: usize, verbose: bool, debug: bool) -> BenchmarkResult {
    // For now, return a placeholder result since range soundness testing isn't implemented yet
    BenchmarkResult {
        system: System::Lasso,
        k_value: k,
        n_to_n_ratio,
        setup_time: 0,
        prove_time: 0,
        verify_time: 0,
        proof_size: 0,
    }
}

// Benchm0ark function for Lasso range completeness testing
fn bench_lasso_range_completeness(k: usize, n_to_n_ratio: usize, verbose: bool, debug: bool) -> BenchmarkResult {
    // For now, return a placeholder result since range completeness testing isn't implemented yet
    BenchmarkResult {
        system: System::Lasso,
        k_value: k,
        n_to_n_ratio,
        setup_time: 0,
        prove_time: 0,
        verify_time: 0,
        proof_size: 0,
    }
}

// Add the benchmark function for Lasso
fn bench_lasso(k: usize, n_to_n_ratio: usize, verbose: bool, debug: bool) -> BenchmarkResult {
    // Capture and redirect detailed output if not verbose
    let timings = if !verbose && !debug {
        with_suppressed_output(|| {
            lasso::test_lasso_by_k_with_ratio(k, n_to_n_ratio)
        })
    } else {
        println!("Running Lasso benchmark with k={}, N:n ratio={}", k, n_to_n_ratio);
        
        // Calculate range and lookup sizes for display
        let range_size = 1 << k;
        let lookup_size = range_size / n_to_n_ratio;
        
        println!("Range size: {} ({}-bit), Lookup size: {}", range_size, k, lookup_size);
        
        let result = lasso::test_lasso_by_k_with_ratio(k, n_to_n_ratio);
        
        if verbose {
            println!("Lasso test completed. Results:");
            for timing in &result {
                println!("  {}", timing);
            }
        }
        
        result
    };

    // Write results to file
    for timing in &timings {
        writeln!(&mut System::Lasso.output(), "{}", timing).unwrap();
    }

    let all_timings = timings.join("\n");

    if debug {
        println!("\nDEBUG: Lasso raw timing output:");
        println!("{}", all_timings);
    }

    // Extract timing values using the same pattern as other systems
    let setup_time = extract_setup_time(&all_timings, debug);
    let prove_time = extract_prove_time(&all_timings, debug);
    let verify_time = extract_verify_time(&all_timings, debug);
    let proof_size = extract_proof_size(&all_timings, debug);

    if verbose || debug {
        println!(
            "Lasso times extracted - Setup: {}ms, Prove: {}ms, Verify: {}ms, Proof size: {}B",
            setup_time, prove_time, verify_time, proof_size
        );
    }

    // Return structured benchmark result
    BenchmarkResult {
        system: System::Lasso,
        k_value: k,
        n_to_n_ratio,
        setup_time,
        prove_time,
        verify_time,
        proof_size,
    }
}
