use itertools::Itertools;
use plonkish_backend::backend::{cq, logupgkr, plookup, lasso, caulk};
use plonkish_backend::halo2_curves::bn256::{Bn256, Fr};
use plonkish_backend::pcs::univariate::UnivariateKzg;
use rayon::prelude::*;
use regex::Regex;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::{
    env::args,
    fmt::Display,
    fs::{create_dir, File, OpenOptions},
    io::Write,
    ops::Range,
    path::Path,
};

// Type alias for Plookup with BN256 curve
type PlookupBn256 = plookup::Plookup<Fr, UnivariateKzg<Bn256>>;

const OUTPUT_DIR: &str = "../target/bench";

fn main() {
    // Parse arguments
    let (systems, k_range, n_to_n_ratios, verbose, output_format, debug) = parse_args();
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
    println!("{}", "=".repeat(60));

    // Use Arc<Mutex<Vec<BenchmarkResult>>> for thread-safe result collection
    let all_results = Arc::new(Mutex::new(Vec::new()));
    let completed_counter = Arc::new(Mutex::new(0));

    // Run benchmarks in parallel using rayon
    benchmark_tasks.par_iter().enumerate().for_each(|(task_index, (system, k, ratio))| {
        // Show progress before starting
        {
            let mut counter = completed_counter.lock().unwrap();
            println!("🔄 [{:3}/{:3}] Starting: {} (k={}, ratio={})", 
                    *counter + 1, total_tasks, system.to_string(), k, ratio);
        }

        let result = system.bench(*k, *ratio, verbose, debug);
        
        // Safely add result to shared collection and update counter
        {
            let mut counter = completed_counter.lock().unwrap();
            all_results.lock().unwrap().push(result);
            *counter += 1;
            println!("✅ [{:3}/{:3}] Completed: {} (k={}, ratio={}) - {:.1}% done", 
                    *counter, total_tasks, system.to_string(), k, ratio, 
                    (*counter as f64 / total_tasks as f64) * 100.0);
        }
    });

    // Extract results from Arc<Mutex<>>
    let final_results = all_results.lock().unwrap().clone();

    // Sort results for consistent display order (by system, then k, then ratio)
    let mut sorted_results = final_results;
    sorted_results.sort_by(|a, b| {
        a.system.to_string()
            .cmp(&b.system.to_string())
            .then(a.k_value.cmp(&b.k_value))
            .then(a.n_to_n_ratio.cmp(&b.n_to_n_ratio))
    });

    println!("\n🎉 All {} benchmark tasks completed successfully!", sorted_results.len());
    println!("{}", "=".repeat(60));

    // Display results based on selected format
    match output_format {
        OutputFormat::Table => display_table_results(&sorted_results),
        OutputFormat::Compact => display_compact_results(&sorted_results),
        OutputFormat::CSV => display_csv_results(&sorted_results),
        OutputFormat::JSON => display_json_results(&sorted_results),
    }

    // Display summary statistics
    display_summary(&sorted_results);
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
    // Create test data based on k
    let table_size = 1 << k.min(10); // 2^k, 最大限制为 2^10 以防止过大
    
    // Calculate lookup size based on the ratio N:n
    let lookup_size = table_size / n_to_n_ratio;
    
    // Create lookup and table vectors
    let table: Vec<Fr> = (1..=table_size).map(|i| Fr::from(i as u64)).collect();
    
    // Create lookup with some repeated values for testing multiplicities
    let mut lookup = Vec::with_capacity(lookup_size);
    for i in 0..lookup_size {
        // Add some repetition pattern - use modulo to create repeating values
        let value = i % (table_size / 4).max(1) + 1;
        lookup.push(Fr::from(value as u64));
    }
    
    // Capture and redirect detailed output if not verbose
    let timings = if !verbose && !debug {
        with_suppressed_output(|| {
            // Convert to logupgkr format using the utility function
            let (m_values, t_values, w_values) = logupgkr::util::convert_to_logupgkr_format(lookup, table);
            
            // Calculate 'a' parameter (can be any value for testing)
            let a = Fr::from(table_size as u64 + 1);
            
            logupgkr::LogupGkr::test_logupgkr(
                m_values,
                t_values,
                w_values,
            )
        })
    } else {
        println!("Running LogupGKR benchmark with k={}, N:n ratio={}", k, n_to_n_ratio);
        println!("Table size: {}, Lookup size: {}", table_size, lookup_size);
        
        // Convert to logupgkr format using the utility function
        let (m_values, t_values, w_values) = logupgkr::util::convert_to_logupgkr_format(lookup, table);
        
        // Calculate 'a' parameter (can be any value for testing)
        let a = Fr::from(table_size as u64 + 1);
        
        println!("Running LogupGKR test with converted data...");
        let result = logupgkr::LogupGkr::test_logupgkr(
            m_values,
            t_values,
            w_values,
        );
        
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
    // Calculate table and lookup sizes based on k and ratio
    let table_size = 1 << k;
    let lookup_size = table_size / n_to_n_ratio;
    
    // Generate table and lookup values using the utility function from CQ
    let (table, lookup) = cq::generate_table_and_lookup(table_size, lookup_size);
    
    // Capture and redirect detailed output if not verbose
    let timings = if !verbose && !debug {
        with_suppressed_output(|| {
            PlookupBn256::test_plookup_by_input(
                table.clone(),
                lookup.clone(),
            )
        })
    } else {
        println!("Running Plookup benchmark with k={}, N:n ratio={}", k, n_to_n_ratio);
        println!("Table size: {}, Lookup size: {}", table_size, lookup_size);
        
        let result = PlookupBn256::test_plookup_by_input(
            table.clone(),
            lookup.clone(),
        );
        
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
        let table_size = 1 << k;
        let lookup_size = table_size / n_to_n_ratio;
        println!("Table size: {}, Lookup size: {}", table_size, lookup_size);
    }

    // Capture and redirect detailed output if not verbose
    let timings = if !verbose && !debug {
        with_suppressed_output(|| plonkish_backend::backend::caulk::Caulk::<Bn256>::test_caulk_by_k_and_ratio(k, n_to_n_ratio))
    } else {
        plonkish_backend::backend::caulk::Caulk::<Bn256>::test_caulk_by_k_and_ratio(k, n_to_n_ratio)
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

fn parse_args() -> (Vec<System>, Range<usize>, Vec<usize>, bool, OutputFormat, bool) {
    let (systems, k_range, n_to_n_ratios, verbose, output_format, debug) =
        args().chain(Some("".to_string())).tuple_windows().fold(
            (Vec::new(), 20..26, vec![2], false, OutputFormat::Table, false),
            |(mut systems, mut k_range, mut n_to_n_ratios, mut verbose, mut output_format, mut debug),
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
                    _ => {}
                }
                (systems, k_range, n_to_n_ratios, verbose, output_format, debug)
            },
        );

    let mut systems = systems.into_iter().sorted().dedup().collect_vec();
    if systems.is_empty() {
        systems = System::all();
    };

    (systems, k_range, n_to_n_ratios, verbose, output_format, debug)
}

fn create_output(systems: &[System]) {
    if !Path::new(OUTPUT_DIR).exists() {
        create_dir(OUTPUT_DIR).unwrap();
    }
    for system in systems {
        File::create(system.output_path()).unwrap();
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
        let range_bits = k.min(8);
        let range_size = 1 << range_bits;
        let lookup_size = range_size / n_to_n_ratio;
        
        println!("Range size: {} ({}-bit), Lookup size: {}", range_size, range_bits, lookup_size);
        
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
