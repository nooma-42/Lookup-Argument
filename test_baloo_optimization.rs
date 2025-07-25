#!/usr/bin/env -S cargo +nightly -Zscript

// Test script to compare legacy vs optimized Baloo implementation performance
//! ```cargo
//! [dependencies]
//! plonkish_backend = { path = "plonkish_backend" }
//! ```

use plonkish_backend::backend::baloo::Baloo;

fn main() {
    println!("🚀 Testing Baloo Implementation Performance Optimization");
    println!("=========================================================");
    
    for k in [8, 10, 12] {
        println!("\n📊 Testing with K = {} (table size = {}, lookup size = {})", k, 2_usize.pow(k), 2_usize.pow(k-1));
        
        // Test legacy implementation
        println!("\n🔄 Legacy Implementation:");
        let legacy_timings = Baloo::test_baloo_by_k(k);
        for timing in &legacy_timings {
            println!("  {}", timing);
        }
        
        // Test optimized implementation  
        println!("\n⚡ Optimized Implementation:");
        let optimized_timings = Baloo::test_baloo_optimized_by_k(k);
        for timing in &optimized_timings {
            println!("  {}", timing);
        }
        
        // Extract and compare proving times
        let legacy_prove_time = extract_prove_time(&legacy_timings);
        let optimized_prove_time = extract_prove_time(&optimized_timings);
        
        if let (Some(legacy), Some(optimized)) = (legacy_prove_time, optimized_prove_time) {
            let speedup = legacy as f64 / optimized as f64;
            println!("\n💨 Performance Improvement:");
            println!("  Legacy prove time: {}ms", legacy);
            println!("  Optimized prove time: {}ms", optimized);
            println!("  Speedup: {:.2}x", speedup);
            
            if speedup > 1.0 {
                println!("  ✅ Optimization successful!");
            } else {
                println!("  ⚠️  No significant improvement");
            }
        }
        
        println!("\n" + &"=".repeat(60));
    }
}

fn extract_prove_time(timings: &[String]) -> Option<u64> {
    for timing in timings {
        if timing.contains("Prove:") || timing.contains("Optimized Prove:") {
            // Extract number from string like "Prove: 123ms" or "Optimized Prove: 123ms"
            let parts: Vec<&str> = timing.split(':').collect();
            if parts.len() >= 2 {
                let time_part = parts[1].trim().replace("ms", "");
                return time_part.parse().ok();
            }
        }
    }
    None
}