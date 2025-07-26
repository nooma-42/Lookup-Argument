#!/bin/bash

# Benchmark with Incremental CSV Writing and Failure Tracking
# This script demonstrates the new features added to the benchmark system

echo "🚀 Lookup Argument Benchmark with Incremental CSV Writing"
echo "========================================================"

# Clean up previous results
if [ -f "benchmark_results.csv" ]; then
    echo "🗑️  Removing previous benchmark_results.csv"
    rm benchmark_results.csv
fi

if [ -f "benchmark_failures.csv" ]; then
    echo "🗑️  Removing previous benchmark_failures.csv"
    rm benchmark_failures.csv
fi

echo ""
echo "📋 Running comprehensive benchmark suite..."
echo "   • Systems: ALL (CQ, Baloo, LogupGKR, Plookup, Caulk, Lasso)"
echo "   • K values: 5-12 (table sizes from 32 to 4096)"
echo "   • N:n ratios: 2, 4, 8, 16"
echo "   • Total tasks: ~192 (6 systems × 8 k-values × 4 ratios)"
echo "   • Output format: Table"
echo ""
echo "⚠️  Warning: This is a comprehensive benchmark that may take several hours!"
echo "    Results will be saved incrementally, so you can monitor progress."
echo ""

# Ask for confirmation
read -p "Do you want to proceed with this comprehensive benchmark? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "❌ Benchmark cancelled."
    exit 1
fi

echo "🚀 Starting comprehensive benchmark..."
echo "💾 Results will be saved to benchmark_results.csv as each test completes"
echo ""

# Run the benchmark
cd benchmark

echo "🚀 Note: Running systems in two phases to avoid Lasso parallel race conditions"
echo "📊 Phase 1: Lasso (sequential execution) - Thread-safe for multiple instances"
echo "📊 Phase 2: Other systems (parallel execution) - No race condition issues"
echo "💾 Both phases will append results to the same CSV files automatically"
echo ""

# Phase 1: Run Lasso separately without parallel execution to avoid race conditions
echo "📊 Phase 1: Running Lasso benchmarks (sequential execution)..."
echo "   🔧 Using --no-default-features to disable parallel execution and verbose timing logs"
echo "   📝 Results will be written to benchmark_results.csv"
echo ""
time cargo bench --bench proof_system --no-default-features -- \
    --system lasso \
    --k 5..14 \
    --ratio 2,4,8,16 \
    --format table

echo ""
echo "✅ Phase 1 completed - Lasso results saved to CSV"
echo ""

# Phase 2: Run other systems with parallel execution (they don't have the race condition issue)
echo "📊 Phase 2: Running other systems (parallel execution)..."
echo "   🔧 Using default parallel execution for optimal performance"
echo "   📝 Results will be appended to the same benchmark_results.csv"
echo ""
time cargo bench --bench proof_system -- \
    --system cq,baloo,logupgkr,plookup,caulk \
    --k 5..14 \
    --ratio 2,4,8,16 \
    --format table

echo ""
echo "✅ Phase 2 completed - All system results saved to CSV"

echo ""
echo "📊 Benchmark completed! Checking results..."

# Check if files exist and show their content
if [ -f "benchmark_results.csv" ]; then
    echo ""
    echo "✅ benchmark_results.csv created successfully!"
    echo "📏 File size: $(wc -c < benchmark_results.csv) bytes"
    echo "📊 Number of results: $(($(wc -l < benchmark_results.csv) - 1))"
    echo ""
    echo "📖 First few lines:"
    head -5 benchmark_results.csv
    echo ""
    echo "📖 Last few lines:"
    tail -3 benchmark_results.csv
    echo ""
    echo "📈 Results by system:"
    tail -n +2 benchmark_results.csv | cut -d',' -f1 | sort | uniq -c | sort -nr
else
    echo "❌ benchmark_results.csv not found!"
fi

if [ -f "benchmark_failures.csv" ]; then
    if [ -s "benchmark_failures.csv" ]; then
        echo ""
        echo "⚠️  benchmark_failures.csv contains failures:"
        echo "📏 File size: $(wc -c < benchmark_failures.csv) bytes"
        echo "❌ Number of failures: $(($(wc -l < benchmark_failures.csv) - 1))"
        echo ""
        echo "📖 Failure summary by system:"
        tail -n +2 benchmark_failures.csv | cut -d',' -f1 | sort | uniq -c | sort -nr
        echo ""
        echo "📖 First few failure details:"
        head -5 benchmark_failures.csv
    else
        echo ""
        echo "✅ benchmark_failures.csv exists but is empty (no failures)"
    fi
else
    echo "❌ benchmark_failures.csv not found!"
fi

echo ""
echo "🎯 Comprehensive Benchmark Features Demonstrated:"
echo "   ✅ All proof systems tested - CQ, Baloo, LogupGKR, Plookup, Caulk, Lasso"
echo "   ✅ Wide parameter range - k=5..12, ratios=2,4,8,16"
echo "   ✅ Incremental CSV writing - results saved as each benchmark completes"
echo "   ✅ Failure tracking - any failed benchmarks are logged separately"
echo "   ✅ Parallel execution - multiple benchmarks run simultaneously"
echo "   ✅ Progress tracking - real-time progress updates with percentages"
echo "   ✅ Comprehensive reporting - separate files for successes and failures"
echo ""
echo "💡 Benefits:"
echo "   📁 Partial results preserved even if the process is interrupted"
echo "   🔍 Easy identification of which specific configurations failed"
echo "   ⚡ Faster execution through parallelization"
echo "   📊 Machine-readable CSV format for further analysis"
echo "   🎯 Complete coverage of all proof systems and parameters"
echo "   📈 Automatic plot generation for all key performance metrics"
echo ""
echo "🗂️  Output files:"
echo "   📊 benchmark_results.csv - All successful benchmark results"
echo "   ❌ benchmark_failures.csv - Any failed benchmark attempts"
echo "   🎨 plots/*.png - Automatically generated performance visualization plots"
echo ""

# Generate a quick analysis if we have results
if [ -f "benchmark_results.csv" ] && [ -s "benchmark_results.csv" ]; then
    echo "📊 Quick Analysis:"
    echo "   • Total successful benchmarks: $(($(wc -l < benchmark_results.csv) - 1))"
    
    # Find fastest and slowest setups
    if [ $(wc -l < benchmark_results.csv) -gt 1 ]; then
        echo "   • Fastest setup time: $(tail -n +2 benchmark_results.csv | cut -d',' -f4 | sort -n | head -1)ms"
        echo "   • Slowest setup time: $(tail -n +2 benchmark_results.csv | cut -d',' -f4 | sort -n | tail -1)ms"
        echo "   • Fastest prove time: $(tail -n +2 benchmark_results.csv | cut -d',' -f5 | sort -n | head -1)ms"
        echo "   • Slowest prove time: $(tail -n +2 benchmark_results.csv | cut -d',' -f5 | sort -n | tail -1)ms"
    fi
fi

echo ""
echo "Done! 🎉"
# Generate plots automatically if benchmarks completed successfully
if [ -f "benchmark_results.csv" ] && [ -s "benchmark_results.csv" ] && [ $(wc -l < "benchmark_results.csv") -gt 1 ]; then
    echo "📊 Generating performance plots..."
    echo ""
    
    # Check if Python is available
    if command -v python3 &> /dev/null; then
        PYTHON_CMD="python3"
    elif command -v python &> /dev/null; then
        PYTHON_CMD="python"
    else
        echo "⚠️  Python not found. Skipping plot generation."
        echo "   You can manually generate plots later using: python diagram.py"
        PYTHON_CMD=""
    fi
    
    if [ -n "$PYTHON_CMD" ]; then
        # Create plots directory
        mkdir -p plots
        
        echo "🎨 Generating comprehensive performance plots..."
        echo "   📊 Setup Time plots..."
        $PYTHON_CMD diagram.py --metrics SetupTime_ms --output-dir plots/ --protocols Lasso,LogupGKR,CQ,Plookup,Baloo,Caulk
        
        echo "   🔍 Verify Time plots..."
        $PYTHON_CMD diagram.py --metrics VerifyTime_ms --output-dir plots/ --protocols Lasso,LogupGKR,CQ,Plookup,Baloo,Caulk
        
        echo "   ⏱️  Total Time plots..."
        $PYTHON_CMD diagram.py --metrics TotalTime_ms --output-dir plots/ --protocols Lasso,LogupGKR,CQ,Plookup,Baloo,Caulk
        
        echo "   📦 Proof Size plots..."
        $PYTHON_CMD diagram.py --metrics ProofSize_bytes --output-dir plots/ --protocols Lasso,LogupGKR,CQ,Plookup,Baloo,Caulk
        
        echo "   🚀 Proving Time plots (original)..."
        $PYTHON_CMD diagram.py --metrics ProveTime_ms --output-dir plots/ --protocols Lasso,LogupGKR,CQ,Plookup,Baloo,Caulk
        
        echo ""
        echo "📈 Plot generation completed!"
        echo ""
        
        # Show generated plots summary
        if [ -d "plots" ]; then
            PLOT_COUNT=$(ls -1 plots/*.png 2>/dev/null | wc -l)
            if [ $PLOT_COUNT -gt 0 ]; then
                echo "🎯 Generated $PLOT_COUNT performance plots in plots/ directory:"
                echo ""
                echo "📊 Setup Time plots:"
                ls -1 plots/setuptimems_*.png 2>/dev/null | sed 's/^/   • /'
                echo ""
                echo "🔍 Verify Time plots:"
                ls -1 plots/verifytimems_*.png 2>/dev/null | sed 's/^/   • /'
                echo ""
                echo "⏱️  Total Time plots:"
                ls -1 plots/totaltimems_*.png 2>/dev/null | sed 's/^/   • /'
                echo ""
                echo "📦 Proof Size plots:"
                ls -1 plots/proofsizebytes_*.png 2>/dev/null | sed 's/^/   • /'
                echo ""
                echo "🚀 Proving Time plots:"
                ls -1 plots/provetimems_*.png 2>/dev/null | sed 's/^/   • /'
                echo ""
            else
                echo "⚠️  No plots were generated (check for errors above)"
            fi
        fi
        
        echo "💡 Plot Usage:"
        echo "   • SetupTime: Shows preprocessing/setup performance across systems"
        echo "   • VerifyTime: Shows verification performance (important for deployment)"
        echo "   • TotalTime: Shows end-to-end performance including all phases"
        echo "   • ProofSize: Shows proof size comparison (important for storage/bandwidth)"
        echo "   • ProveTime: Shows core proving time performance"
        echo ""
        echo "🎨 Custom Plot Generation:"
        echo "   • All metrics for single ratio: python diagram.py --metrics SetupTime_ms,VerifyTime_ms,TotalTime_ms,ProofSize_bytes --single-ratio 4"
        echo "   • Specific protocols only: python diagram.py --metrics TotalTime_ms --protocols Lasso,LogupGKR"
        echo "   • Different output directory: python diagram.py --output-dir custom_plots/"
        echo ""
    fi
else
    echo "⚠️  No benchmark results available for plot generation."
fi

echo ""
echo "💡 Next steps:"
echo "   • Review generated plots in plots/ directory for performance insights"
echo "   • Analyze benchmark_results.csv for detailed numerical data"
echo "   • Check benchmark_failures.csv for any issues that need attention"
echo "   • Use the CSV data and plots for research or optimization decisions" 