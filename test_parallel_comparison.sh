#!/bin/bash

# Test script for comparing parallel vs non-parallel performance in Lasso
# This script runs benchmarks with and without the parallel feature enabled

echo "🧪 Testing Lasso Parallel vs Non-Parallel Performance"
echo "======================================================"

# Configuration
K_VALUES="8..10"
RATIOS="2,4"
SYSTEM="lasso"

echo "📋 Test Configuration:"
echo "  - System: $SYSTEM"
echo "  - K values: $K_VALUES"
echo "  - N:n ratios: $RATIOS"
echo ""

# Test 1: With parallel feature enabled (default)
echo "🚀 Test 1: Running with parallel algorithms ENABLED"
echo "Command: cargo bench --bench proof_system -- --system $SYSTEM --k $K_VALUES --ratio $RATIOS --format compact"
echo ""

time cargo bench --bench proof_system -- \
    --system $SYSTEM \
    --k $K_VALUES \
    --ratio $RATIOS \
    --format compact \
    --verbose

echo ""
echo "=============================================="
echo ""

# Test 2: With parallel feature disabled
echo "🐌 Test 2: Running with parallel algorithms DISABLED"
echo "Command: cargo bench --bench proof_system --no-default-features --features timer -- --system $SYSTEM --k $K_VALUES --ratio $RATIOS --format compact"
echo ""

time cargo bench --bench proof_system \
    --no-default-features \
    --features timer \
    -- \
    --system $SYSTEM \
    --k $K_VALUES \
    --ratio $RATIOS \
    --format compact \
    --verbose

echo ""
echo "✅ Parallel vs Non-Parallel comparison completed!"
echo ""
echo "💡 Tips:"
echo "  - Compare the 'Prove' times between the two runs"
echo "  - Parallel version should be faster on multi-core systems"
echo "  - Check CPU utilization during the runs using 'htop' or 'top'"
echo ""
echo "📁 Detailed results are saved in ../target/bench/Lasso" 