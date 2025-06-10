use plonkish_backend::{logging, log_info, log_debug};

fn main() {
    logging::init_logging();
    
    println!("=== Lookup Argument Logging Demo ===");
    println!("Current log level: {:?}", logging::log_level());
    println!();
    
    log_info!("This is an INFO level message - shows important information");
    log_debug!("This is a DEBUG level message - shows detailed debugging info");
    
    println!();
    println!("To control log levels, set the LOOKUP_LOG_LEVEL environment variable:");
    println!("  LOOKUP_LOG_LEVEL=silent");
    println!("  LOOKUP_LOG_LEVEL=info");
    println!("  LOOKUP_LOG_LEVEL=debug");
    println!();
    println!("Example usage:");
    println!("  LOOKUP_LOG_LEVEL=debug cargo run --example logging_demo");
    println!("  LOOKUP_LOG_LEVEL=silent cargo test");
} 