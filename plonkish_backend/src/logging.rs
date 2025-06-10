use std::sync::OnceLock;
use std::env;

/// Log levels for controlling output verbosity
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum LogLevel {
    /// Show no debug messages
    Silent = 0,
    /// Show only important information (default)
    Info = 1,
    /// Show detailed debug information
    Debug = 2,
}

/// Global log level instance
static LOG_LEVEL: OnceLock<LogLevel> = OnceLock::new();

/// Initialize the logging system
/// Reads from environment variable LOOKUP_LOG_LEVEL:
/// - "silent" or "0" -> Silent
/// - "info" or "1" -> Info (default)
/// - "debug" or "2" -> Debug
pub fn init_logging() {
    LOG_LEVEL.get_or_init(|| {
        match env::var("LOOKUP_LOG_LEVEL").as_deref() {
            Ok("silent") | Ok("0") => LogLevel::Silent,
            Ok("debug") | Ok("2") => LogLevel::Debug,
            _ => LogLevel::Info, // Default level
        }
    });
}

/// Get the current log level
pub fn log_level() -> LogLevel {
    *LOG_LEVEL.get_or_init(|| LogLevel::Info)
}

/// Check if info level logging is enabled
pub fn is_info_enabled() -> bool {
    log_level() >= LogLevel::Info
}

/// Check if debug level logging is enabled
pub fn is_debug_enabled() -> bool {
    log_level() >= LogLevel::Debug
}

/// Macro for info level logging
#[macro_export]
macro_rules! log_info {
    ($($arg:tt)*) => {
        if $crate::logging::is_info_enabled() {
            println!("[INFO] {}", format!($($arg)*));
        }
    };
}

/// Macro for debug level logging
#[macro_export]
macro_rules! log_debug {
    ($($arg:tt)*) => {
        if $crate::logging::is_debug_enabled() {
            println!("[DEBUG] {}", format!($($arg)*));
        }
    };
}

/// Conditional debug println - only prints if debug level is enabled
#[macro_export]
macro_rules! debug_println {
    ($($arg:tt)*) => {
        if $crate::logging::is_debug_enabled() {
            println!($($arg)*);
        }
    };
}

/// Initialize logging with explicit level (for testing)
pub fn set_log_level(level: LogLevel) {
    // This won't work if already initialized, but useful for tests
    let _ = LOG_LEVEL.set(level);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_log_levels() {
        // Test the ordering
        assert!(LogLevel::Silent < LogLevel::Info);
        assert!(LogLevel::Info < LogLevel::Debug);
        
        // Test environment variable parsing
        env::set_var("LOOKUP_LOG_LEVEL", "debug");
        // Note: In real usage, you'd need to reset the OnceLock for this to work
        
        env::set_var("LOOKUP_LOG_LEVEL", "info");
        env::set_var("LOOKUP_LOG_LEVEL", "silent");
    }
} 