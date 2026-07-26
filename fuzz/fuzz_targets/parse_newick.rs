#![no_main]

//! Throw arbitrary text at the Newick parser.
//!
//! The parser is hand-written, scans bytes, and is deliberately strict, which
//! is exactly the combination that hides an index that runs off the end or a
//! loop that fails to advance. Any input at all must come back as `Ok` or
//! `Err`, never a panic and never a hang.

use libfuzzer_sys::fuzz_target;

use distree::parser::parse_newick;

fuzz_target!(|data: &[u8]| {
    // The parser takes &str; anything that is not UTF-8 is rejected before it
    // ever gets there, so there is nothing to learn from feeding it here.
    let Ok(text) = std::str::from_utf8(data) else {
        return;
    };

    // Very long inputs only slow the fuzzer down; the interesting behaviour is
    // in the structure, not the size.
    if text.len() > 64 * 1024 {
        return;
    }

    let _ = parse_newick(text);
});
