use std::num::ParseFloatError;

use crate::tree::Node;

/// Temporary recursive structure used during Newick parsing.
/// Converted to a flat `Vec<Node>` by [`flatten_raw`] after parsing.
pub struct RawNode {
    pub name: Option<String>,
    pub length: f64,
    pub children: Vec<RawNode>,
}

impl Drop for RawNode {
    /// Dismantle the subtree iteratively.
    ///
    /// The compiler-generated drop glue recurses once per level, which
    /// overflows the stack on the deeply nested trees the parser and the
    /// flattener were made iterative to support.
    fn drop(&mut self) {
        let mut pending = std::mem::take(&mut self.children);
        while let Some(mut node) = pending.pop() {
            pending.append(&mut node.children);
        }
    }
}

/// Parse a complete Newick tree from a string.
///
/// Strips whitespace outside of quoted labels before parsing.
/// Accepts trees with or without a trailing semicolon.
/// Returns an error on malformed input or empty input.
pub fn parse_newick(input: &str) -> Result<RawNode, String> {
    let trimmed = input.trim();
    if trimmed.is_empty() {
        return Err("Empty input: no Newick tree found.".to_string());
    }
    // Pre-process: strip whitespace outside of single/double quotes
    let cleaned = strip_whitespace_outside_quotes(trimmed);
    let bytes = cleaned.as_bytes();
    let mut pos = 0;
    let root = parse_subtree_iterative(bytes, &mut pos)?;

    // Only an optional ';' and trailing comments may follow the tree. Anything
    // else means the input is not a single well-formed tree, and silently
    // ignoring it would hand back a matrix for something the user never asked
    // for.
    skip_ignorable_bytes(bytes, &mut pos)?;
    let mut ended_with_semicolon = false;
    if pos < bytes.len() && bytes[pos] == b';' {
        pos += 1;
        ended_with_semicolon = true;
    }
    skip_ignorable_bytes(bytes, &mut pos)?;
    if pos < bytes.len() {
        if ended_with_semicolon && bytes[pos] == b'(' {
            return Err(format!(
                "Unexpected content at position {}: the file appears to hold more than one tree. \
                 distree processes a single Newick tree; split the file first.",
                pos
            ));
        }
        return Err(format!(
            "Unexpected content after the end of the tree at position {}: '{}'.",
            pos,
            bytes[pos] as char
        ));
    }

    Ok(root)
}

/// Strip whitespace characters outside of quoted regions.
///
/// Works on raw bytes, which is safe for UTF-8 input: every byte of a
/// multi-byte sequence has the high bit set and therefore never matches an
/// ASCII quote or whitespace byte. Only whole ASCII bytes are dropped, so the
/// remaining bytes are still valid UTF-8.
fn strip_whitespace_outside_quotes(s: &str) -> String {
    let src = s.as_bytes();
    let mut out: Vec<u8> = Vec::with_capacity(src.len());
    // Last non-whitespace byte emitted. A quote only opens a quoted label
    // where a label may begin: at the very start of the tree, or right after
    // '(', ',' or ')'. Anywhere else an apostrophe is an ordinary character of
    // an unquoted label ("O'Brien") or of a comment ("[don't]").
    let mut previous: Option<u8> = None;
    let mut i = 0;

    while i < src.len() {
        let b = src[i];
        let opens_label = matches!(b, b'\'' | b'"')
            && matches!(previous, None | Some(b'(') | Some(b',') | Some(b')'));

        if opens_label {
            // Copy the quoted label verbatim so its whitespace survives,
            // treating a doubled quote as an escaped literal rather than as a
            // close followed by a re-open.
            out.push(b);
            i += 1;
            while i < src.len() {
                if src[i] == b {
                    if src.get(i + 1) == Some(&b) {
                        out.extend_from_slice(&[b, b]);
                        i += 2;
                        continue;
                    }
                    out.push(b);
                    i += 1;
                    break;
                }
                out.push(src[i]);
                i += 1;
            }
            previous = Some(b);
        } else if b.is_ascii_whitespace() {
            i += 1;
        } else {
            out.push(b);
            previous = Some(b);
            i += 1;
        }
    }

    String::from_utf8(out).expect("dropping ASCII whitespace preserves UTF-8 validity")
}

/// Iterative Newick parser — no stack overflow on deep trees.
fn parse_subtree_iterative(bytes: &[u8], pos: &mut usize) -> Result<RawNode, String> {
    // We use an explicit stack to avoid recursion.
    // Each frame represents a node being constructed.
    struct Frame {
        node: RawNode,
    }

    let mut stack: Vec<Frame> = Vec::new();

    loop {
        // A subtree may be preceded by comments: the rooting marker some
        // programs emit ("[&R] (A,B);") or per-branch metadata ("(A,[&x=1]B)").
        skip_ignorable_bytes(bytes, pos)?;

        // Decide what to parse at current position
        let current_node = if *pos < bytes.len() && bytes[*pos] == b'(' {
            // Internal node: push frame, advance past '('
            *pos += 1;
            skip_whitespace_bytes(bytes, pos);
            stack.push(Frame {
                node: RawNode { name: None, length: 0.0, children: Vec::new() },
            });
            continue; // loop back to parse first child
        } else {
            // Leaf node
            let name = parse_label_bytes(bytes, pos)?;
            skip_comments_bytes(bytes, pos)?;
            let length = parse_optional_length(bytes, pos)?;
            skip_comments_bytes(bytes, pos)?;
            let mut leaf = RawNode { name: None, length, children: Vec::new() };
            if !name.is_empty() {
                leaf.name = Some(name);
            }
            leaf
        };

        // Now we have a completed node. Feed it back up the stack.
        let mut completed = current_node;
        loop {
            if let Some(frame) = stack.last_mut() {
                frame.node.children.push(completed);
                // Check what comes next: ',' means more children, ')' means done
                if *pos < bytes.len() && bytes[*pos] == b',' {
                    *pos += 1;
                    skip_whitespace_bytes(bytes, pos);
                    break; // break inner loop, continue outer to parse next child
                } else if *pos < bytes.len() && bytes[*pos] == b')' {
                    *pos += 1;
                    // Internal node complete — read its label and length
                    skip_comments_bytes(bytes, pos)?;
                    let name = parse_label_bytes(bytes, pos)?;
                    skip_comments_bytes(bytes, pos)?;
                    let length = parse_optional_length(bytes, pos)?;
                    skip_comments_bytes(bytes, pos)?;

                    let mut frame = stack.pop().unwrap();
                    if !name.is_empty() {
                        frame.node.name = Some(name);
                    }
                    frame.node.length = length;
                    completed = frame.node;
                    // continue inner loop to feed this up further
                } else if *pos >= bytes.len() {
                    // End of input while still inside parentheses: the tree is
                    // truncated. Closing it silently would emit a matrix whose
                    // branch lengths are quietly wrong.
                    return Err(format!(
                        "Unexpected end of input: {} unclosed '(' remain. The tree is truncated.",
                        stack.len()
                    ));
                } else {
                    let ch = bytes[*pos] as char;
                    return Err(format!("Expected ',' or ')', found '{}' at position {}", ch, *pos));
                }
            } else {
                // Stack empty — this is the root
                return Ok(completed);
            }
        }
    }
}

fn skip_whitespace_bytes(bytes: &[u8], pos: &mut usize) {
    while *pos < bytes.len() && bytes[*pos].is_ascii_whitespace() {
        *pos += 1;
    }
}

/// Skip any run of whitespace and '[...]' comments.
fn skip_ignorable_bytes(bytes: &[u8], pos: &mut usize) -> Result<(), String> {
    loop {
        let before = *pos;
        skip_whitespace_bytes(bytes, pos);
        skip_comments_bytes(bytes, pos)?;
        if *pos == before {
            return Ok(());
        }
    }
}

/// Parse a node label from bytes. Supports single-quoted and double-quoted labels.
/// Returns an error if a quoted label is not properly closed.
fn parse_label_bytes(bytes: &[u8], pos: &mut usize) -> Result<String, String> {
    if *pos >= bytes.len() {
        return Ok(String::new());
    }

    // Handle quoted labels (single or double quotes)
    // Per the Newick spec, quotes inside a quoted label are escaped by doubling:
    // 'it''s a name' → it's a name
    let ch = bytes[*pos];
    if ch == b'\'' || ch == b'"' {
        let open_pos = *pos;
        *pos += 1; // consume opening quote
        let mut label: Vec<u8> = Vec::new();
        let mut closed = false;
        while *pos < bytes.len() {
            if bytes[*pos] == ch {
                // Check for escaped quote (doubled)
                if *pos + 1 < bytes.len() && bytes[*pos + 1] == ch {
                    label.push(ch);
                    *pos += 2; // skip both quotes
                } else {
                    *pos += 1; // consume closing quote
                    closed = true;
                    break;
                }
            } else {
                label.push(bytes[*pos]);
                *pos += 1;
            }
        }
        if !closed {
            return Err(format!(
                "Unclosed quote starting at position {}.",
                open_pos
            ));
        }
        return label_from_utf8(label, open_pos);
    }

    // Unquoted label: read until delimiter
    let start = *pos;
    while *pos < bytes.len() {
        let c = bytes[*pos];
        if c == b':'
            || c == b','
            || c == b')'
            || c == b'('
            || c == b';'
            || c == b'['
            || c.is_ascii_whitespace()
        {
            break;
        }
        *pos += 1;
    }
    label_from_utf8(bytes[start..*pos].to_vec(), start)
}

/// Turn the raw bytes of a label into a `String`.
///
/// Labels are scanned byte-wise because every Newick delimiter is ASCII, but
/// the bytes in between may encode any UTF-8 text and must be decoded as such
/// rather than reinterpreted one byte at a time.
fn label_from_utf8(bytes: Vec<u8>, start: usize) -> Result<String, String> {
    String::from_utf8(bytes).map_err(|_| {
        format!(
            "Label starting at position {} is not valid UTF-8.",
            start
        )
    })
}

/// Skip '[...]' comment blocks (NHX annotations, BEAST metadata, etc.).
/// Handles nested brackets. Errors if a comment is never closed, since the
/// rest of the tree would otherwise be swallowed as comment text.
fn skip_comments_bytes(bytes: &[u8], pos: &mut usize) -> Result<(), String> {
    while *pos < bytes.len() && bytes[*pos] == b'[' {
        let open_pos = *pos;
        *pos += 1;
        let mut depth: usize = 1;
        while *pos < bytes.len() && depth > 0 {
            match bytes[*pos] {
                b'[' => depth += 1,
                b']' => depth -= 1,
                _ => {}
            }
            *pos += 1;
        }
        if depth > 0 {
            return Err(format!(
                "Unclosed comment starting at position {}: no matching ']'.",
                open_pos
            ));
        }
    }
    Ok(())
}

/// Parse optional ":length" — returns 0.0 if no colon present.
fn parse_optional_length(bytes: &[u8], pos: &mut usize) -> Result<f64, String> {
    if *pos < bytes.len() && bytes[*pos] == b':' {
        *pos += 1;
        parse_length_bytes(bytes, pos)
    } else {
        Ok(0.0)
    }
}

/// Parse a floating-point branch length (supports scientific notation and negative values).
/// Returns an error if no valid numeric characters follow the colon.
fn parse_length_bytes(bytes: &[u8], pos: &mut usize) -> Result<f64, String> {
    let start = *pos;
    while *pos < bytes.len() {
        let c = bytes[*pos];
        if c.is_ascii_digit() || c == b'.' || c == b'e' || c == b'E' || c == b'+' || c == b'-' {
            *pos += 1;
        } else {
            break;
        }
    }
    if start == *pos {
        return Err(format!(
            "Expected a numeric branch length at position {}, found '{}'.",
            pos,
            bytes.get(*pos).map(|&b| b as char).unwrap_or('\0')
        ));
    }
    let numstr = std::str::from_utf8(&bytes[start..*pos])
        .map_err(|_| "Invalid UTF-8 in branch length".to_string())?;
    let value = numstr
        .parse::<f64>()
        .map_err(|e: ParseFloatError| {
            format!("Failed to parse branch length '{}': {}", numstr, e)
        })?;
    // Rust parses an exponent past f64's range as infinity rather than as an
    // error, and an infinite branch length turns the whole matrix into inf with
    // NaN down the diagonal, since inf - inf is NaN.
    if !value.is_finite() {
        return Err(format!(
            "Branch length '{}' at position {} is not a finite number.",
            numstr, start
        ));
    }
    Ok(value)
}

/// Flatten a `RawNode` tree into a flat `Vec<Node>` iteratively (no stack overflow).
///
/// Returns the index of the root node in `nodes`.
/// Parent/child relationships are set up correctly for LCA queries.
pub fn flatten_raw(raw: &RawNode, parent: Option<usize>, nodes: &mut Vec<Node>) -> usize {
    // Stack of (raw_node_ref, parent_index)
    let mut stack: Vec<(&RawNode, Option<usize>)> = vec![(raw, parent)];
    let root_idx = nodes.len();

    while let Some((raw_node, par)) = stack.pop() {
        let idx = nodes.len();
        nodes.push(Node {
            name: raw_node.name.clone(),
            length: raw_node.length,
            parent: par,
            children: Vec::new(),
        });
        if let Some(p) = par {
            nodes[p].children.push(idx);
        }
        // Push children in reverse order so they're processed left-to-right
        for child in raw_node.children.iter().rev() {
            stack.push((child, Some(idx)));
        }
    }

    root_idx
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_simple_tree() {
        let raw = parse_newick("(A:1.0,B:2.0);").unwrap();
        assert_eq!(raw.children.len(), 2);
        assert_eq!(raw.children[0].name.as_deref(), Some("A"));
        assert_eq!(raw.children[0].length, 1.0);
        assert_eq!(raw.children[1].name.as_deref(), Some("B"));
        assert_eq!(raw.children[1].length, 2.0);
    }

    #[test]
    fn test_parse_internal_labels() {
        let raw = parse_newick("((A:1.0,B:2.0)inner:0.5,C:3.0);").unwrap();
        assert_eq!(raw.children.len(), 2);
        assert_eq!(raw.children[0].name.as_deref(), Some("inner"));
        assert_eq!(raw.children[0].length, 0.5);
        assert_eq!(raw.children[0].children.len(), 2);
    }

    #[test]
    fn test_parse_no_branch_lengths() {
        let raw = parse_newick("(A,B,(C,D));").unwrap();
        assert_eq!(raw.children.len(), 3);
        assert_eq!(raw.children[0].length, 0.0);
        assert_eq!(raw.children[2].children.len(), 2);
    }

    #[test]
    fn test_parse_nhx_comments() {
        let raw = parse_newick("((A:0.1[&&NHX:S=human],B:0.2[&&NHX:S=mouse]):0.3,C:0.4);").unwrap();
        assert_eq!(raw.children.len(), 2);
        let inner = &raw.children[0];
        assert_eq!(inner.children[0].name.as_deref(), Some("A"));
        assert_eq!(inner.children[0].length, 0.1);
        assert_eq!(inner.children[1].name.as_deref(), Some("B"));
        assert_eq!(inner.children[1].length, 0.2);
        assert_eq!(inner.length, 0.3);
    }

    #[test]
    fn test_parse_quoted_labels_single() {
        let raw = parse_newick("('Taxon A':1.0,'Taxon B':2.0);").unwrap();
        assert_eq!(raw.children[0].name.as_deref(), Some("Taxon A"));
        assert_eq!(raw.children[1].name.as_deref(), Some("Taxon B"));
    }

    #[test]
    fn test_parse_quoted_labels_double() {
        let raw = parse_newick(r#"("Taxon A":1.0,"Taxon B":2.0);"#).unwrap();
        assert_eq!(raw.children[0].name.as_deref(), Some("Taxon A"));
        assert_eq!(raw.children[1].name.as_deref(), Some("Taxon B"));
    }

    #[test]
    fn test_non_ascii_labels_preserved() {
        let raw = parse_newick("((Señor_Ñu:1.0,'β strain':2.0):0.5,日本株:3.0);").unwrap();
        let inner = &raw.children[0];
        assert_eq!(inner.children[0].name.as_deref(), Some("Señor_Ñu"));
        assert_eq!(inner.children[1].name.as_deref(), Some("β strain"));
        assert_eq!(raw.children[1].name.as_deref(), Some("日本株"));
    }

    #[test]
    fn test_parse_scientific_notation() {
        let raw = parse_newick("(A:1.5e-3,B:2.0E+1);").unwrap();
        assert!((raw.children[0].length - 0.0015).abs() < 1e-10);
        assert!((raw.children[1].length - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_flatten() {
        let raw = parse_newick("((A:1.0,B:2.0):0.5,C:3.0);").unwrap();
        let mut nodes = Vec::new();
        let root = flatten_raw(&raw, None, &mut nodes);
        assert_eq!(root, 0);
        assert_eq!(nodes.len(), 5);
        assert!(nodes[0].parent.is_none());
        assert_eq!(nodes[0].children.len(), 2);
    }

    #[test]
    fn test_bracket_in_label_position() {
        let raw = parse_newick("((A:0.1[comment],B:0.2):0.3[more],C:0.4[x]);").unwrap();
        assert_eq!(raw.children[0].children[0].name.as_deref(), Some("A"));
        assert_eq!(raw.children[1].name.as_deref(), Some("C"));
    }

    #[test]
    fn test_whitespace_in_newick() {
        let raw = parse_newick("( A:1.0 , B:2.0 ) ;").unwrap();
        assert_eq!(raw.children.len(), 2);
        assert_eq!(raw.children[0].name.as_deref(), Some("A"));
        assert_eq!(raw.children[1].name.as_deref(), Some("B"));
    }

    #[test]
    fn test_newlines_in_newick() {
        let raw = parse_newick("((A:1.0,\nB:2.0):0.5,\nC:3.0);\n").unwrap();
        assert_eq!(raw.children.len(), 2);
    }

    #[test]
    fn test_deep_tree_no_stack_overflow() {
        // Build a caterpillar tree with 5,000 levels using push (no recursive format!)
        let mut tree = String::with_capacity(200_000);
        for _ in 0..5_000 {
            tree.push('(');
        }
        tree.push_str("A:0.1");
        for i in 1..=5_000 {
            tree.push_str(&format!(",T{}:0.1):0.01", i));
        }
        tree.push(';');
        let raw = parse_newick(&tree).unwrap();
        // Should not stack overflow; just verify it parsed
        let mut nodes = Vec::new();
        let root = flatten_raw(&raw, None, &mut nodes);
        assert!(nodes.len() > 5_000);
        assert!(nodes[root].parent.is_none());
    }

    #[test]
    fn test_leading_rooting_comment() {
        // IQ-TREE, MrBayes and BEAST prefix the tree with a rooting marker
        let raw = parse_newick("[&R] ((A:1.0,B:2.0):0.5,C:3.0);").unwrap();
        assert_eq!(raw.children.len(), 2);
        assert_eq!(raw.children[1].name.as_deref(), Some("C"));
        assert_eq!(raw.children[0].children[0].name.as_deref(), Some("A"));
    }

    #[test]
    fn test_comment_before_label() {
        let raw = parse_newick("(A:1.0,[&x=1]B:2.0);").unwrap();
        assert_eq!(raw.children[1].name.as_deref(), Some("B"));
        assert_eq!(raw.children[1].length, 2.0);
    }

    #[test]
    fn test_comment_before_subtree() {
        let raw = parse_newick("([&clade=1](A:1.0,B:2.0):0.5,C:3.0);").unwrap();
        assert_eq!(raw.children.len(), 2);
        assert_eq!(raw.children[0].children.len(), 2);
        assert_eq!(raw.children[0].length, 0.5);
    }

    #[test]
    fn test_deeply_nested_tree_drops_without_overflow() {
        // A ladder 200,000 levels deep: parsing, flattening and *dropping*
        // the parse tree all have to stay iterative
        let depth = 200_000;
        let mut tree = String::with_capacity(2 * depth + 8);
        for _ in 0..depth {
            tree.push('(');
        }
        tree.push_str("A:0.1");
        for _ in 0..depth {
            tree.push(')');
        }
        tree.push(';');

        let raw = parse_newick(&tree).unwrap();
        let mut nodes = Vec::new();
        let root = flatten_raw(&raw, None, &mut nodes);
        assert_eq!(nodes.len(), depth + 1);
        assert!(nodes[root].parent.is_none());
        drop(raw);
    }

    #[test]
    fn test_nested_brackets() {
        let raw = parse_newick("((A:0.1[&rate=0.5[inner]],B:0.2):0.3,C:0.4);").unwrap();
        assert_eq!(raw.children[0].children[0].name.as_deref(), Some("A"));
    }

    #[test]
    fn test_negative_branch_length_parsed() {
        let raw = parse_newick("(A:-0.5,B:2.0);").unwrap();
        assert!((raw.children[0].length - (-0.5)).abs() < 1e-10);
    }

    #[test]
    fn test_empty_input_error() {
        assert!(parse_newick("").is_err());
        assert!(parse_newick("   ").is_err());
        assert!(parse_newick("\n\t").is_err());
    }

    #[test]
    fn test_empty_branch_length_error() {
        // ":" followed by non-numeric should error, not panic
        let result = parse_newick("(A:,B:2.0);");
        assert!(result.is_err(), "Expected error for empty branch length");
    }

    #[test]
    fn test_escaped_single_quotes_in_labels() {
        // Newick spec: doubled quotes inside a quoted label are a single literal quote
        // 'it''s' -> it's
        let raw = parse_newick(r"('it''s':1.0,B:2.0);").unwrap();
        assert_eq!(raw.children[0].name.as_deref(), Some("it's"));
        assert_eq!(raw.children[1].name.as_deref(), Some("B"));
    }

    #[test]
    fn test_escaped_double_quotes_in_labels() {
        // "a ""name""" -> a "name"
        let input = r#"("a ""name""":1.0,B:2.0);"#;
        let raw = parse_newick(input).unwrap();
        assert_eq!(raw.children[0].name.as_deref(), Some(r#"a "name""#));
    }

    #[test]
    fn test_apostrophe_in_unquoted_label() {
        // An apostrophe mid-label must not be read as an opening quote
        let raw = parse_newick("(O'Brien:1.0 , B:2.0);").unwrap();
        assert_eq!(raw.children[0].name.as_deref(), Some("O'Brien"));
        assert_eq!(raw.children[1].name.as_deref(), Some("B"));
    }

    #[test]
    fn test_apostrophe_in_comment() {
        let raw = parse_newick("(A:1.0[don't panic],B:2.0);").unwrap();
        assert_eq!(raw.children[0].name.as_deref(), Some("A"));
        assert_eq!(raw.children[1].length, 2.0);
    }

    #[test]
    fn test_escaped_quote_keeps_inner_whitespace() {
        // The doubled quote must not end the label, or the space after it
        // would be stripped as if it were outside the quotes
        let raw = parse_newick(r"('it''s a name':1.0,B:2.0);").unwrap();
        assert_eq!(raw.children[0].name.as_deref(), Some("it's a name"));
    }

    #[test]
    fn test_quoted_internal_node_label() {
        let raw = parse_newick("((A:1.0,B:2.0)'my clade':0.5,C:3.0);").unwrap();
        assert_eq!(raw.children[0].name.as_deref(), Some("my clade"));
    }

    #[test]
    fn test_no_trailing_semicolon() {
        // Trees without semicolons are valid in many tools
        let result = parse_newick("(A:1.0,B:2.0)");
        // Should parse successfully
        assert!(result.is_ok());
        let raw = result.unwrap();
        assert_eq!(raw.children.len(), 2);
    }

    #[test]
    fn test_unclosed_quote_error() {
        let result = parse_newick("('unclosed:1.0,B:2.0);");
        assert!(result.is_err(), "Unclosed quote should error");
        let msg = result.err().expect("should be Err");
        assert!(msg.contains("Unclosed quote"), "Error message: {}", msg);
    }

    #[test]
    fn test_unmatched_paren_error() {
        // A truncated tree must be rejected: closing it silently would give a
        // matrix whose branch lengths are wrong without any warning.
        let err = parse_newick("((A:1.0,B:2.0)").err().expect("truncated tree should error");
        assert!(err.contains("truncated"), "Error message: {}", err);
        assert!(parse_newick("((A:1.0,B:2.0):0.5,(C:1.0,D:2.0").is_err());
    }

    #[test]
    fn test_trailing_content_error() {
        let err = parse_newick("(A:1.0,B:2.0);garbage").err().expect("trailing junk should error");
        assert!(err.contains("Unexpected content"), "Error message: {}", err);
    }

    #[test]
    fn test_multiple_trees_error() {
        let err = parse_newick("(A:1,B:2);\n(C:1,D:2);\n")
            .err().expect("a multi-tree file should error");
        assert!(err.contains("more than one tree"), "Error message: {}", err);
    }

    #[test]
    fn test_free_text_is_not_a_leaf() {
        // Plain prose used to parse as one enormous leaf label
        assert!(parse_newick("not a tree at all (((").is_err());
    }

    #[test]
    fn test_unclosed_comment_error() {
        let err = parse_newick("((A:1.0,B:2.0):0.5,C:3.0)[unclosed;")
            .err().expect("unclosed comment should error");
        assert!(err.contains("Unclosed comment"), "Error message: {}", err);
    }

    #[test]
    fn test_trailing_comment_allowed() {
        let raw = parse_newick("(A:1.0,B:2.0); [ generated by iqtree ]").unwrap();
        assert_eq!(raw.children.len(), 2);
    }
}
