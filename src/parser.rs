use std::iter::Peekable;
use std::num::ParseFloatError;
use std::str::Chars;

use crate::tree::Node;

/// Temporary structure used during Newick parsing.
pub struct RawNode {
    pub name: Option<String>,
    pub length: f64,
    pub children: Vec<RawNode>,
}

/// Recursively parse a Newick subtree and return a `RawNode`.
pub fn parse_subtree(chars: &mut Peekable<Chars>) -> Result<RawNode, String> {
    let mut node = RawNode {
        name: None,
        length: 0.0,
        children: Vec::new(),
    };

    if let Some(&c) = chars.peek() {
        if c == '(' {
            chars.next();
            loop {
                let child = parse_subtree(chars)?;
                node.children.push(child);
                match chars.peek() {
                    Some(',') => {
                        chars.next();
                        continue;
                    }
                    Some(')') => {
                        chars.next();
                        break;
                    }
                    Some(other) => {
                        return Err(format!("Expected ',' or ')', found '{}'", other));
                    }
                    None => return Err("Unexpected end of input after '('".to_string()),
                }
            }

            // After ')', skip any NHX/bracket comments
            skip_comments(chars);

            // Optional internal node label
            if let Some(&c2) = chars.peek() {
                if c2 != ':' && c2 != ',' && c2 != ')' && c2 != ';' {
                    let name = parse_label(chars);
                    if !name.is_empty() {
                        node.name = Some(name);
                    }
                }
            }

            // Skip comments after label too
            skip_comments(chars);

            if let Some(&':') = chars.peek() {
                chars.next();
                let length = parse_length(chars)?;
                node.length = length;
            }

            // Skip comments after branch length
            skip_comments(chars);

            return Ok(node);
        }
    }

    // Leaf node
    let name = parse_label(chars);
    if !name.is_empty() {
        node.name = Some(name);
    }

    // Skip comments after label
    skip_comments(chars);

    if let Some(&':') = chars.peek() {
        chars.next();
        let length = parse_length(chars)?;
        node.length = length;
    }

    // Skip comments after branch length
    skip_comments(chars);

    Ok(node)
}

/// Parse a node label. Supports single-quoted labels (strips quotes).
pub fn parse_label(chars: &mut Peekable<Chars>) -> String {
    // Handle single-quoted labels
    if let Some(&c) = chars.peek() {
        if c == '\'' {
            chars.next(); // consume opening quote
            let mut label = String::new();
            while let Some(&ch) = chars.peek() {
                if ch == '\'' {
                    chars.next(); // consume closing quote
                    break;
                }
                label.push(ch);
                chars.next();
            }
            return label;
        }
    }

    let mut label = String::new();
    while let Some(&c) = chars.peek() {
        if c == ':' || c == ',' || c == ')' || c == ';' || c == '[' || c.is_whitespace() {
            break;
        }
        label.push(c);
        chars.next();
    }
    label
}

/// Skip '[...]' comment blocks (NHX annotations, BEAST metadata, etc.).
fn skip_comments(chars: &mut Peekable<Chars>) {
    while let Some(&'[') = chars.peek() {
        chars.next();
        let mut depth = 1;
        while depth > 0 {
            match chars.next() {
                Some('[') => depth += 1,
                Some(']') => depth -= 1,
                None => break,
                _ => {}
            }
        }
    }
}

/// Parse a floating-point branch length (supports scientific notation).
pub fn parse_length(chars: &mut Peekable<Chars>) -> Result<f64, String> {
    let mut numstr = String::new();
    while let Some(&c) = chars.peek() {
        if c.is_ascii_digit() || c == '.' || c == 'e' || c == 'E' || c == '+' || c == '-' {
            numstr.push(c);
            chars.next();
        } else {
            break;
        }
    }
    numstr
        .parse::<f64>()
        .map_err(|e: ParseFloatError| format!("Failed to parse branch length '{}': {}", numstr, e))
}

/// Recursively flatten a `RawNode` into a flat `Vec<Node>`, returning the index of the new node.
pub fn flatten_raw(raw: &RawNode, parent: Option<usize>, nodes: &mut Vec<Node>) -> usize {
    let idx = nodes.len();
    nodes.push(Node {
        name: raw.name.clone(),
        length: raw.length,
        parent,
        children: Vec::new(),
    });
    if let Some(p) = parent {
        nodes[p].children.push(idx);
    }
    for child in &raw.children {
        flatten_raw(child, Some(idx), nodes);
    }
    idx
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_simple_tree() {
        let input = "(A:1.0,B:2.0);";
        let mut chars = input.chars().peekable();
        let raw = parse_subtree(&mut chars).unwrap();
        assert_eq!(raw.children.len(), 2);
        assert_eq!(raw.children[0].name.as_deref(), Some("A"));
        assert_eq!(raw.children[0].length, 1.0);
        assert_eq!(raw.children[1].name.as_deref(), Some("B"));
        assert_eq!(raw.children[1].length, 2.0);
    }

    #[test]
    fn test_parse_internal_labels() {
        let input = "((A:1.0,B:2.0)inner:0.5,C:3.0);";
        let mut chars = input.chars().peekable();
        let raw = parse_subtree(&mut chars).unwrap();
        assert_eq!(raw.children.len(), 2);
        assert_eq!(raw.children[0].name.as_deref(), Some("inner"));
        assert_eq!(raw.children[0].length, 0.5);
        assert_eq!(raw.children[0].children.len(), 2);
    }

    #[test]
    fn test_parse_no_branch_lengths() {
        let input = "(A,B,(C,D));";
        let mut chars = input.chars().peekable();
        let raw = parse_subtree(&mut chars).unwrap();
        assert_eq!(raw.children.len(), 3);
        assert_eq!(raw.children[0].length, 0.0);
        assert_eq!(raw.children[2].children.len(), 2);
    }

    #[test]
    fn test_parse_nhx_comments() {
        let input = "((A:0.1[&&NHX:S=human],B:0.2[&&NHX:S=mouse]):0.3,C:0.4);";
        let mut chars = input.chars().peekable();
        let raw = parse_subtree(&mut chars).unwrap();
        assert_eq!(raw.children.len(), 2);
        let inner = &raw.children[0];
        assert_eq!(inner.children[0].name.as_deref(), Some("A"));
        assert_eq!(inner.children[0].length, 0.1);
        assert_eq!(inner.children[1].name.as_deref(), Some("B"));
        assert_eq!(inner.children[1].length, 0.2);
        assert_eq!(inner.length, 0.3);
    }

    #[test]
    fn test_parse_quoted_labels() {
        let input = "('Taxon A':1.0,'Taxon B':2.0);";
        let mut chars = input.chars().peekable();
        let raw = parse_subtree(&mut chars).unwrap();
        assert_eq!(raw.children[0].name.as_deref(), Some("Taxon A"));
        assert_eq!(raw.children[1].name.as_deref(), Some("Taxon B"));
    }

    #[test]
    fn test_parse_scientific_notation() {
        let input = "(A:1.5e-3,B:2.0E+1);";
        let mut chars = input.chars().peekable();
        let raw = parse_subtree(&mut chars).unwrap();
        assert!((raw.children[0].length - 0.0015).abs() < 1e-10);
        assert!((raw.children[1].length - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_flatten() {
        let input = "((A:1.0,B:2.0):0.5,C:3.0);";
        let mut chars = input.chars().peekable();
        let raw = parse_subtree(&mut chars).unwrap();
        let mut nodes = Vec::new();
        let root = flatten_raw(&raw, None, &mut nodes);
        assert_eq!(root, 0);
        // root -> (inner, C), inner -> (A, B) = 5 nodes
        assert_eq!(nodes.len(), 5);
        assert!(nodes[0].parent.is_none());
        assert_eq!(nodes[0].children.len(), 2);
    }

    #[test]
    fn test_bracket_in_label_position() {
        let input = "((A:0.1[comment],B:0.2):0.3[more],C:0.4[x]);";
        let mut chars = input.chars().peekable();
        let raw = parse_subtree(&mut chars).unwrap();
        assert_eq!(raw.children[0].children[0].name.as_deref(), Some("A"));
        assert_eq!(raw.children[1].name.as_deref(), Some("C"));
    }
}
