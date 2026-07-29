/// Internal representation of a phylogenetic tree node stored in a flat `Vec<Node>`.
///
/// Indices into the enclosing `Vec` are used instead of pointers.
pub struct Node {
    /// Leaf label, or `None` for internal nodes.
    pub name: Option<String>,
    /// Branch length from this node to its parent (0.0 if absent or for the root).
    pub length: f64,
    /// Index of the parent node, or `None` for the root.
    pub parent: Option<usize>,
    /// Indices of child nodes in DFS left-to-right order.
    pub children: Vec<usize>,
}
