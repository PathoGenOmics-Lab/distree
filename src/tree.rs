/// Internal representation of a tree node.
pub struct Node {
    pub name: Option<String>,
    pub length: f64,
    pub parent: Option<usize>,
    pub children: Vec<usize>,
}
