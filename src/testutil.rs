//! Helpers shared by the test modules: a deterministic PRNG and a random tree
//! generator, so the randomised checks stay reproducible without pulling in a
//! dependency.

/// xorshift64.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(seed)
    }

    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    /// A number in `0..n`.
    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    /// A branch length, zero one time in seven: collapsing weakly supported
    /// branches leaves plenty of zero-length ones in real trees.
    pub fn length(&mut self) -> f64 {
        if self.below(7) == 0 {
            return 0.0;
        }
        0.01 + (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64 * 2.0
    }
}

/// Join random groups of 2 or 3 until one tree remains, so the sample covers
/// polytomies as well as strictly bifurcating trees.
pub fn random_newick(rng: &mut Rng, n_leaves: usize) -> String {
    let mut pool: Vec<String> = (0..n_leaves)
        .map(|i| format!("L{}:{:.6}", i, rng.length()))
        .collect();
    while pool.len() > 1 {
        let arity = if pool.len() > 2 && rng.below(4) == 0 { 3 } else { 2 };
        let children: Vec<String> = (0..arity)
            .map(|_| pool.swap_remove(rng.below(pool.len())))
            .collect();
        pool.push(format!("({}):{:.6}", children.join(","), rng.length()));
    }
    format!("{};", pool.pop().unwrap())
}
