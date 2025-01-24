#[cfg(test)]
mod tests {
    use ntt::{polymul, polymul_ntt};

    #[test]
    fn test_polymul_ntt() {
        let p: u64 = 17; // Prime modulus
        let root: u64 = 3; // Primitive root of unity
        let n: usize = 8; // Length of the NTT (must be a power of 2)

        // Input polynomials (padded to length `n`)
        let mut a = vec![1, 2, 3, 4];
        let mut b = vec![5, 6, 7, 8];
        a.resize(n, 0);
        b.resize(n, 0);

        // Perform the standard polynomial multiplication
        let c_std = polymul(&a, &b, n as u64, p);
        
        // Perform the NTT-based polynomial multiplication
        let c_fast = polymul_ntt(&a, &b, n, p, root);

        // Ensure both methods produce the same result
        assert_eq!(c_std, c_fast, "The results of polymul and polymul_ntt do not match");
    }
}