#[cfg(test)]
mod tests {
    use ntt::{omega, polymul, polymul_ntt};

    #[test]
    fn test_polymul_ntt() {
        let p: i64 = 17; // Prime modulus
        let root: i64 = 3; // Primitive root of unity
        let n: usize = 8;  // Length of the NTT (must be a power of 2)
        let omega = omega(root, p, n); // n-th root of unity

        // Input polynomials (padded to length `n`)
        let mut a = vec![1, 2, 3, 4];
        let mut b = vec![5, 6, 7, 8];
        a.resize(n, 0);
        b.resize(n, 0);

        // Perform the standard polynomial multiplication
        let c_std = polymul(&a, &b, n as i64, p);
        
        // Perform the NTT-based polynomial multiplication
        let c_fast = polymul_ntt(&a, &b, n, p, omega);

        // Ensure both methods produce the same result
        assert_eq!(c_std, c_fast, "The results of polymul and polymul_ntt do not match");
    }

    #[test]
    fn test_polymul_ntt_square_modulus() {
        let modulus: i64 = 17*17; // Prime modulus
        let root: i64 = 3; // Primitive root of unity
        let n: usize = 8;  // Length of the NTT (must be a power of 2)
        let omega = omega(root, modulus, n); // n-th root of unity

        // Input polynomials (padded to length `n`)
        let mut a = vec![1, 2, 3, 4];
        let mut b = vec![5, 6, 7, 8];
        a.resize(n, 0);
        b.resize(n, 0);

        // Perform the standard polynomial multiplication
        let c_std = polymul(&a, &b, n as i64, modulus);
        
        // Perform the NTT-based polynomial multiplication
        let c_fast = polymul_ntt(&a, &b, n, modulus, omega);

        // Ensure both methods produce the same result
        assert_eq!(c_std, c_fast, "The results of polymul and polymul_ntt do not match");
    }

    #[test]
    fn test_polymul_ntt_prime_power_modulus() {
        let modulus: i64 = (17 as i64).pow(4); // modulus p^k
        let root: i64 = 3; // Primitive root of unity
        let n: usize = 8;  // Length of the NTT (must be a power of 2)
        let omega = omega(root, modulus, n); // n-th root of unity

        // Input polynomials (padded to length `n`)
        let mut a = vec![1, 2, 3, 4];
        let mut b = vec![4, 5, 6, 7];
        a.resize(n, 0);
        b.resize(n, 0);

        // Perform the standard polynomial multiplication
        let c_std = polymul(&a, &b, n as i64, modulus);
        
        // Perform the NTT-based polynomial multiplication
        let c_fast = polymul_ntt(&a, &b, n, modulus, omega);

        // Ensure both methods produce the same result
        assert_eq!(c_std, c_fast, "The results of polymul and polymul_ntt do not match");
    }

    #[test]
    fn test_polymul_ntt_twice_prime_power_modulus() {
        let modulus: i64 = 2*(17 as i64).pow(4); // modulus 2*p^k
        let root: i64 = 3; // Primitive root of unity
        let n: usize = 8;  // Length of the NTT (must be a power of 2)
        let omega = omega(root, modulus, n); // n-th root of unity

        // Input polynomials (padded to length `n`)
        let mut a = vec![1, 2, 3, 4];
        let mut b = vec![4, 5, 6, 7];
        a.resize(n, 0);
        b.resize(n, 0);

        // Perform the standard polynomial multiplication
        let c_std = polymul(&a, &b, n as i64, modulus);
        
        // Perform the NTT-based polynomial multiplication
        let c_fast = polymul_ntt(&a, &b, n, modulus, omega);

        // Ensure both methods produce the same result
        assert_eq!(c_std, c_fast, "The results of polymul and polymul_ntt do not match");
    }
}
