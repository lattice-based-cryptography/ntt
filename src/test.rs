#[cfg(test)]
mod tests {
    use ntt::{omega, polymul, polymul_ntt};

    #[test]
    fn test_polymul_ntt() {
        let p: i64 = 17; // Prime modulus
        let n: usize = 8;  // Length of the NTT (must be a power of 2)
        let omega = omega(p, n); // n-th root of unity

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
        let cases = [
            (17*17, 4),        // small square modulus
            (12289*12289, 512) // large square modulus
        ];
        
        for &(modulus, n) in &cases {
            let omega = omega(modulus, 2*n); // n-th root of unity
            let mut a: Vec<i64> = (0..n).map(|x| x as i64).collect();
            let mut b: Vec<i64> = (0..n).map(|x| x as i64).collect();
            a.resize(2*n, 0);
            b.resize(2*n, 0);
        
            let c_std = polymul(&a, &b, 2*n as i64, modulus);
            let c_fast = polymul_ntt(&a, &b, 2*n, modulus, omega);
        
            assert_eq!(c_std, c_fast, "The results of polymul and polymul_ntt do not match for modulus {} and n {}", modulus, n);
        }

    }

    #[test]
    fn test_polymul_ntt_prime_power_modulus() {
        let modulus: i64 = (17 as i64).pow(4); // modulus p^k
        let n: usize = 8;  // Length of the NTT (must be a power of 2)
        let omega = omega(modulus, n); // n-th root of unity

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
    fn test_polymul_ntt_non_prime_power_modulus() {
        let moduli = [17*41, 17*73, 17*41*73]; // Different moduli to test
        let n: usize = 8;  // Length of the NTT (must be a power of 2)
    
        for &modulus in &moduli {
            let omega = omega(modulus, n);
            
            let mut a = vec![1, 2, 3, 4];
            let mut b = vec![4, 5, 6, 7];
            a.resize(n, 0);
            b.resize(n, 0);
    
            let c_std = polymul(&a, &b, n as i64, modulus);
            let c_fast = polymul_ntt(&a, &b, n, modulus, omega);
    
            assert_eq!(c_std, c_fast, "Failed for modulus {}", modulus);
        }
    }
    
}
