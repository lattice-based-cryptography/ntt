// Modular arithmetic functions
fn mod_add(a: u64, b: u64, p: u64) -> u64 {
    (a + b) % p
}

fn mod_mul(a: u64, b: u64, p: u64) -> u64 {
    (a * b) % p
}

pub fn mod_exp(mut base: u64, mut exp: u64, p: u64) -> u64 {
    let mut result = 1;
    base %= p;
    while exp > 0 {
        if exp % 2 == 1 {
            result = mod_mul(result, base, p);
        }
        base = mod_mul(base, base, p);
        exp /= 2;
    }
    result
}

fn mod_inv(a: u64, p: u64) -> u64 {
    mod_exp(a, p - 2, p) // Using Fermat's Little Theorem
}

// Forward transform using NTT
pub fn ntt(a: &[u64], omega: u64, n: usize, p: u64) -> Vec<u64> {
    let mut result = vec![0; n];
    for k in 0..n {
        let mut value = 0;
        for j in 0..n {
            value = mod_add(value, mod_mul(a[j], mod_exp(omega, (j * k) as u64, p), p), p);
        }
        result[k] = value;
    }
    result
}

// Inverse transform using INTT
pub fn intt(a: &[u64], omega: u64, n: usize, p: u64) -> Vec<u64> {
    let omega_inv = mod_inv(omega, p);
    let n_inv = mod_inv(n as u64, p);
    let mut result = vec![0; n];
    for k in 0..n {
        let mut value = 0;
        for j in 0..n {
            value = mod_add(value, mod_mul(a[j], mod_exp(omega_inv, (j * k) as u64, p), p), p);
        }
        result[k] = mod_mul(value, n_inv, p);
    }
    result
}

// Naive polynomial multiplication
pub fn polymul(a: &Vec<u64>, b: &Vec<u64>, n: u64, p: u64) -> Vec<u64> {
    let mut result = vec![0; n as usize];
    for i in 0..a.len() {
        for j in 0..b.len() {
            result[(i + j) % n as usize] = mod_add(result[(i + j) % n as usize], mod_mul(a[i], b[j], p), p);
        }
    }
    result
}

/// Multiply two polynomials using NTT (Number Theoretic Transform)
/// 
/// # Arguments
/// 
/// * `a` - First polynomial (as a vector of coefficients).
/// * `b` - Second polynomial (as a vector of coefficients).
/// * `n` - Length of the polynomials and the NTT (must be a power of 2).
/// * `p` - Prime modulus for the operations.
/// * `root` - Primitive root of unity modulo `p`.
///
/// # Returns
/// A vector representing the polynomial product modulo `p`.
pub fn polymul_ntt(a: &[u64], b: &[u64], n: usize, p: u64, root: u64) -> Vec<u64> {
    // Compute n-th root of unity (omega = root^((p - 1) / n) % p)
    let omega = mod_exp(root, (p - 1) / n as u64, p);

    // Step 1: Perform the NTT (forward transform) on both polynomials
    let a_ntt = ntt(a, omega, n, p);
    let b_ntt = ntt(b, omega, n, p);
    
    // Step 2: Perform pointwise multiplication in the NTT domain
    let c_ntt: Vec<u64> = a_ntt
        .iter()
        .zip(b_ntt.iter())
        .map(|(x, y)| mod_mul(*x, *y, p)) // pointwise multiplication
        .collect();

    // Step 3: Apply the inverse NTT to get the result
    let c = intt(&c_ntt, omega, n, p);
    
    c
}

