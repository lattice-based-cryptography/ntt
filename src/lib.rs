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
