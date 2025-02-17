use reikna::totient::totient;
use reikna::factor::quick_factorize;
use std::collections::HashMap;

/// Modular arithmetic functions using i64
fn mod_add(a: i64, b: i64, p: i64) -> i64 {
    (a + b) % p
}

/// Modular multiplication
fn mod_mul(a: i64, b: i64, p: i64) -> i64 {
    (a * b) % p
}

/// Modular exponentiation
/// # Arguments
///
/// * `base` - Base of the exponentiation.
/// * `exp` - Exponent.
/// * `p` - Prime modulus for the operations.
///
/// # Returns
/// The result of the exponentiation modulo `p`.
fn mod_exp(mut base: i64, mut exp: i64, p: i64) -> i64 {
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

/// Extended Euclidean algorithm
fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        (a, 1, 0)  // gcd, x, y
    } else {
        let (gcd, x1, y1) = extended_gcd(b, a % b);
        (gcd, y1, x1 - (a / b) * y1)
    }
}

/// Compute the modular inverse of a modulo modulus
fn mod_inv(a: i64, modulus: i64) -> i64 {
    let (gcd, x, _) = extended_gcd(a, modulus);
    if gcd != 1 {
        panic!("{} and {} are not coprime, no inverse exists", a, modulus);
    }
    (x % modulus + modulus) % modulus  // Ensure a positive result
}

/// Compute n-th root of unity (omega) for p not necessarily prime
/// # Arguments
///
/// * `modulus` - Modulus. n must divide each prime power factor.
/// * `n` - Order of the root of unity.
/// 
/// # Returns
/// The n-th root of unity modulo `modulus`.
pub fn omega(modulus: i64, n: usize) -> i64 {
    let factors = factorize(modulus as i64);
    if factors.len() == 1 {
        let (p, e) = factors.into_iter().next().unwrap();
        let root = primitive_root(p, e); // primitive root mod p
        let grp_size = totient(modulus as u64) as i64;
        assert!(grp_size % n as i64 == 0, "{} does not divide {}", n, grp_size);
        return mod_exp(root, grp_size / n as i64, modulus) // order of mult. group is Euler's totient function
    }
    else {
        return root_of_unity(modulus, n as i64)
    }
}

/// Forward transform using NTT, output bit-reversed
/// # Arguments
///
/// * `a` - Input vector.
/// * `omega` - Primitive root of unity modulo `p`.
/// * `n` - Length of the input vector and the result.
/// * `p` - Prime modulus for the operations.
///
/// # Returns
pub fn ntt(a: &[i64], omega: i64, n: usize, p: i64) -> Vec<i64> {
    let mut result = a.to_vec();
    let mut step = n/2;
	while step > 0 {
		let w_i  = mod_exp(omega, (n/(2*step)).try_into().unwrap(), p);
		for i in (0..n).step_by(2*step) { 
			let mut w = 1;
			for j in 0..step {
				let u = result[i+j];
				let v = result[i+j+step];
				result[i+j] = mod_add(u,v,p);
				result[i+j+step] = mod_mul(mod_add(u,p-v,p),w,p);
				w = mod_mul(w,w_i,p);
			}
		}
		step/=2;
	}
	result
}

/// Inverse transform using INTT, input bit-reversed
/// # Arguments
/// 
/// * `a` - Input vector (bit-reversed).
/// * `omega` - Primitive root of unity modulo `p`.
/// * `n` - Length of the input vector and the result.
/// * `p` - Prime modulus for the operations.
///
/// # Returns
/// A vector representing the inverse NTT of the input vector.
pub fn intt(a: &[i64], omega: i64, n: usize, p: i64) -> Vec<i64> {
    let omega_inv = mod_inv(omega, p);
    let n_inv = mod_inv(n as i64, p);
    let mut result = a.to_vec();
    let mut step = 1;
	while step < n {  
		let w_i = mod_exp(omega_inv, (n/(2*step)).try_into().unwrap(), p);
		for i in (0..n).step_by(2*step) { 
			let mut w = 1;
			for j in 0..step {
				let u = result[i+j];
				let v = mod_mul(result[i+j+step],w,p);
				result[i+j] = mod_add(u,v,p);
				result[i+j+step] = mod_add(u,p-v,p);
				w = mod_mul(w,w_i,p);
			}
		}
		step*=2;
	}
	result
		.iter()
        .map(|x| mod_mul(*x,n_inv,p))
        .collect()
}

/// Naive polynomial multiplication
/// # Arguments
///
/// * `a` - First polynomial (as a vector of coefficients).
/// * `b` - Second polynomial (as a vector of coefficients).
/// * `n` - Length of the polynomials and the result.
/// * `p` - Prime modulus for the operations.
///
/// # Returns
/// A vector representing the polynomial product modulo `p`.
pub fn polymul(a: &Vec<i64>, b: &Vec<i64>, n: i64, p: i64) -> Vec<i64> {
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
pub fn polymul_ntt(a: &[i64], b: &[i64], n: usize, p: i64, omega: i64) -> Vec<i64> {

    // Step 1: Perform the NTT (forward transform) on both polynomials
    let a_ntt = ntt(a, omega, n, p);
    let b_ntt = ntt(b, omega, n, p);
    
    // Step 2: Perform pointwise multiplication in the NTT domain
    let c_ntt: Vec<i64> = a_ntt
        .iter()
        .zip(b_ntt.iter())
        .map(|(x, y)| mod_mul(*x, *y, p)) // pointwise multiplication
        .collect();

    // Step 3: Apply the inverse NTT to get the result
    let c = intt(&c_ntt, omega, n, p);
    
    c
}

/// Compute the prime factorization of `n` (with multiplicities)
/// Uses reikna::quick_factorize internally
/// # Arguments
/// 
/// * `n` - Number to factorize.
/// 
/// # Returns
/// A HashMap with the prime factors of `n` as keys and their multiplicities as values.
fn factorize(n: i64) -> HashMap<i64, u32> {
    let mut factors = HashMap::new();
    for factor in quick_factorize(n as u64) {
        *factors.entry(factor as i64).or_insert(0) += 1;
    }
    factors
}

/// Fast computation of a primitive root mod p^e
pub fn primitive_root(p: i64, e: u32) -> i64 {
    let g = primitive_root_mod_p(p);
    let mut g_lifted = g; // Lift it to p^e
    for _ in 1..e {
        if mod_exp(g_lifted, p-1, p.pow(e)) == 1 {
            g_lifted += p.pow(e - 1);
        }
    }
    g_lifted
}

/// Finds a primitive root modulo a prime p
/// # Arguments
///
/// * `p` - Prime modulus.
///
/// # Returns
/// A primitive root modulo `p`.
fn primitive_root_mod_p(p: i64) -> i64 {
    let phi = p - 1;
    let factors = factorize(phi); // Reusing factorize to get both prime factors and multiplicities
    for g in 2..p {
        // Check if g is a primitive root by checking mod_exp conditions with all prime factors of phi
        if factors.iter().all(|(&q, _)| mod_exp(g, phi / q, p) != 1) {
            return g;
        }
    }
    0 // Should never happen
}

/// the Chinese remainder theorem for two moduli
/// # Arguments
///
/// * `a1` - First residue.
/// * `n1` - First modulus.
/// * `a2` - Second residue.
/// * `n2` - Second modulus.
///
/// # Returns
/// The solution to the system of congruences x = a1 (mod n1) and x = a2 (mod n2).
pub fn crt(a1: i64, n1: i64, a2: i64, n2: i64) -> i64 {
    let n = n1 * n2;
    let m1 = mod_inv(n1, n2); // Inverse of n1 mod n2
    let m2 = mod_inv(n2, n1); // Inverse of n2 mod n1
    let x = (a1 * m2 * n2 + a2 * m1 * n1) % n;
    if x < 0 { x + n } else { x }
}

/// computes an n^th root of unity modulo a composite modulus
/// note we require that an n^th root of unity exists for each multiplicative group modulo p^e
/// use the CRT isomorphism to pull back each n^th root of unity to the composite modulus
/// for the NTT, we require than a 2n^th root of unity exists
/// # Arguments
///
/// * `modulus` - Modulus. n must divide each prime power factor.
/// * `n` - Order of the root of unity.
///
/// # Returns
/// The n-th root of unity modulo `modulus`.
pub fn root_of_unity(modulus: i64, n: i64) -> i64 {
    let factors = factorize(modulus);
    let mut result = 1;
    for (&p, &e) in factors.iter() {
		let omega = omega(p.pow(e), n.try_into().unwrap()); // Find primitive nth root of unity mod p^e
        result = crt(result, modulus / p.pow(e), omega, p.pow(e)); // Combine with the running result using CRT
	}
	result
}

/// ensure the root of unity satisfies sum_{j=0}^{n-1} omega^{jk} = 0 for 1 \le k < n
/// # Arguments
///
/// * `omega` - n-th root of unity.
/// * `n` - Order of the root of unity.
/// * `modulus` - Modulus.
///
/// # Returns
/// True if the root of unity satisfies the condition.
pub fn verify_root_of_unity(omega: i64, n: i64, modulus: i64) -> bool {
    assert!(mod_exp(omega, n, modulus as i64) == 1, "omega is not an n-th root of unity");
    assert!(mod_exp(omega, n/2, modulus as i64) == modulus-1, "omgea^(n/2) != -1 (mod modulus)");
    true
}

