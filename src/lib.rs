use reikna::totient::totient;
use reikna::factor::quick_factorize;
use std::collections::HashMap;

fn gcd(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        let temp = b;
        b = a % b;
        a = temp;
    }
    a.abs()
}

/// Compute LCM of two numbers.
fn lcm(a: i64, b: i64) -> i64 {
    (a * b) / gcd(a, b)
}

// Modular arithmetic functions using i64
fn mod_add(a: i64, b: i64, p: i64) -> i64 {
    (a + b) % p
}

fn mod_mul(a: i64, b: i64, p: i64) -> i64 {
    (a * b) % p
}

pub fn mod_exp(mut base: i64, mut exp: i64, p: i64) -> i64 {
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

//compute the modular inverse of a modulo p using Fermat's little theorem, p not necessarily prime
fn mod_inv(a: i64, p: i64) -> i64 {
    assert!(gcd(a, p) == 1, "{} and {} are not coprime", a, p);
    mod_exp(a, totient(p as u64) as i64 - 1, p) // order of mult. group is Euler's totient function
}

// Compute n-th root of unity (omega) for p not necessarily prime
pub fn omega(root: i64, p: i64, n: usize) -> i64 {
    let grp_size = totient(p as u64) as i64;
    assert!(grp_size % n as i64 == 0, "{} does not divide {}", n, grp_size);
    mod_exp(root, grp_size / n as i64, p) // order of mult. group is Euler's totient function
}

// Forward transform using NTT, output bit-reversed 
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

// Inverse transform using INTT, input bit-reversed 
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

// Naive polynomial multiplication
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

/// Compute the prime factorization of `n` (with multiplicities).
fn factorize(n: i64) -> HashMap<i64, u32> {
    let mut factors = HashMap::new();
    for factor in quick_factorize(n as u64) {
        *factors.entry(factor as i64).or_insert(0) += 1;
    }
    factors
}

/// Fast computation of a primitive root mod p^e
pub fn primitive_root(p: i64, e: u32) -> i64 {

    // Find a primitive root mod p
    let g = find_primitive_root_mod_p(p);

    // Lift it to p^e
    let mut g_lifted = g;
    for _ in 1..e {
        if g_lifted.pow((p - 1) as u32) % p.pow(e) == 1 {
            g_lifted += p.pow(e - 1);
        }
    }
    g_lifted
}

/// Finds a primitive root modulo a prime p
fn find_primitive_root_mod_p(p: i64) -> i64 {
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

/// Finds an element in (Z/NZ)* whose order is divisible by `n`.
pub fn find_cyclic_subgroup(modulus: i64, n: i64) -> (i64, i64) {
    // Check if n is a power of 2
    if n == 0 || (n & (n - 1)) != 0 {
        panic!("n must be a power of 2");
    }

    // Factorize modulus (assuming a function exists)
    let factors = factorize(modulus);
    let mut generators = Vec::new();
    let mut orders = Vec::new();

    // Loop through factors to find generators and orders
    for (&p, &e) in &factors {
        let phi = (p - 1) * p.pow(e - 1); // Euler's totient function
        let g = primitive_root(p, e);
        generators.push(g);
        orders.push(phi);
    }

    let mut chosen_element = 1;
    let mut chosen_order = 1;

    // Loop through generators and orders to find element with required order
    for (&g, &k) in generators.iter().zip(orders.iter()) {
        // Calculate required order
        let required_order = lcm(k, n); // Least common multiple
        let exponent = required_order / gcd(k, required_order); // Adjust exponent
        chosen_element = (chosen_element * mod_exp(g, exponent, modulus)) % modulus; // mod_exp computes power mod modulus
        chosen_order = lcm(chosen_order, k / gcd(k, n)); // Adjust chosen order
    }

    // Assert the order is divisible by n
    assert_eq!(chosen_order % n, 0);

    (chosen_element, chosen_order)
}

