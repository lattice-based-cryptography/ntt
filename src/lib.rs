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

fn mod_inv(a: i64, p: i64) -> i64 {
    mod_exp(a, p - 2, p) // Using Fermat's Little Theorem
}

// Compute n-th root of unity (omega = root^((p - 1) / n) % p)
pub fn omega(root: i64, p: i64, n: usize) -> i64{
    mod_exp(root, (p - 1) / n as i64, p)
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
