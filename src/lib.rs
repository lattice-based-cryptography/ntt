pub fn mod_exp(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    let mut result = 1;
    base %= modulus;
    while exp > 0 {
        if exp % 2 == 1 {
            result = (result * base) % modulus;
        }
        base = (base * base) % modulus;
        exp /= 2;
    }
    result
}

pub fn ntt(input: &[u64], root: u64, modulus: u64) -> Vec<u64> {
    let n = input.len();
    let mut output = input.to_vec();
    let mut step = 1;

    while step < n {
        let w = mod_exp(root, (modulus - 1) / (2 * step as u64), modulus);
        for i in (0..n).step_by(2 * step) {
            let mut w_i = 1;
            for j in 0..step {
                let u = output[i + j];
                let v = (output[i + j + step] * w_i) % modulus;
                output[i + j] = (u + v) % modulus;
                output[i + j + step] = (u + modulus - v) % modulus;
                w_i = (w_i * w) % modulus;
            }
        }
        step *= 2;
    }
    output
}

pub fn intt(input: &[u64], root: u64, modulus: u64) -> Vec<u64> {
    let n = input.len() as u64;
    let n_inv = mod_exp(n, modulus - 2, modulus); // Modular multiplicative inverse of n
    println!("n_inv = {}",n_inv);
    let root_inv = mod_exp(root, modulus - 2, modulus);

    let mut output = ntt(input, root_inv, modulus);
    for x in output.iter_mut() {
        *x = (*x * n_inv) % modulus;
    }
    output
}