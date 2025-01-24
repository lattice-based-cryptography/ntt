mod test;

use ntt::{ntt, intt, mod_exp, polymul, polymul_ntt};

fn main() {
    let p: u64 = 17; // Prime modulus
    let root: u64 = 3;     // Primitive root of unity for the modulus
    let n: usize = 8;      // Length of the NTT (must be a power of 2)

    // Compute n-th root of unity: ω = g^((p - 1) / n) % p
    let omega = mod_exp(root, (p - 1) / n as u64, p);

    // Input polynomials (padded to length `n`)
    let mut a = vec![1, 2, 3, 4];
    let mut b = vec![5, 6, 7, 8];
    a.resize(n, 0);
    b.resize(n, 0);

    // Perform the forward NTT
    let a_ntt = ntt(&a, omega, n, p);
    let b_ntt = ntt(&b, omega, n, p);

    // Perform the inverse NTT on the transformed A for verification
    let a_ntt_intt = intt(&a_ntt, omega, n, p);

    // Pointwise multiplication in the NTT domain
    let c_ntt: Vec<u64> = a_ntt
        .iter()
        .zip(b_ntt.iter())
        .map(|(x, y)| (x * y) % p)
        .collect();

    // Inverse NTT to get the polynomial product
    let c = intt(&c_ntt, omega, n, p);

    let c_std = polymul(&a,&b,n as u64,p);
    let c_fast = polymul_ntt(&a,&b,n,p,root);

    // Output the results
    println!("Polynomial A: {:?}", a);
    println!("Polynomial B: {:?}", b);
    println!("Transformed A: {:?}", a_ntt);
    println!("Transformed B: {:?}", b_ntt);
    println!("Recovered A: {:?}", a_ntt_intt);
    println!("Pointwise Product in NTT Domain: {:?}", c_ntt);
    println!("Resultant Polynomial (c): {:?}", c);
    println!("Standard polynomial mult. result: {:?}", c_std);
    println!("Polynomial multiplication method using NTT: {:?}", c_fast);
}

