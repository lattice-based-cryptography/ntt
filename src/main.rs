use ntt::{ntt,intt};

fn main() {
    let modulus: u64 = 65536; // 2^16
    let root: u64 = 3; // Example root of unity for 2^16
    let n: usize = 8; // Length of the NTT (power of 2)

    // Input polynomials (padded to length `n`)
    let mut a = vec![1, 2, 3, 4];
    let mut b = vec![5, 6, 7, 8];
    a.resize(n, 0);
    b.resize(n, 0);

    // Perform the forward NTT
    let a_ntt = ntt(&a, root, modulus);
    let b_ntt = ntt(&b, root, modulus);

    // Pointwise multiplication in the NTT domain
    let c_ntt: Vec<u64> = a_ntt
        .iter()
        .zip(b_ntt.iter())
        .map(|(x, y)| (x * y) % modulus)
        .collect();

    // Inverse NTT to get the polynomial product
    let c = intt(&c_ntt, root, modulus);

    // Output the results
    println!("Polynomial A: {:?}", a);
    println!("Polynomial B: {:?}", b);
    println!("Transformed A: {:?}", a_ntt);
    println!("Transformed B: {:?}", b_ntt);
    println!("Pointwise Product in NTT Domain: {:?}", c_ntt);
    println!("Resultant Polynomial (c): {:?}", c);
}
