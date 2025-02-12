mod test;

use reikna::totient::totient;
use ntt::{omega, ntt, intt , polymul, polymul_ntt, find_cyclic_subgroup, primitive_root, mod_exp, crt};

fn main() {
    let p: i64 = 17; // Prime modulus
    let root: i64 = 3; // Primitive root of unity for the modulus
    let n: usize = 8;  // Length of the NTT (must be a power of 2)
    let omega = omega(root, p, n); // n-th root of unity: root^((p - 1) / n) % p

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
    let c_ntt: Vec<i64> = a_ntt
        .iter()
        .zip(b_ntt.iter())
        .map(|(x, y)| (x * y) % p)
        .collect();

    // Inverse NTT to get the polynomial product
    let c = intt(&c_ntt, omega, n, p);

    let c_std = polymul(&a, &b, n as i64, p);
    let c_fast = polymul_ntt(&a, &b, n, p, omega);

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

    let modulus = 45; // Example modulus
    let n = 4;  // Must be a power of 2
    let (g, g_order) = find_cyclic_subgroup(modulus, n);
    let omega = mod_exp(g, g_order / n, modulus);
    let root = primitive_root(23, 2);
    println!("Totient of {}: {}", modulus, totient(modulus as u64));
    println!("Primitive root: {}", root);
    println!("(g, order) = {}, {}", g, g_order);
    println!("omega: {}", omega);
    println!("omega^n: {}", mod_exp(omega, n, modulus));

    let x1 = crt(2, 3, 3, 5);
    println!("Expected: 8, Computed: {}", x1);
    let x3 = crt(4, 7, 5, 11);
    println!("Expected: 60, Computed: {}", x3);
}
