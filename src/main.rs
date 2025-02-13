mod test;

use reikna::totient::totient;
use ntt::{ntt, intt , polymul, polymul_ntt, mod_exp, mod_inv, verify_root_of_unity};

fn main() {
    let p: i64 = 17; // Prime modulus
    let n: usize = 8;  // Length of the NTT (must be a power of 2)
    let omega = ntt::omega(p, n); // n-th root of unity: root^((p - 1) / n) % p

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
    println!("verify omega = {}", verify_root_of_unity(omega, n as i64, p));
    println!("Polynomial A: {:?}", a);
    println!("Polynomial B: {:?}", b);
    println!("Transformed A: {:?}", a_ntt);
    println!("Transformed B: {:?}", b_ntt);
    println!("Recovered A: {:?}", a_ntt_intt);
    println!("Pointwise Product in NTT Domain: {:?}", c_ntt);
    println!("Resultant Polynomial (c): {:?}", c);
    println!("Standard polynomial mult. result: {:?}", c_std);
    println!("Polynomial multiplication method using NTT: {:?}", c_fast);

    //test the composite modulus case
    let modulus = 51; // Example modulus
    let n = 8;  // Must be a power of 2
    let omega = ntt::omega(modulus, n); // n-th root of unity
    println!("Totient of {}: {}", modulus, totient(modulus as u64));
    println!("omega: {}", omega);
    println!("verify omega = {}", verify_root_of_unity(omega, n as i64, modulus));
    (0..=n).for_each(|i| println!("omega^{}: {}", i, mod_exp(omega, i as i64, modulus)));
    println!("n^-1 = {}", mod_inv(n as i64, modulus));
    let a_ntt = ntt(&a, omega, n, modulus);
    let b_ntt = ntt(&b, omega, n, modulus);
    let a_ntt_intt = intt(&a_ntt, omega, n, modulus);
    let b_ntt_intt = intt(&b_ntt, omega, n, modulus);
    let c_ntt: Vec<i64> = a_ntt
        .iter()
        .zip(b_ntt.iter())
        .map(|(x, y)| (x * y) % modulus)
        .collect();
    let c = intt(&c_ntt, omega, n, modulus);
    let c_std = polymul(&a, &b, n as i64, modulus);
    let c_fast = polymul_ntt(&a, &b, n, modulus, omega);
    println!("A: {:?}", a);
    println!("Transformed A: {:?}", a_ntt);
    println!("Transformed B: {:?}", b_ntt);
    println!("Recovered A: {:?}", a_ntt_intt);
    println!("Recovered B: {:?}", b_ntt_intt);
    println!("Pointwise Product in NTT Domain: {:?}", c_ntt);
    println!("Standard polynomial mult. result: {:?}", c_std);
    println!("Resultant Polynomial (c): {:?}", c);
    println!("Polynomial multiplication method using NTT: {:?}", c_fast);

}
