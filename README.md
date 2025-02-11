# ntt
![CI](https://github.com/jacksonwalters/ntt/actions/workflows/ci.yml/badge.svg)
![MIT License](https://img.shields.io/badge/License-MIT-brightgreen)
[![crates.io](https://img.shields.io/crates/v/ntt.svg)](https://crates.io/crates/ntt)


Implementation of the number theoretic transform (NTT) in Rust.

This is a discrete Fourier transform over a finite field of prime order p rather than over the complex numbers.

We allow the case of the ring Z/NZ where $N = p^k$. In this case the multiplicative group is cyclic and the order is given by the Euler totient function, $\phi(p^k) = p^k(p-1)$.

The array size is `n` and must be a power of two for the divide-and-conquer algorithm to work, and `n` must also divide $\phi(p^k)$ to have an `n`th root of unity `omega`. This is equivalent to $n|p-1$ for $p > 2$.

Note if $N=2p^k$, the multiplicative group is still cyclic and $\phi(2p^k) = p^k(p-1)$, but $gcd(n,N)=2$, so `n` is not invertible modulo $N$. 
