# ntt
![CI](https://github.com/jacksonwalters/ntt/actions/workflows/ci.yml/badge.svg)
![MIT License](https://img.shields.io/badge/License-MIT-brightgreen)
[![crates.io](https://img.shields.io/crates/v/ntt.svg)](https://crates.io/crates/ntt)


Implementation of the number theoretic transform (NTT) in Rust.

This is a discrete Fourier transform over a finite field of prime order p rather than over the complex numbers.

We generally allow the case of the ring Z/NZ where $N = p^k$ or $2 p^k$. In this case the multiplicative group is cyclic and the order is given by the Euler totient function, $\phi(p^k) = \phi(2*p^k) = p^k(p-1)$. 

The array size is `n` and must divide $\phi(p^k)$. If `n` is a power of 2 and $p > 2$, then one must have $n|p-1$.
