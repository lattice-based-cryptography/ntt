# ntt
![CI](https://github.com/jacksonwalters/ntt/actions/workflows/ci.yml/badge.svg)
![MIT License](https://img.shields.io/badge/License-MIT-brightgreen)
[![crates.io](https://img.shields.io/crates/v/ntt.svg)](https://crates.io/crates/ntt)


Implementation of the number theoretic transform (NTT) in Rust.

The array size `n` must be a power of two. 

We allow composite moduli `modulus` for prime powers `p^k` or moduli such that `n` divides `phi(p^e)` for each prime factor p, where `phi` is the Euler totient.
