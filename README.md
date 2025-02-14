# ntt
![CI](https://github.com/jacksonwalters/ntt/actions/workflows/ci.yml/badge.svg)
![MIT License](https://img.shields.io/badge/License-MIT-brightgreen)
[![crates.io](https://img.shields.io/crates/v/ntt.svg)](https://crates.io/crates/ntt)


Implementation of the number theoretic transform (NTT) in Rust.

The NTT is a DFT over the ring Z/mZ. We use a fast divide-and conquer algorithm. The array size `n` must be a power of two. 

We allow composite moduli as long as `n` divides `phi(p^e)` for each prime factor p of the modulus, where `phi` is the Euler totient.
