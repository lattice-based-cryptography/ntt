# ntt
![CI](https://github.com/jacksonwalters/ntt/actions/workflows/ci.yml/badge.svg)
![MIT License](https://img.shields.io/badge/License-MIT-brightgreen)
[![crates.io](https://img.shields.io/crates/v/ntt.svg)](https://crates.io/crates/ntt)


Implementation of the number theoretic transform (NTT) in Rust.

This is a discrete Fourier transform over a finite field of prime order p rather than over the complex numbers.

We also allow the case of the ring Z_{p^2}, the integers modulo p^2. In this case the mult. group has order p^2-p.
