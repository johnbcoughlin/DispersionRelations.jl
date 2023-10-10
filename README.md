# DispersionRelations

[![Build Status](https://github.com/johnbcoughlin/DispersionRelations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/johnbcoughlin/DispersionRelations.jl/actions/workflows/CI.yml?query=branch%3Amain)

DispersionRelations.jl is a small package for estimating linear dispersion relations from data.
It is intended for use in computational physics, where it is often necessary to compare the
measured growth rate of some instability or damping process to linear theory.

The package consists of two functions:
- `fit_pure_growth_rate(t, E)`, which fits `E` as a purely growing exponential function of `t`, and
- `fit_complex_frequency(t, E)`, which fits an exponential with complex frequency `ω`.

![Collisional Landau damping example](images/electric_energy.png)
*Estimating collisional Landau damping rates with `fit_complex_frequency`*
