Optimization for Shannon diversity in AMPL
================
Derek Corcoran
17/10, 2022

-   [1 Introduction](#1-introduction)
-   [2 Final AMPL model](#2-final-ampl-model)

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

The goal of this repository is to explain and prepare models for the
optimization of geographical areas in order to maximize diversity
(Shannon or other indexes) using AMPL.

# 1 Introduction

Fore every cell $V_i$ in a spatial raster there are $L$ different
possible landuses.

$$\begin{equation*}
\begin{aligned}
& \underset{x}{\text{maximize}}
& & H=-\sum_{i=1}^{SP} P_i\,log_2\,P_i \\
& \text{subject to}
& & f_i(x) \leq b_i, \; i = 1, \ldots, SP.
\end{aligned}
\end{equation*}$$

<details>
<summary>

Show code

</summary>

``` r
mtcars %>%
    ggplot(aes(x = hp, y = mpg)) + geom_point()
```

</details>

# 2 Final AMPL model

``` r
set V;   # vertex set of the spatial graph
set SP; #Set of Species
param L; # Landuse types
param budget >= 0;
param u {SP, V,1..L} ; # node presence absence pred for species SP, in landuse in time L
param c {V,1..L} ; # node carbon content, in landuse in time L
param m {V}; # node cost (all 1) this is for calculating the budget
var y {v in V,l in 1..L} >= 0; # Desition on which landuse to use for node V
var z {v in V} >= 0; # Desition to preserve or not node

maximize ShanonDiv: sum {s in SP, v in V, l in 0..L} m[v]*y[s,v,l]*c[v,l];

subj to MaxBudget {v in V, l in 1..L}:
  budget <= sum {(v) in V} m[v]*z[v];
```
