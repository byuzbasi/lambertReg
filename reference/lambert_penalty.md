# The Fixed Lambert Penalty and Scalar Update

Evaluate the bounded Lambert penalty or its exact scalar least-squares
update with fixed shape \\c=1\\.

## Usage

``` r
lambert_penalty(t, lambda)
lambert_threshold(z, lambda)
```

## Arguments

- t:

  Finite, nonnegative coefficient magnitudes.

- z:

  Finite signed scalar scores.

- lambda:

  A finite positive penalty level.

## Details

Write \\W_0\\ for the principal real Lambert function, characterized by
\\W_0(u)\exp(W_0(u))=u\\ for \\u\geq0\\. For a magnitude \\t\\, let
\\w=W_0(t/\lambda)\\. The penalty is
\$\$p\_\lambda(t)=\frac{\lambda^2}{4} \[e^{2w}(1+2w-2w^2)-1\],\quad
0\leq t\leq e\lambda.\$\$ For \\t\>e\lambda\\, it is constant at
\\\lambda^2(e^2-1)/4\\. The expression is evaluated with
cancellation-resistant arithmetic near zero.

The scalar update minimizes \\(z-b)^2/2+p\_\lambda(\|b\|)\\ over \\b\in
R\\. It is zero when \\\|z\|\leq\lambda\\, equals
\\z\log(\|z\|/\lambda)\\ for \\\lambda\<\|z\|\<e\lambda\\, and equals
\\z\\ when \\\|z\|\geq e\lambda\\. This normalization requires unit
coordinate curvature, obtained by training-only RMS scaling in
`lambert`.

## Value

An unnamed numeric vector of the same length as `t` or `z`, containing
penalty values or scalar updates, respectively. A scalar `lambda` is
shared by all elements; elementwise penalty-level vectors are not
supported. Zero-length input vectors return zero-length vectors.

## See also

[`lambert`](https://byuzbasi.github.io/lambertReg/reference/lambert.md)

## Examples

``` r
lambert_penalty(c(0, 1, exp(1), 4), lambda = 1)
#> [1] 0.0000000 0.9088501 1.5972640 1.5972640
lambert_threshold(c(-4, -1, 0, 1, 2, 4), lambda = 1)
#> [1] -4.000000  0.000000  0.000000  0.000000  1.386294  4.000000
```
