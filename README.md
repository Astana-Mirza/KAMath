# KAMath
A program for finding Gröbner bases and solving polynomial equations in Galois fields
of characteristic of 2.

Qt and boost are required to build KAMath.

# Description
Solving of systems of polynomial equations can be done with Gröbner bases.
In short, Gröbner basis is an equivalent system of pilynomials. Its benefit is in
the fact that one of equations will contain only one variable or equation `1 = 0`
if the system can't be solved.

System of polynomials can be reduced for simplicity. For example, system
```
x^2y + xy^2z + z^3 = 0
xy = 0
```
can be reduced to
```
z^3 = 0
xy = 0
```
That is done with polynomial division and taking remainder.

## Finding Gröbner bases
KAMath finds Gröbner bases for systems of polynomials over Zn and GF(2^n),
where n is a natural number. For polynomials over Zn, n should be prime,
so Zn would be a field. (Otherwise, divison can be impossible for some elements of Zn.)

For Zn, number n must be specified in appropriate input field.

For GF(2^n), an irreducible polynomial with variable `x` must be specified.
In system of polynomials coefficients from GF(2^n) should be written like this:
`a^k` (or zero), where `a` is a primitive element of this Galios field,
and `k` is some natural power.

Polynomials must be written in appropriate input field starting from new line.

Monomial ordering, used to sort terms in polynomials, can be specified
(default if lexicographical, LEX). Available orderings: LEX, INVLEX, GRLEX and GREVLEX.

Algorithm for finding Gröbner basis can be specified, default is brute-force Buchberger
algorithm. Also, an improved Buchberger algorithm is available (with Buchberger criteria).

The features are available in "System, Z_n" and "System, GF(2^n)" modes.

## Simplifying solution of polynomial equations over GF(2^n)
All elements of field GF(2^n) can be represented in basis `1, a, a^2, ..., a^{n-1}`.
In polynomial over GF(2^n), variables can be decomposed in this basis like this:
```
x = x_0 + a * x_1 + a^2 * x_3 + ... + a^{n-1} * x_{n-1}
```
where all `x_k` are in Z2 field. Sibstituting all variables with these expressions
and composing expressions, multiplied by `a` in every power from `1` to `n-1`,
we will get a system of `n-1` equations over Z2.

KAMath makes this substitution and finds Gröbner basis of system over Z2.
This essentially helps in solution of polynomial equations over GF(2^n).

The feature is available in "Equation, GF(2^n)" mode.
