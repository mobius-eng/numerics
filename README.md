# numerics

Numerical analysis in Clojure. The library's main focus is on small size problems: with
vector sizes not larger than 100.

## Motivation

I needed a numerical analysis library for my PhD project. While there are
core.matrix and Incanter, they main focus is on working with data and machine learning,
whereas basic numerical tools are absent:

+ Iterative linear solvers
+ Nonlinear solvers
+ ODE solvers
+ Integration methods
+ etc.

This library is an attempt to close this gap. It uses core.matrix as a foundation for all
vector operations.

## Implementation features

Some algorithms are implemented in Java if there is clear performance benefit
(so far, Neville's polynomial interpolation), but if there is no performance penalty or it
is rather small, it is implemented in Clojure directly (most notably, fixed point
algorithm).

At the moment implementation depends quite deeply on vectorz, as the code matures, I hope
to make it more abstract and to be able to use, for example, clatrix as well (although at the
time of writing its development seemed to be stalled).

## What is implemented so far

+ BiCG stabilised method for linear systems, accepting both matrix and
  linear operator (function) as an input.
+ Neville's algorithm for polynomial interpolation.
+ Fixed point algorithm as a basis for other iterative methods.
+ Newton-Raphson algorithm for solving nonlinear equations (including vector equations).
+ ODE: so far in a rudimentary form:
    - Classic fixed step Runge-Kutta of 4th order.
    - Variable step embedded Runge-Kutta of 4-5th order with Cash-Karp coefficients.
    - Midpoint method: extension to Bulirsch-Stoer is coming.
    - Crank-Nicolson: 2nd order implicit method with step control

## Nearest plans

+ ODE: add Bulirsch-Stoer method.
+ ODE: add methods for stiff problems.
+ Add integration module: Simpson method and Gaussian quadratures.
+ Conditional numbers and eigenvalues-eigenvectors for matrices.
+ Optimisation: simulated annealing method.

## Status: not ready yet, proceed with caution

This is not even alpha release yet. The test suite is (very) incomplete.
API may change without version change at least until the project matures enough to reach
version 0.0.5!

## License

Copyright © 2015 Alexey Cherkaev

Distributed under the Lesser General Public License v3.
