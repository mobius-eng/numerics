# ODE

Namespace `numerics.ode` contains functions that help to solve ODE: `Du(t) = f(t,u(t))`
with initial value `u(t0) = u0`.

Symbols `rk4`, `rk45ck` and `cn` provide a particular method implementation:
Runge-Kutta 4th order without error control, Runge-Kutta 4th order with error control
(using embedded 5th order formulae due to Cash-Karp) and Crank-Nicolson with error control
using half step with predictor-corrector. Notice, `cn` requiers a spatial Jacobian of `f`.

Main driver is `ode-solve`. See `(doc ode-solve)` for more information.