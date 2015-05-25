# Nonlinear methods

The only function exported from `numerics.nonlinear` is `fsolve`, which solve nonlinear
equation of the form `f(x)=0`. At the moment `fsolve` uses Newton-Raphson method. If
Jacobian is not provided, it is calculated numerically. To control the solver (i.e.
how many iterations it is allowed to make or what is the tolerance) use `criteria`
parameter. [More on criteria](criteria.md)

Namespace `numerics.nonlinear.newton-raphson` provides Newton-Raphson method
implementation. This method uses `fixed-point` algorithm (from `numerics.fixed-point`)
on a function `g(x)=x - f(x)/Df(x)`, where `Df` is a Jacobian, to solve for `x=g(x)`.

The result of `fsolve` is a variable of a type `Criteria`. To check its status you
can use `finished?`, `failed?` or `continue?` functions. The value can be extracted using
`:value` keyword: `(:value y)`. More details can be found in `numerics.criteria`.
