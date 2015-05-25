# Linear methods

`numerics.linear` exports the following symbols:
- `op*` and `op*!` are generalised matrix-to-vector multiplications: applications of a
  linear operator to a vector. The former is pure, the latter modifies its first argument.
- `op-inv` is the application of an inverted linear operator. In short, it solves @@Ax=b@@
  for @@x@@: @@x = A^-1(b)@@. It is a generic operator: @@A@@ does not have to be a matrix,
  it can also be any kind of linear operator, including `lin-operator` type.
- Functions to work with functional linear operator: `lin-operator`, `lin-operator-full`
  (constructors) and `lin-operator?`. Linear operator consists of two functions: one
  purely maps a vector to vector, another maps the first argument by modifying its second
  argument. It is possible to create just a pure version of a linear operator
  (`lin-operator`), in which case destructive operations, such as `op*!` will wrap pure
  function call.

Functions `op*`, `op*!` and `op-inv` are wrappers over generic functions defined in
protocols `PLinOperator`, `PLinOperator!` and `PInvLinOperator`. Definitions of these
protocols can be found in `numerics.linear.protocols`. These protocols at the moment
are implemented for basic `vectorz` matrices and vectors and a functional linear operator
(see below).

`numerics.linear.bicg-stab` contains the implementation of BiCG stabilised method for
generic linear operators. This method is used by functional linear operators to apply
the inversed operator, as it helps to avoid the instantiation of a matrix.

Namespace `numerics.linear.operator` contains the definition of the function linear
operator (symbols from there are re-exported in `numerics.linear`). 

