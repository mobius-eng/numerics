;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.nonlinear
  (:refer-clojure :exclude [+ - * / ==])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [numerics.nonlinear.newton-raphson :as newton]
            [numerics.criteria :as c]
            [numerics.diff :refer [D]]
            [numerics.linear :refer :all]))

(defn fsolve
  "Solve non-linear equation f(x) = 0. Arguments:
  f                     : function f(x)
  x0                    : initial guess
  jac (optional)        : Jacobian of f(x), it must be a function of x, returning
                          a linear operator
  criteria (optional)   : criteria that determines the behavior of the solver,
                          see numerics.criteria for more details
  lin-solver (optional) : linear solver to be used to solve linear equations arising
                          in the process of solving non-linear equation, defaults to
                          op-inv for a particular type of a Jacobian"
  ([f x0]
   (let [criteria (c/combine (c/close-enough #(< (distance % %2) (* 1.0e-9 (length %))))
                             (c/max-count 50))]
     (fsolve criteria f (D f) op-inv x0)))
  ([f jac x0]
   (let [criteria (c/combine (c/close-enough #(< (distance % %2) (* 1.0e-9 (length %))))
                             (c/max-count 50))]
     (fsolve criteria f jac op-inv x0)))
  ([criteria f jac x0] (fsolve criteria f jac op-inv x0))
  ([criteria f jac lin-solver x0]
   (newton/solve criteria f jac lin-solver x0)))
