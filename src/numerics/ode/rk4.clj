;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.ode.rk4
  (:refer-clojure :exclude [+ - * / ==])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [numerics.ode.utils :refer :all]
            [numerics.criteria :as c]))

(defn step
  "Fixed step-size 4th order Runge-Kutta method. Solves
  D[u](t) = f(t,u)
  where x is a scalar and y can be either scalar or a vector value.
  Parameters:
  f        : function f(x,y)
  ode      : ODESolution record
  step     : step size"
  [f {:keys [t0 u0 h-next] :as ode0}]
  (let [h (double h-next)
        t_half (double (+ t0 (* 0.5 h)))
        f0 (f t0 u0)
        k1 (* h f0)
        k2 (* h (f t_half (+ u0 (* 0.5 k1))))
        k3 (* h (f t_half (+ u0 (* 0.5 k2))))
        k4 (* h (f (+ t0 h) (+ u0 k3)))
        u1 (+ u0 (* (/ 6.0) (+ k1 k4)) (* (/ 3.0) (+ k2 k3)))]
    (c/finished
      (ode-update ode0
                  :h-next h-next
                  :h-did h-next
                  :u1 u1
                  :f1 (f (+ t0 h-next) u1)))))
