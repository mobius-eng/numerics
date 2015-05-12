;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.ode
  (:refer-clojure :exclude [+ - * / ==])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [numerics.linear :refer :all]
            [numerics.criteria :as criteria]
            [numerics.fixed-point :refer :all]
            [numerics.newton-raphson :as newton]))

;;; Steppers
(defn rk4-step
  "Fixed step-size 4th order Runge-Kutta method. Solves
  D[u](t) = f(t,u)
  where x is a scalar and y can be either scalar or a vector value.
  Parameters:
  f        : function f(x,y)
  u0       : initial value
  t0       : initial point
  step     : step size"
  [f [t0 u0] step]
  (let [h (double step)
        t_half (double (+ t0 (* 0.5 h)))
        f0 (f t0 u0)
        k1 (* h f0)
        k2 (* h (f t_half (+ u0 (* 0.5 k1))))
        k3 (* h (f t_half (+ u0 (* 0.5 k2))))
        k4 (* h (f (+ t0 h) (+ u0 k3)))]
    (+ u0 (* (/ 6.0) (+ k1 k4)) (* (/ 3.0) (+ k2 k3)))))

(defn reduce-indexed-transient [f init coll]
  (let [n (count coll)]
    (loop [result init i 0]
      (cond (reduced? result) @result
            (== i n) result
            :else (recur (f result i (coll i)) (inc i))))))

(let [a [0 0.2 0.3 0.6 1.0 0.875]
      c [(/ 37.0 378) 0 (/ 250 621) (/ 125.0 594) 0 (/ 512.0 1771)]
      c* [(/ 2825.0 27648) 0 (/ 18575.0 48384) (/ 13525.0 55296) (/ 277.0 14336) 0.25]
      dc (mapv - c c*)
      b  [[0 0 0 0 0]
          [0.2 0 0 0 0]
          [(/ 3.0 40) (/ 9.0 40) 0 0 0]
          [0.3 -0.9 1.2 0 0]
          [(/ -11.0 54) 2.5 (/ -70.0 27) (/ 35.0 27) 0]
          [(/ 1631.0 55296) (/ 175.0 512) (/ 575.0 13824) (/ 44275.0 110592) (/ 253.0 4096)]]]
  (defn rk45ck-step
    "Runge-Kutta 4th and 5th order by Cash-Karp"
    [f [t0 u0] step]
    (let [h (double step)]
      (loop [k [(* h (f t0 u0))] i 1]
        (cond (== i 6) (let [ck (map * c k)
                             err (reduce + (map * dc k))]
                         [(reduce + u0 ck) (Math/sqrt (dot err err))])
              :else (recur (conj k (* h (f (+ t0 (* (a i) h))
                                            (reduce-indexed-transient
                                              (fn [u j kj] (+ u (* ((b i) j) kj)))
                                              u0
                                              k))))
                           (inc i)))))))

(defn rk45ck-attempt
  "Runge-Kutta Cash-Karp method: attempt to advance with step h
  Arguments:
  f        : function f(t,u)
  u0       : initial value
  t0       : initial point
  step     : attempting step
  tol      : realtive tolerance
  u-scale  : scaling factors for u"
  [f [t0 u0] step tol u-scale]
  (let [safety 0.9
        p-grow -0.2
        p-shrink -0.25
        errcon 1.89e-4
        improve (fn [{:keys [h-next]}]
                  (let [[u delta] (rk45ck-step f [t0 u0] h-next)
                        err (/ (emax (emap (fn [x y] (abs (/ x y))) delta u-scale)) tol)]
                    (if (<= err 1.0)
                      (let [new-h-next (if (> err errcon)
                                         (* safety h-next (Math/pow err p-grow))
                                         (* 5.0 h-next))]
                        {:h-did h-next
                         :h-next new-h-next
                         :error err
                         :u u})
                      (let [h-attempt (* safety h-next (Math/pow err p-shrink))
                            new-h-next (if (> h-next 0.0)
                                         (max h-attempt (* 0.1 h-next))
                                         (min h-attempt (* 0.1 h-next)))]
                        {:h-did h-next
                         :h-next new-h-next
                         :error err
                         :u u}))))
        succeeded-criteria (fn [_] (fn [{:keys [error] :as x}]
                                     (if (< error 1.0)
                                       (criteria/finished x)
                                       (criteria/continue x))))
        failed-criteria (fn [_] (fn [{:keys [h-next] :as x}]
                                 (if (== t0 (+ t0 h-next))
                                   (assoc (criteria/failed x)
                                     :info "RK45CK: cannot control the error: step = 0")
                                   (criteria/continue x))))
        overall-criteria (criteria/combine-criteria succeeded-criteria failed-criteria)]
    (fixed-point overall-criteria improve {:h-next step})))
