;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.ode.rk45ck
  (:refer-clojure :rename {+ num+, - num-, * num*, / num-div, == num=})
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [numerics.utils :refer :all]
            [numerics.fixed-point :refer :all]
            [numerics.criteria :as c]
            [numerics.ode.utils :refer :all]))

(let [a  [0 0.2 0.3 0.6 1.0 0.875]
      c  [(/ 37.0 378) 0 (/ 250 621) (/ 125.0 594) 0 (/ 512.0 1771)]
      c* [(/ 2825.0 27648) 0 (/ 18575.0 48384) (/ 13525.0 55296) (/ 277.0 14336) 0.25]
      dc (mapv - c c*)
      b  [[0 0 0 0 0]
          [0.2 0 0 0 0]
          [(/ 3.0 40) (/ 9.0 40) 0 0 0]
          [0.3 -0.9 1.2 0 0]
          [(/ -11.0 54) 2.5 (/ -70.0 27) (/ 35.0 27) 0]
          [(/ 1631.0 55296) (/ 175.0 512) (/ 575.0 13824) (/ 44275.0 110592) (/ 253.0 4096)]]]
  (defn simple-step
    "Embedded Runge-Kutta 4th and 5th order by Cash-Karp: step and error calculation"
    [f {:keys [t0 u0 h-next]}]
    (let [h (double h-next)
          k (reduce (fn [k i]
                      (conj! k (* h (f (num+ t0 (num* (a i) h))
                                       (reduce-indexed-transient
                                         (fn [u j kj] (+ u (* ((b i) j) kj)))
                                         u0
                                         k)))))
                    (transient [(* h (f t0 u0))])
                    (range 1 6))
          kk (persistent! k)
          ck (map * c kk)
          err (reduce + (map * dc kk))]
      [(reduce + u0 ck) err])))

(defn step
  "Runge-Kutta Cash-Karp method: attempt to advance with step h
  Arguments:
  f        : function f(t,u)
  ode0     : ODESolution initial value
  step     : attempting step
  tol      : realtive tolerance
  u-scale  : scaling factors for u"
  [f ode0 tol u-scale]
  (let [safety 0.9
        p-grow -0.2
        p-shrink -0.25
        errcon 1.89e-4
        improve (fn [{:keys [h-next] :as ode}]
                  (let [[u delta] (simple-step f ode)
                        err (/ (emax (emap (fn [x y] (abs (/ x y))) delta u-scale)) tol)
                        new-h-next (cond (< err errcon) (* 5.0 h-next)
                                         (< err 1.0) (* safety h-next (Math/pow err p-grow))
                                         :else (max (* safety h-next (Math/pow err p-shrink))
                                                    (* 0.1 h-next)))]
                    (ode-update ode
                                :h-did h-next
                                :h-next new-h-next
                                :error err
                                :u1 u)))
        overall-criteria (c/combine succeeded-criteria failed-with-zero-step)
        sol (fixed-point overall-criteria improve ode0)]
    (cond (c/finished? sol) (assoc-in sol
                                      [:value :f1]
                                      (f (->> sol :value :t1) (->> sol :value :u1)))
          :else sol)))


;
;(loop [k [(* h (f t0 u0))] i 1]
;  (cond (== i 6) (let [ck (map * c k)
;                       err (reduce + (map * dc k))]
;                   [(reduce + u0 ck) err])
;        :else (recur (conj k (* h (f (+ t0 (* (a i) h))
;                                     (reduce-indexed-transient
;                                       (fn [u j kj] (+ u (* ((b i) j) kj)))
;                                       u0
;                                       k))))
;                     (inc i))))