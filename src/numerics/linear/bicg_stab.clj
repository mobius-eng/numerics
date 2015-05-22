;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.linear.bicg-stab
  (:refer-clojure :exclude [+ - * / ==])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [numerics.linear.protocols :refer :all]
            [numerics.criteria :as criteria]
            [numerics.fixed-point :refer :all]
            [numerics.utils :refer [almost-zero?]]))

(defn solve [criteria A b x0]
  (let [n  (ecount b)
        z (zero-vector n)
        r0 (- b (apply-operator A x0))
        y0 {:rho   1.0
            :alpha 1.0
            :w     1.0
            :v     (mutable z)
            :p     (mutable z)
            :t     (mutable z)
            :r     (mutable r0)
            :x     (mutable x0)}
        f  (fn [{:keys [rho alpha w v p r t x]}]
             (let [new-rho (let [x (double (dot r0 r))]
                             (if (almost-zero? x)
                               (throw (Exception. "BiCGStab: |rho| = 0"))
                               x))
                   beta      (* (/ new-rho rho) (/ alpha w))]
               (doto p
                 (add-scaled! v (- w))
                 (scale! beta)
                 (add! r))
               (apply-operator! A p v)
               (let [r0*v (dot r0 v)
                     new-alpha (cond (almost-zero? r0*v) 0.0
                                     :else (/ new-rho r0*v))]
                 (add-scaled! r v (- new-alpha))
                 (apply-operator! A r t)
                 (let [t-2 (length-squared t)]
                   (if (almost-zero? t-2)
                     (do
                       (add-scaled! x p new-alpha)
                       (apply-operator! A x r)
                       (sub! r b)
                       (negate! r)
                       {:rho new-rho :alpha new-alpha :w w :v v :p p :r r :x x :t t})
                     (let [new-w (/ (dot t r) t-2)]
                       (doto x
                         (add-scaled! p new-alpha)
                         (add-scaled! r new-w))
                       (add-scaled! r t (- new-w))
                       {:rho new-rho :alpha new-alpha :w new-w :v v :p p :r r :x x :t t}))))))
        q ((criteria y0) y0)]
    (cond (criteria/failed? q) q
          (criteria/finished? q) q
          (criteria/continue? q) (if (almost-zero? (length-squared b))
                          (assoc (criteria/finished (assoc y0
                                                      :x
                                                      (zero-vector n)))
                            :info "BiCGStab: RHS is zero")
                          (fixed-point criteria f y0)))))

(defn simple-criteria [error-tolerance max-count]
  (criteria/combine (criteria/value-finished (fn [x]
                                                (let [r2 (dot (:r x) (:r x))]
                                                  (< r2 (* error-tolerance error-tolerance)))))
                    (criteria/max-count max-count)))
