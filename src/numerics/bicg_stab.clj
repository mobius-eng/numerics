;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.bicg-stab
  (:refer-clojure :exclude [+ - * / ==])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [numerics.linear :refer :all]
            [numerics.criteria :refer :all]
            [numerics.fixed-point :refer :all]))

(def ^:dynamic *rho-min* 1.0e-12)

(defn solve [criteria A b x0]
  (let [n  (ecount b)
        z (zero-vector n)
        r0 (- b (op* A x0))
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
                             (if (< (Math/abs x) *rho-min*)
                               (throw (Exception. "BiCGStab: |rho| = 0"))
                               x))
                   beta      (* (/ new-rho rho) (/ alpha w))]
               (doto p
                 (add-scaled! v (- w))
                 (scale! beta)
                 (add! r))
               (apply-operator! A p v)
               (let [r0*v (dot r0 v)
                     new-alpha (cond (< (Math/abs r0*v) *rho-min*) 0.0
                                     :else (/ new-rho r0*v))]
                 (add-scaled! r v (- new-alpha))
                 (apply-operator! A r t)
                 (let [t-2 (length-squared t)]
                   (if (< t-2 *rho-min*)
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
    (cond (failed? q) q
          (finished? q) q
          (continue? q) (if (< (length-squared b) *rho-min*)
                          (assoc (finished (assoc y0
                                                  :x
                                                  (zero-vector n)))
                            :info "BiCGStab: RHS is zero")
                          (fixed-point criteria f y0)))))

(defn simple-criteria [error-tolerance max-count]
  (combine-criteria (value-finished-criteria (fn [x]
                                               (let [r2 (dot (:r x) (:r x))]
                                                 (< r2 (* error-tolerance error-tolerance)))))
                    (max-count-criteria max-count)))
