(ns numerics.bicg-stab
  (:refer-clojure :exclude [+ - * / ==])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [numerics.linear :refer :all]
            [numerics.criteria :refer :all]
            [numerics.fixed-point :refer :all]))

(def ^:dynamic *rho-min* 1.0e-12)
;; (def ^:dynamic *sqrt-rho-min* (Math/sqrt *rho-min*))

(defn solve [criteria A b x0]
  (let [vsquare (fn [v] (dot v v))
        n  (ecount b)
        r0 (v- b (op* A x0))
        y0 {:rho   1.0
            :alpha 1.0
            :w     1.0
            :v     (zero-vector n)
            :p     (zero-vector n)
            :r     r0
            :x     x0}
        f  (fn [{:keys [rho alpha w v p r x]}]
             (let [new-rho   (let [x ^double (dot r0 r)]
                               (if (< (Math/abs x) *rho-min*)
                                 (throw (Exception. "BiCGStab: |rho| = 0"))
                                 x))
                   beta      (* (/ new-rho rho) (/ alpha w))
                   new-p     (+ r (* beta (- p (* w v))))
                   new-v     (op* A new-p)
                   r0*new-v  (dot r0 new-v)
                   new-alpha (if (< (Math/abs r0*new-v) *rho-min*)
                               0.0
                               (/ new-rho r0*new-v))
                   s         (- r (v* new-alpha new-v))
                   t         (op* A s)
                   t-2       (vsquare t)]
               (if (< t-2 *rho-min*)
                 (let [new-x (+ x (* new-alpha new-p))
                       new-r (- b (op* A new-x))]
                   {:rho new-rho :alpha new-alpha :w w :v new-v
                    :p new-p :r new-r :x new-x})
                 (let [new-w (/ (dot t s) t-2)
                       new-x (+ x (* new-alpha new-p) (* new-w s))
                       new-r (- s (* new-w t))]
                   {:rho new-rho :alpha new-alpha :w new-w :v new-v
                    :p new-p :r new-r :x new-x}))))
        q ((criteria y0) y0)]
    (cond (failed? q) q
          (finished? q) q
          (continue? q) (if (< (vsquare b) *rho-min*)
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
