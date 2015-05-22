;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.ode.crank-nicolson
  (:refer-clojure :exclude [+ - * / ==])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [numerics.linear :refer :all]
            [numerics.nonlinear :refer [fsolve]]
            [numerics.criteria :as c]
            [numerics.ode.utils :refer :all]
            [numerics.fixed-point :refer :all]
            [numerics.utils :refer [average]]))

(defn direct-step [u0 f0 f1 ^double h]
  (+ u0 (* (* 0.5 h) (+ f0 f1))))

(def ^:dynamic *crank-nicolson-solver-criteria*
  (c/combine (c/close-enough
               (fn [x y]
                 (<= (distance x y) (* 1.0e-9 (length x)))))
             (c/max-count 20)))

(defn simple-step [f jac {:keys [t0 u0 f0 h-next]}]
  (let [h-2 (double (* 0.5 h-next))
        t   (double (+ t0 h-next))
        g   (fn [y]
              (- y (* h-2 (f t y)) u0 (* h-2 f0)))
        n   (int (ecount u0))
        dg  (if (== n 1)
              (fn [u] (+ (* (- h-2) (jac t u)) 1.0))
              (fn [u] (+ (* (- h-2) (jac t u)) (identity-matrix (ecount u0)))))]
    (fsolve *crank-nicolson-solver-criteria* g dg u0)))

(defn step [f jac {:keys [t0 u0 f0] :as ode0} tolerance u-scale]
  (let [safety       1.2
        errmin       0.008
        improve      (fn [{:keys [h-next] :as ode}]
                       (let [y (simple-step f jac ode)]
                         (if (c/finished? y)
                           (let [h-2              (* 0.5 h-next)
                                 u1               (:value y)
                                 f1               (f (+ t0 h-next) u1)
                                 f-half-predictor (average f0 f1)
                                 u-half-predictor (direct-step u0 f0 f-half-predictor h-2)
                                 f-half-corrector (f (+ t0 h-2) u-half-predictor)
                                 u-half-corrector (direct-step u0 f0 f-half-corrector h-2)
                                 u1-refined       (direct-step u-half-corrector
                                                                              f-half-corrector
                                                                              f1
                                                                              h-2)
                                 delta            (* (* safety 0.75) (- u1 u1-refined))
                                 err              (/ (emax (emap (fn [x y] (abs (/ x y))) delta u-scale)) tolerance)]
                             (cond (< err errmin) (ode-update ode
                                                              :h-did h-next
                                                              :h-next (* 5.0 h-next)
                                                              :u1 u1
                                                              :f1 f1
                                                              :error err)
                                   (<= err 1.0) (ode-update ode
                                                            :u1 u1
                                                            :h-did h-next
                                                            :h-next (* h-next (Math/pow err -0.3333))
                                                            :error err
                                                            :f1 f1)
                                   :else (ode-update ode
                                                     :u1 u1
                                                     :h-did h-next
                                                     :h-next (min (* h-next (Math/pow err -0.3333))
                                                                  (* 0.1 h-next))
                                                     :error err
                                                     :f1 f1)))
                           (ode-update ode0
                                       :u1 u0
                                       :h-did 0.0
                                       :h-next (* 0.1 h-next)
                                       :error 100.0
                                       :f1 f0))))
        check-status (c/combine succeeded-criteria failed-with-zero-step)]
    (fixed-point check-status improve ode0)))
