;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.ode-test
  (:refer-clojure :exclude [+ - * / ==])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [numerics.utils :refer :all]
            [numerics.ode :refer :all]
            [numerics.ode.crank-nicolson :as cn]
            [numerics.ode.utils :refer :all]
            [numerics.criteria :as c]
            [expectations :refer :all]))

(defn f-scalar-linear [_ x] (- x))
(def ode0 (-> (ode-initial 0.0 1.0 -1.0)
              (ode-request-step 0.1)))

(def y-scalar-linear-4 (ode-step rk4 f-scalar-linear ode0))
(def y-scalar-linear-45 (ode-step rk45ck f-scalar-linear ode0 1.0e-5))
(def y-scalar-linear-answer (Math/exp -0.1))

(expect (c/finished? y-scalar-linear-4))
(expect (< (distance (:u1 (:value y-scalar-linear-4)) y-scalar-linear-answer) 1.0e-5))
(expect (c/finished? y-scalar-linear-45))
(expect (< (distance y-scalar-linear-answer (:u1 (:value y-scalar-linear-45))) 1.0e-5))

(let [y-45 (ode-step rk45ck f-scalar-linear (ode-request-step ode0 0.5) 1.0e-6)
      y (:value y-45)]
  (expect (c/finished? y-45))
  (expect (< (:error y) 1.0))
  (expect (<= (:h-did y) 0.5))
  (expect (< (abs (- (Math/exp (- (:h-did y))) (:u1 y))) 1.0e-6)))

(let [y-scalar-cn (ode-step (cn (fn [_ _] -1)) f-scalar-linear ode0 1.0e-7 1.0)
      y (:value y-scalar-cn)]
  (expect (c/finished? y-scalar-cn))
  (expect (< (:error y) 1.0))
  (expect (< (:h-did y) 0.1))
  (expect (> (:h-did y) 0.0))
  (expect (< (abs (- (Math/exp (- (:h-did y))) (:u1 y))) 1.0e-7)))

(let [y-scalar-cn (ode-step (cn (numeric-jac f-scalar-linear)) f-scalar-linear ode0 1.0e-7 1.0)
      y (:value y-scalar-cn)]
  (expect (c/finished? y-scalar-cn))
  (expect (< (:error y) 1.0))
  (expect (< (:h-did y) 0.1))
  (expect (> (:h-did y) 0.0))
  (expect (< (abs (- (Math/exp (- (:h-did y))) (:u1 y))) 1.0e-7)))

(let [times        (inclusive-range 0.0 0.1 1.0)
      y-scalar-sol (ode-solve f-scalar-linear 0.0 1.0 times)
      exact-sol    (emap #(Math/exp (- %)) times)]
  (expect (== (ecount times) 11))
  (expect (== (ecount times) (ecount y-scalar-sol)))
  (expect (< (distance y-scalar-sol exact-sol) 1.0e-6)))

(let [times        (inclusive-range 0.0 0.1 1.0)
      y-scalar-sol (ode-solve f-scalar-linear
                              0.0
                              1.0
                              times
                              (cn (numeric-jac f-scalar-linear))
                              1.0e-9)
      exact-sol    (emap #(Math/exp (- %)) times)]
  (expect (== (ecount times) 11))
  (expect (== (ecount times) (ecount y-scalar-sol)))
  (println y-scalar-sol)
  (println exact-sol)
  (expect (< (emax (- y-scalar-sol exact-sol)) (* 1.0e-8 (emax exact-sol)))))
