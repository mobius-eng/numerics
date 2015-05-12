;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.ode-test
  (:refer-clojure :exclude [+ - * / ==])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [numerics.ode :as ode]
            [numerics.criteria :as criteria]
            [expectations :refer :all]))

(defn f-scalar-linear [_ x] (- x))

(def y-scalar-linear-4 (ode/rk4-step f-scalar-linear [0.0 1.0] 0.1))
(def y-scalar-linear-45 (ode/rk45ck-attempt f-scalar-linear [0.0 1.0] 0.1 1.0e-5 1.0))
(def y-scalar-linear-answer (Math/exp -0.1))

(expect (< (abs (- y-scalar-linear-4 y-scalar-linear-answer)) 1.0e-5))
(expect (criteria/finished? y-scalar-linear-45))
(expect (< (abs (- y-scalar-linear-answer (:u (:value y-scalar-linear-45)))) 1.0e-5))

(let [y-45 (ode/rk45ck-attempt f-scalar-linear [0.0 1.0] 0.5 1.0e-6 1.0)
      y (:value y-45)]
  (expect (criteria/finished? y-45))
  (expect (< (:error y) 1.0))
  (expect (<= (:h-did y) 0.5))
  (expect (< (abs (- (Math/exp (- (:h-did y))) (:u y))) 1.0e-6)))