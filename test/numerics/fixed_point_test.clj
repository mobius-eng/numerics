;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.fixed-point-test
  (:refer-clojure :exclude [+ - * / ==])
  (:require [numerics.criteria :as c]
            [numerics.fixed-point :refer :all]
            [numerics.newton-raphson :as newton]
            [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [expectations :refer :all]))

;;; Fixed-point algorithm for scalars
(let [x0 2.0
      status (c/combine-criteria (c/close-enough-criteria (fn [x y] (< (Math/abs (- x y)) 1.0e-8)))
                                 (c/max-count-criteria 50))
      f (fn [x] (Math/sqrt x))
      x (fixed-point status f x0)]
  (expect (c/finished? x))
  (expect (< (Math/abs (- (:value x) 1.0)) 1.0e-7)))


;;; Newton-Raphson for scalars
(let [x0 2.0
      status1 (c/combine-criteria (c/close-enough-criteria (fn [x y] (< (abs (- x y)) 1.0e-8)))
                                  (c/max-count-criteria 30))
      status2 (c/combine-criteria (c/close-enough-criteria (fn [x y] (< (abs (- x y)) 1.0e-8)))
                                  (c/max-count-criteria 10))
      f (fn [x] (- x (Math/sqrt x)))
      df (fn [x] (- 1.0 (/ 0.5 (Math/sqrt x))))
      x1 (newton/solve status1 f x0)
      x2 (newton/solve status2 f df x0)]
  (expect (c/finished? x1))
  (expect (< (Math/abs (- (:value x1) 1.0)) 1.0e-7))
  (expect (c/finished? x2))
  (expect (< (Math/abs (- (:value x2) 1.0)) 1.0e-7)))
