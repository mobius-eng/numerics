;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.linear-test
  (:refer-clojure :exclude [+ - * / ==])
  (:require [numerics.criteria :as c]
            [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [clojure.core.matrix.linear :as linear]
            [numerics.linear :refer :all]
            [numerics.bicg-stab :as bicg]
            [expectations :refer :all]))

(set-current-implementation :vectorz)

;; Linear operation coverage for core.matrix
(let [v (array [1 2 3])
      m (matrix [[1 2 3] [4 5 6] [7 9 9]])
      b (array [14 32 52])
      x (mutable (zero-vector 3))]
  (expect (mmul m v) (op* m v))
  (do
    (apply-operator! m v x)
    (expect (mmul m v) x))
  (expect (linear/solve m b) (apply-inversed m b)))

;; Simple BiCGStab
(let [m (matrix [[1 2 3] [4 5 6] [7 9 9]])
      b (array [14 32 52])
      criteria (bicg/simple-criteria 1.0e-8 3)
      x (bicg/solve criteria m b (array [1 1 1]))
      true-x (linear/solve m b)]
  (expect (c/finished? x))
  (expect (< (length (- (:x (:value x)) true-x)) 1.0e-7)))


(let [m (identity-matrix 10)
      b (array [1 2 3 4 5 6 7 8 9 10])
      x0 (array (repeat 10 7))
      criteria (bicg/simple-criteria 1.0e-8 3)
      x (bicg/solve criteria m b x0)]
  (expect (c/finished? x))
  (expect (< (length (- b (:x (:value x)))) 1.0e-7)))