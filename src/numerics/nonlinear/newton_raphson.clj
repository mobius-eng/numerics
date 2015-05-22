;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.nonlinear.newton-raphson
  (:refer-clojure :exclude [+ - * / ==])
  (:require [numerics.fixed-point :refer :all]
            [numerics.criteria :as c]
            [numerics.linear :refer :all]
            [numerics.diff :refer :all]
            [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]))

(defn solve
  [criteria f jac lin-solver x0]
  (let [g (fn [x]
            (let [fx (f x)]
              (- x (lin-solver (jac x) fx fx))))]
    (fixed-point criteria g x0)))

