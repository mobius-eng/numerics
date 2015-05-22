;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.ode.midpoint
  (:refer-clojure :rename {+ num+, - num-, * num*, / num-div, == num=})
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [numerics.utils :refer :all]
            [numerics.fixed-point :refer :all]
            [numerics.criteria :as criteria]))


(defn modified-midpoint
  "Modified midpoint method"
  [f f0 y0 x0 step num-steps]
  (let [h (/ step num-steps)
        z0 y0]
    (loop [z-1 z0
           z   (+ z0 (* h f0))
           m 1]
      (cond (== m num-steps) (* 0.5 (+ z z-1 (* h (f (+ x0 step) z))))
            :else (recur z
                         (+ z-1 (* (* 2.0 h) (f (+ x0 (* h m)) z)))
                         (inc m))))))