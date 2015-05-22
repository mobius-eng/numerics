;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.diff
  (:refer-clojure :rename {+ num+, - num-, * num*,
                           / num-div, == num=})
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]))

(def ^:dynamic *dx* 1.0e-6)

(defn D [f]
  "Numerical differentiation operator on function f."
  (fn [x]
    (let [dx (double *dx*)
          dx2 (double (num+ dx dx))]
      (cond (number? x) (num-div (num- (f (num+ x dx)) (f (num- x dx))) dx2)
            :else (let [n (long (ecount x))
                        x- (mutable (zero-vector n))
                        x+ (mutable (zero-vector n))
                        ind (array (range 0 n))
                        m (mutable (new-matrix n n))]
                    (loop [i 0]
                      (cond (== i n) m
                            :else (do
                                    (emap! #(if (num= %3 i) (num- %2 dx) %2) x- x ind)
                                    (emap! #(if (num= %3 i) (num+ %2 dx) %2) x+ x ind)
                                    (let [fx- (f x-)
                                          fx+ (mutable (f x+))]
                                      (set-column! m i (scale! (sub! fx+ fx-)
                                                               (/ dx2)))
                                      (recur ^long (inc i)))))))))))
