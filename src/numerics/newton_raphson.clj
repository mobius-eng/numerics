;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.newton-raphson
  (:refer-clojure :exclude [+ - * / ==])
  (:require [numerics.fixed-point :refer :all]
            [numerics.linear :refer :all]
            [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]))

(def ^:dynamic *dx* 1.0e-3)

(defn D [f]
  (fn [x]
    (cond (number? x) (/ (- (f (+ x *dx*)) (f (- x *dx*)))
                         (+ *dx* *dx*))
          :else (let [n (ecount x)
                      m (new-matrix n n)]
                  (loop [i 0]
                    (cond (== i n) m
                          :else (let [x- (mset x
                                               i
                                               (- (mget x i) *dx*))
                                      x+ (mset x
                                               i
                                               (+ (mget x i) *dx*))
                                      fx- (f x-)
                                      fx+ (f x+)]
                                  (set-column! m
                                               i
                                               (/ (- fx+ fx-)
                                                  (+ *dx* *dx*)))
                                  (recur (inc i)))))))))

(defn solve
  ([criteria f jac x0]
   (let [g (fn [x]
             (- x (apply-inversed (jac x) (f x))))]
     (fixed-point criteria g x0)))
  ([criteria f x0] (solve criteria f (D f) x0)))