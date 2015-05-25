;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.utils
  (:refer-clojure :exclude [+ - * / ==])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]))

(defmacro unless
  "if test is false evaluate the body, otherwise returns nil"
  [test & body]
  `(if (not ~test)
     (do
       ~@body)
     nil))

(defn average
  "Computes average of numbers and vectors"
  ([x] x)
  ([x y] (* 0.5 (+ x y)))
  ([x y & more]
   (let [[s n] (reduce (fn [[si i] z] [(+ si z) (inc i)]) [(+ x y) 2] more)]
     (/ s n))))

(defn cube
  "Cube of a number or a vector (elemnt-wise)"
  [x] (* x x x))

(def ^:dynamic *zero-tolerance* 1.0e-12)

(defn almost-zero? [x]
  (< (Math/abs x) *zero-tolerance*))


(defn reduce-indexed-transient [f init coll]
  (let [n (count coll)]
    (loop [result init i 0]
      (cond (reduced? result) @result
            (== i n) result
            :else (recur (f result i (coll i)) (inc i))))))

(defn linspace
  "Like range but: (i) makes a core.matrix vector; (ii) includes the last point if possible"
  ([start end] (linspace start 1 end))
  ([start step end]
    (let [n (quot (+ step (- end start)) step)]
      (compute-matrix [n] #(+ start (* % step))))))