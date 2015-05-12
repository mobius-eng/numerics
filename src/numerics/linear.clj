;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.linear
  (:refer-clojure :exclude [+ - * / ==])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [clojure.core.matrix.linear :as linear]))

(set-current-implementation :vectorz)

;;; Linear operator protocols
(defprotocol PLinOperator
  (apply-operator [A u])
  (apply-operator! [A u dest]))

(defprotocol PInvLinOperator
  (apply-inversed [A u] [A u x0]))

;;; Functions: shorthands
(defn op* [A u] (apply-operator A u))

(defn op-inv [A u] (apply-inversed A u))

;;; Implementations
(extend-protocol PLinOperator
  java.lang.Number
  (apply-operator [A u] (* A u))
  (apply-operator! [_ _ _]
    (throw (Exception. "apply-operator!: unsupported for numbers")))
  mikera.matrixx.Matrix
  (apply-operator [A u] (mmul A u))
  (apply-operator! [A u dest]
    (let [[rows cols] (shape A)]
      (loop [i 0]
        (cond (== i rows) nil
              :esle (do
                      (mset! dest i 0)
                      (loop [j 0]
                        (cond (== j cols) :done
                              :else (do
                                      (mset! dest
                                             i
                                             (+ (mget dest i)
                                                (* (mget A i j) (mget u j))))
                                      (recur (inc j)))))
                      (recur (inc i)))))))
  mikera.matrixx.impl.IdentityMatrix
  (apply-operator [A u] u)
  (apply-operator! [A u dest]
    (let [n (ecount u)]
      (loop [i 0]
        (cond (== i n) nil
              :else (do
                      (mset! dest i (mget u i))
                      (recur (inc i))))))))

(extend-protocol PInvLinOperator
  java.lang.Number
  (apply-inversed
    ([A u] (/ u A))
    ([A u _] (/ u A)))
  mikera.matrixx.Matrix
  (apply-inversed
    ([A u] (linear/solve A u))
    ([A u _] (linear/solve A u)))
  mikera.matrixx.impl.IdentityMatrix
  (apply-inversed
    ([A u] u)
    ([A u _] u)))
