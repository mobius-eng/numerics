;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.linear.protocols
  (:refer-clojure :exclude [+ - * / ==])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [numerics.utils :refer :all]
            [clojure.core.matrix.linear :as linear]))

(set-current-implementation :vectorz)

;;; Linear operator protocols
(defprotocol PLinOperator
  (apply-operator [A u]
    "Apply linear operator A to u"))

(defprotocol PLinOperator!
  (apply-operator! [A u dest]
    "Apply linear operator A to u, putting result into dest")
  (apply-operator-add! [A u dest]
    "Apply operator A to u and add the result to dest"))

(defprotocol PInvLinOperator
  (apply-inversed [A u] [A u x0]
    "Apply inversed A to u, if provided, use x0 as initial approximation"))

;;; Implementations
(extend-protocol PLinOperator
  java.lang.Number
  (apply-operator [A u] (* A u))
  mikera.matrixx.AMatrix
  (apply-operator [A u] (mmul A u))
  mikera.matrixx.impl.IdentityMatrix
  (apply-operator [_ u] u))

(extend-protocol PLinOperator!
  mikera.matrixx.AMatrix
  (apply-operator! [A u dest]
    (let [[rows cols] (shape A)]
      (loop [i 0]
        (unless (== i rows)
                (mset! dest i 0)
                (loop [j 0]
                  (unless (== j cols)
                          (mset! dest
                                 i
                                 (+ (mget dest i)
                                    (* (mget A i j) (mget u j))))
                          (recur (inc j))))
                (recur (inc i))))))
  (apply-operator-add! [A u dest]
    (let [[rows cols] (shape A)]
      (loop [i 0]
        (unless (== i rows)
                (loop [j 0]
                  (unless (== j cols)
                          (mset! dest
                                 i
                                 (+ (mget dest i)
                                    (* (mget A i j) (mget u j))))
                          (recur (inc j))))
                (recur (inc i))))))
  mikera.matrixx.impl.IdentityMatrix
  (apply-operator! [_ u dest]
    (let [n (ecount u)]
      (loop [i 0]
        (unless (== i n)
                (mset! dest i (mget u i))
                (recur (inc i))))))
  (apply-operator-add! [_ u dest]
    (add! dest u)
    nil))

(extend-protocol PInvLinOperator
  java.lang.Number
  (apply-inversed
    ([A u] (/ u A))
    ([A u _] (/ u A)))
  mikera.matrixx.AMatrix
  (apply-inversed
    ([A u] (linear/solve A u))
    ([A u _] (linear/solve A u)))
  mikera.matrixx.impl.IdentityMatrix
  (apply-inversed
    ([_ u] u)
    ([_ u _] u)))
