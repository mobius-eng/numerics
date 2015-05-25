;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.linear.operator
  (:refer-clojure :exclude [+ - * / ==])
  (:require [numerics.fixed-point :refer :all]
            [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [clojure.core.matrix.linear :as linear]
            [numerics.linear.protocols :refer :all]
            [numerics.linear.bicg-stab :as bicg]
            [numerics.utils :refer [unless]]))

(def ^:dynamic *lin-operator-tol* 1.0e-8)
(def ^:dynamic *lin-operator-max-count* 20)

(deftype LinOperator [f f!]
  clojure.lang.IFn
  (invoke [_ x] (f x))
  (applyTo [this args] (clojure.lang.AFn/applyToHelper this args))
  clojure.core.matrix.protocols.PMatrixScaling
  (scale [m a]
    (LinOperator. (fn [u] (* ((.f m) u) a))
                  (fn [u dest]
                    ((.f! m) u dest)
                    (scale! dest a)
                    nil)))
  (pre-scale [m a]
    (LinOperator. (fn [u] (* a ((.f m) u)))
                  (fn [u dest]
                    ((.f! m) u dest)
                    (scale! dest a)
                    nil)))
  clojure.core.matrix.protocols.PMatrixAdd
  (matrix-add [m a]
    (LinOperator. (fn [u] (+ ((.f m) u) (apply-operator a u)))
                  (fn [u dest]
                    ((.f! m) u dest)
                    (apply-operator-add! a u dest))))
  (matrix-sub [m a]
    (LinOperator. (fn [u] (- ((.f m) u) (apply-operator a u)))
                  (fn [u dest]
                    ((.f! m) u dest)
                    (sub! dest (apply-operator a u))
                    nil)))
  clojure.core.matrix.protocols.PNegation
  (negate [m]
    (LinOperator. (fn [u] (- ((.f m) u)))
                  (fn [u dest]
                    ((.f! m) u dest)
                    (negate! dest)
                    nil)))
  PLinOperator
  (apply-operator [A u] ((.f A) u))
  PLinOperator!
  (apply-operator! [A u dest]
    (if-let [f! (.f! A)]
      (f! u dest)
      (let [v ((.f A) u)
            n (ecount u)]
        (loop [i 0]
          (unless (== i n)
                  (mset! dest i (mget v i))
                  (recur ^long (inc i)))))))
  PInvLinOperator
  (apply-inversed [A u]
    (bicg/solve (bicg/simple-criteria *lin-operator-tol*
                                      *lin-operator-max-count*)
                A
                u
                u))
  (apply-inversed [A u x0]
    (bicg/solve (bicg/simple-criteria *lin-operator-tol*
                                      *lin-operator-max-count*)
                A
                u
                x0)))

(defn lin-operator [f] (LinOperator. f nil))
(defn lin-operator-full [f f!] (LinOperator. f f!))
(defn lin-operator? [f] (instance? LinOperator f))
