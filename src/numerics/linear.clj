(ns numerics.linear
  (:refer-clojure :exclude [+ - * / ==])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [clojure.core.matrix.linear :as linear))

;;; Linear operator protocols
(defprotocol PLinOperator
  (apply-operator [A u]))

(defprotocol PInvLinOperator
  (apply-inversed [A u] [A u x0]))

;;; Functions: shorthands
(defn op* [A u] (apply-operator A u))

(defn op-inv [A u] (apply-inversed A u))

;;; Implementations
(extend-protocol PLinOperator
  java.lang.Number
  (apply-operator [A u] (* A u))
  mikera.matrixx.Matrix
  (apply-operator [A u] (mmul A u))
  mikera.matrixx.impl.IdentityMatrix
  (apply-operator [A u] u))

(extend-protocol PInvLinOperator
  java.lang.Number
  (apply-inversed [A u] (/ u A))
  (apply-inversed [A u _] (/ u A))
  mikera.matrixx.Matrix
  (apply-inversed [A u] (linear/solve A u))
  (apply-inversed [A u _] (linear/solve A u))
  mikera.matrixx.impl.IdentityMatrix
  (apply-inversed [A u] u)
  (apply-inversed [A u _] u))
