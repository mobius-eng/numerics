(ns numerics.lin-operator
  (:refer-clojure :exclude [+ - * / ==])
  (:require [numerics.fixed-point :refer :all]
            [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [clojure.core.matrix.linear :as linear]
            [numerics.linear :refer :all]
            [numerics.bicg-stab :as bicg]))

(def ^:dynamic *lin-operator-tol* 1.0e-8)
(def ^:dynamic *lin-operator-max-count* 20)

(deftype LinOperator [f]
  clojure.lang.IFn
  (invoke [_ x] (f x))
  (applyTo [this args] (clojure.lang.AFn/applyToHelper this args))
  clojure.core.matrix.protocols.PMatrixScaling
  (scale [m a]
    (LinOperator. (fn [u] (* ((.f m) u) a))))
  (pre-scale [m a]
    (LinOperator. (fn [u] (* a ((.f m) u)))))
  clojure.core.matrix.protocols.PMatrixAdd
  (matrix-add [m a]
    (LinearOperator. (fn [u] (+ ((.f m) u) (op* a u)))))
  (matrix-sub [m a]
    (LinearOperator. (fn [u] (- ((.f m) u) (op* a u)))))
  clojure.core.matrix.protocols.PNegate
  (negate [m]
    (LinOperator. (fn [u] (- ((.f m) u)))))
  PLinOperator
  (apply-operator [A u] ((.f A) u))
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

(defn lin-operator [f] (LinOperator. f))
