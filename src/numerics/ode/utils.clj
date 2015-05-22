;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.ode.utils
  (:refer-clojure :rename {+ num+, - num-, * num*, / num-div, == num=})
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [numerics.criteria :as criteria]))

(defrecord ODESolution [t0 u0 f0 h-did h-next t1 u1 f1 error])

(defn ode-initial [t0 u0 f0]
  (->ODESolution t0 u0 f0 0.0 0.0 t0 nil nil nil))

(defn ode-update [ode-solution & {:keys [h-did h-next u1 f1 error]}]
  (assoc ode-solution
    :h-did h-did
    :h-next h-next
    :u1 u1
    :f1 f1
    :error error
    :t1 (num+ (:t0 ode-solution) h-did)))

(defn ode-request-step [ode step] (assoc ode :h-next step))

(defn ode-advance
  "Time step advance moving from t0 to t1. Invalidates all the values at t1"
  [{:keys [t1 u1 f1 h-next]}]
  (->ODESolution t1 u1 f1 0.0 h-next t1 nil nil nil))

(def succeeded-criteria
  (criteria/value-finished
    (fn [{:keys [error]}] (<= error 1.0))
    "Error < Tolerance"))

(def failed-with-zero-step
  (criteria/value-failed
    (fn [{:keys [t0 h-next]}] (num= t0 (num+ t0 h-next)))
    "ODE: cannot control the error, step = 0"))

(defprotocol PODEStepper
  (ode-step
    [method f ode]
    [method f ode tolerance]
    [method f ode tolerance u-scale]))