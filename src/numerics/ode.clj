;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.ode
  (:refer-clojure :rename {+ num+, - num-, * num*, / num-div, == num=})
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [numerics.ode.rk4 :as rk4]
            [numerics.ode.rk45ck :as rk45ck]
            [numerics.ode.crank-nicolson :as cn]
            [numerics.ode.utils :refer :all]
            [numerics.diff :refer [D]]
            [numerics.utils :refer :all]
            [numerics.criteria :as c]))

(def ^{:doc "4th order Runge-Kutta with constant step (no error control)"} rk4
  (reify numerics.ode.utils.PODEStepper
    (ode-step [_ f ode] (rk4/step f ode))
    (ode-step [_ f ode _] (rk4/step f ode))
    (ode-step [_ f ode _ _] (rk4/step f ode))))

(def ^{:doc "Embedded Runge-Kutta method of 4th and 5th oders due to Cash-Karp"} rk45ck
  (reify numerics.ode.utils.PODEStepper
    (ode-step [_ f ode]
      (let [tol 1.0e-7
            u-scale (emap (fn [_] 1.0) (:u0 ode))]
        (rk45ck/step f ode tol u-scale)))
    (ode-step [_ f ode tolerance]
      (let [u-scale (emap (fn [_] 1.0) (:u0 ode))]
        (rk45ck/step f ode tolerance u-scale)))
    (ode-step [_ f ode tolerance u-scale] (rk45ck/step f ode tolerance u-scale))))

(defn cn [jac]
  "Crank-Nicolson method with variable step and error control"
  (reify numerics.ode.utils.PODEStepper
    (ode-step [_ f ode]
      (let [u-scale (emap (constantly 1.0) (:u0 ode))
            tolerance 1.0e-7]
        (cn/step f jac ode tolerance u-scale)))
    (ode-step [_ f ode tolerance]
      (let [u-scale (emap (constantly 1.0) (:u0 ode))]
        (cn/step f jac ode tolerance u-scale)))
    (ode-step [_ f ode tolerance u-scale]
      (cn/step f jac ode tolerance u-scale))))

(defn numeric-jac [f]
  "Returns a Jacobian of function f(t,x): linear operator df(t,x)"
  (fn [t x] ((D (partial f t)) x)))

(defn- times-to-compute [t0 t-out]
  (if (almost-zero? (num- (first t-out) t0))
    (rest t-out)
    (seq t-out)))

(defn- close-target-time? [ode t-target]
  (< (num- t-target (num+ (:t0 ode) (:h-next ode))) 0.0))

(defn ode-solve
  ([f t0 u0 t-out]
    (ode-solve f t0 u0 t-out rk45ck))
  ([f t0 u0 t-out method]
    (let [tolerance 1.0e-7]
      (ode-solve f t0 u0 t-out method tolerance)))
  ([f t0 u0 t-out method local-tolerance]
    (let [u-scale (emap (constantly 1.0) u0)]
      (ode-solve f t0 u0 t-out method local-tolerance u-scale)))
  ([f t0 u0 t-out method local-tolerance u-scale]
    (let [t-out (times-to-compute t0 t-out)
          h-init (double (num- ^double (first t-out) t0))
          ode0 (-> (ode-initial t0 u0 (f t0 u0)) (ode-request-step h-init))]
      (loop [ode ode0 u-out [u0] time-out t-out]
        (if (empty? time-out)
          (array u-out)
          (let [target-t (double (first time-out))
                ode1 (if (close-target-time? ode target-t)
                       (assoc ode :h-next (num- target-t (:t0 ode)))
                       ode)
                new-res (ode-step method f ode1 local-tolerance u-scale)
                new-ode (:value new-res)]
            (if (c/finished? new-res)
              (let [t-achieved (double (:t1 new-ode))]
                (if (almost-zero? (num- t-achieved target-t))
                  (recur (-> new-ode
                             ode-advance
                             (update-in [:h-next]
                                        #(Math/max ^double (:h-next ode) ^double %)))
                         (conj u-out (:u1 new-ode))
                         (rest time-out))
                  (recur (-> new-ode ode-advance) u-out time-out)))
              (Exception. (str "ode-solve: cannot advance in " new-res)))))))))