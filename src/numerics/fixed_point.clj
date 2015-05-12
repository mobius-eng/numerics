;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.fixed-point
  (:refer-clojure)
  (:require [numerics.criteria :refer :all]))

(defn fixed-point [check-status f x0]
  (let [status-f (check-status x0)]
    (loop [x (-> x0 f status-f)]
      (cond (continue? x) (recur (-> x :value f status-f))
            :else x))))

