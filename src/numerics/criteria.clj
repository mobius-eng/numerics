;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.criteria
  (:refer-clojure)
  (:import clojure.lang.Keyword))

(defrecord Criteria [^Keyword status value])

(defn continue [x] (->Criteria ::continue x))
(defn failed [x] (->Criteria ::failed x))
(defn finished [x] (->Criteria ::finished x))

(defn continue? [x] (= (:status x) ::continue))
(defn failed? [x] (= (:status x) ::failed))
(defn finished? [x] (= (:status x) ::finished))

(defn value-finished-criteria [is-finished?]
  (fn [_]
    (fn [x]
      (if (is-finished? x)
        (finished x)
        (continue x)))))

(defn value-failed-criteria [is-failed?]
  (fn [_]
    (fn [x]
      (if (is-failed? x)
        (failed x)
        (continue x)))))

(defn close-enough-criteria [close-enough?]
  (fn [x0]
    (let [x (atom x0)]
      (fn [y]
        (if (close-enough? @x y)
          (finished y)
          (do
            (reset! x y)
            (continue y)))))))

(defn max-count-criteria [max-count]
  (fn [_]
    (let [left (atom max-count)]
      (fn [y]
        (swap! left dec)
        (if (or (zero? @left) (neg? @left))
          (assoc (failed y) :info (format "Exceeded max iteration count %d" max-count))
          (continue y))))))

(defn combine-criteria [& cs]
  (fn [x0]
    (let [cs-f (map #(% x0) cs)]
      (fn [x]
        (reduce (fn [r c]
                  (let [y (c (:value r))]
                    (cond (continue? y) y
                          :else (reduced y))))
                (continue x)
                cs-f)))))
