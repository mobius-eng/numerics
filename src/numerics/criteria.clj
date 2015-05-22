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

(defn value-finished
  ([is-finished? info]
   (fn [_]
     (fn [x]
       (if (is-finished? x)
         (if info
           (update-in (finished x)
                      [:info]
                      (fn [prev-info]
                        (str prev-info "\n" info)))
           (finished x))
         (continue x)))))
  ([is-finished?]
    (value-finished is-finished? nil)))

(defn value-failed
  ([is-failed? info]
   (fn [_]
     (fn [x]
       (if (is-failed? x)
         (if info
           (update-in (failed x)
                      [:info]
                      (fn [prev-info]
                        (str prev-info "\n" info)))
           (failed x))
         (continue x)))))
  ([is-failed?]
    (value-failed is-failed? nil)))

(defn log-value [msg-fn]
  (fn [_]
    (fn [x]
      (println (msg-fn x))
      (continue x))))

(defn close-enough [close-enough?]
  (fn [x0]
    (let [x (atom x0)]
      (fn [y]
        (if (close-enough? @x y)
          (finished y)
          (do
            (reset! x y)
            (continue y)))))))

(defn max-count [max-count]
  (fn [_]
    (let [left (atom max-count)]
      (fn [y]
        (swap! left dec)
        (if (or (zero? @left) (neg? @left))
          (update-in (failed y)
                     [:info]
                     #(str % "\n" (format "Exceeded max iteration count %d"
                                          max-count)))
          (continue y))))))

(defn combine [& cs]
  (fn [x0]
    (let [cs-f (map #(% x0) cs)]
      (fn [x]
        (reduce (fn [r c]
                  (let [y (c (:value r))]
                    (cond (continue? y) y
                          :else (reduced y))))
                (continue x)
                cs-f)))))
