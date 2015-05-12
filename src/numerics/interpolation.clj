;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.interpolation
  (:refer-clojure :exclude [+ - * / ==])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all])
  (:import numerics.Polint))

(defn lagrange-polint [points]
  (let [xs (mapv first points)
        ys (mapv second points)
        xs-denom (mapv (fn [x y]
                         (/ y
                            (reduce (fn [r xn]
                                       (if (== x xn)
                                         r
                                         (* r (- x xn))))
                                    1
                                    xs)))
                       xs
                       ys)
        numerator-f (fn [i x]
                      (reduce-kv (fn [r j xj]
                                   (if (== i j)
                                     r
                                     (* r (- x xj))))
                                 1
                                 xs))]
    (fn [x]
      (reduce +
              (map-indexed (fn [i dn] (* (numerator-f i x) dn))
                           xs-denom)))))

(defn find-nearest-index
  ([v x [lower upper]]
   (if (== (- upper lower) 1)
     (if (< (abs (- x (v upper))) (abs (- x (v lower))))
       upper
       lower)
     (let [middle (quot (+ lower upper) 2)]
       (cond (< x (v middle)) (recur v x [lower middle])
             (> x (v middle)) (recur v x [middle upper])
             :else middle))))
  ([v x] (find-nearest-index v x [0 (dec (count v))])))

(defn neville-polint
  "Neville's algorithm for polynomial interpolation.
  Returns a function of double"
  [points]
  (let [xa (->> points (mapv first) double-array)
        ya (->> points (mapv second) double-array)]
    (Polint/neville xa ya)))

