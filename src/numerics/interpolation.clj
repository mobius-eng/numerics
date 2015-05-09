(ns numerics.interpolation
  (:refer-clojure :exclude [+ - * / ==])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]))

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

(defn neville-polint
  "Neville\'s algorithm for polynomial interpolation"
  [points]
  (let [xa (mapv first points)
        ya (mapv second points)
        n  (count xa)]
    (fn [x]
      (let [[x-dif x-ind] (reduce-kv (fn [[min-diff ind] i xx]
                                       (let [diff (abs (- x xx))]
                                         (if (< diff min-dff)
                                           [diff i]
                                           [min-diff ind])))
                                     xa)]
        (loop [m 1
               y (ya (dec x-ind))
               c ya
               d ya]
          (if (== m n)
            (do-something)
            (let [[new-c new-d] (loop [i 0
                   new-c []
                   new-d []]
              (if (== i (- n m))
                [new-c new-d]
                (let [ho (- (xa i) x)
                      hp (- (xa (+ i m)) x)
                      w (- (c (inc i)) (d i))
                      den (- ho hp)
                      v (/ w den)]
                  (recur (inc i) (conj new-c (* hp v)) (conj new-d (* ho v))))))])
            ))))))






