;
; Copyright (c) 2015 Alexey Cherkaev.
; This file is licensed under Lesser General Public License v.3
; The text of the license can be found at http://www.gnu.org/licenses/lgpl-3.0.txt
;

(ns numerics.linear
  (:refer-clojure :exclude [+ - * / ==])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [numerics.linear.protocols :refer :all]
            [numerics.linear.operator]
            [potemkin :refer [import-vars]]))

(defn op*
  "Apply linear operator A to u"
  [A u] (apply-operator A u))

(defn op*!
  "Apply A to u, putting the result into dest"
  [dest A u] (apply-operator! A u dest))

(defn op-inv
  ([A u] (apply-inversed A u))
  ([A u x0] (apply-inversed A u x0)))

(import-vars
  [numerics.linear.operator
   lin-operator
   lin-operator-full
   lin-operator?])