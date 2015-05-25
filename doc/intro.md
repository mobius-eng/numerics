# Introduction to numerics


## `core.matrix` dependency

`numerics` relies heavily on `core.matrix`. To use `numerics` you will also need to
import `core.matrix`. You namespace declaration might look as follows:
```
(ns my-project-namespace
  (:refer-clojure :exclude [+ - * / ==])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [numerics.linear :refer :all]
            [numerics.nonlinear :refer :all]
            ...))
```


TODO: write [great documentation](http://jacobian.org/writing/what-to-write/)
