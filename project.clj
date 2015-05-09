(defproject numerics "0.1.0-SNAPSHOT"
  :description "Numeric algrithms in Clojure"
  :url "http://example.com/FIXME"
  :license {:name "LGPL3"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :plugins [[lein-gorilla "0.3.4"]
            [lein-autoexpect "1.4.2"]]
  :dependencies [[org.clojure/clojure "1.6.0"]
                 [net.mikera/core.matrix "0.34.0"]
                 [net.mikera/vectorz-clj "0.29.0"]]
  :profiles {:dev {:dependencies [[expectations "2.0.9"]]}})
