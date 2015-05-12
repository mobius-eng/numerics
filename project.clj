(defproject numerics "0.1.0-SNAPSHOT"
  :description "Numeric algrithms for Clojure"
  :url "https://github.com/mobius-eng/numerics"
  :license {:name "LGPL3"
            :url "http://www.gnu.org/licenses/lgpl-3.0.txt"}
  :plugins [[lein-gorilla "0.3.4"]
            [lein-autoexpect "1.4.2"]]
  :dependencies [[org.clojure/clojure "1.6.0"]
                 [net.mikera/core.matrix "0.34.0"]
                 [net.mikera/vectorz-clj "0.29.0"]]
  :profiles {:dev {:dependencies [[expectations "2.0.9"]]}}
  :java-source-paths ["java"]
  :source-paths ["src"]
  :prep-tasks [["compile" "numerics.criteria"]
               "javac" "compile"])
