(in-ns 'bioclojure.core)

(load "utils")

(def genome "ACGTTGCATGTCGCATGATGCATGAGAGCT")

(defn run-p26
  "Problem 2.6
  	A Branch-and-Bound Algorithm for Cyclopeptide Sequencing"
  [& args]
  (println "Inside problem 2.6 with args:" args))

(defn run-kmers
  "Find all kmers and count them"
  [k & args]
  (println k genome))