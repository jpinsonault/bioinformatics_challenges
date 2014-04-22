(ns bioclojure.core
  (:gen-class))
(def problem-ns *ns*)
(load "p26")
(load "utils")

(defn get-function [arg-string]
  (let [func-string (str "run-" arg-string)]
    (ns-resolve problem-ns (symbol func-string))))

(defn print-usage
  "Print usage"
  []
  (println "Usage: lein run problem_name args..."))

(defn -main
  "Entry point for program. First argument is the name of the function to run (p26, etc),
  the rest are for arguments to pass to that function"
  [& args]
  ; If the function exists
  (if-let [function (get-function (first args))]
    (apply function (rest args))
    (print-usage)))