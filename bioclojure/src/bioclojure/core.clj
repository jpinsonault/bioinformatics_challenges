(ns bioclojure.core
  (:gen-class))

(load "p26")
(load "utils")

(defn call [this & that]
  ((ns-resolve *ns* (symbol this))))

(defn -main
  "Entry point for program. First argument is the name of the function to run (p26, etc),
  the rest are for arguments to pass to that function"
  [& args]
  (call (str "run-" (first args))))
