(in-ns 'bioclojure.core)

(load "utils")

(def genome "ACGTTGCATGTCGCATGATGCATGAGAGCT")

(defn max-value
    "Returns the max value in a map"
    [input-map]
    (apply max (vals input-map)))

(defn run-kmers
  "Find all kmers and count them"
  [k-str & args]
  (let [k (read-string k-str)
        counts-map (frequencies (partition k 1 genome))
        max-count (max-value counts-map)]

    (println 
         (apply str (interpose " " (map 
            #(apply str %)
            (keys
                ; Find all keys in the map who's values are max-count
                (filter 
                    (fn [[k v]] (= v max-count))
                    counts-map))))))))