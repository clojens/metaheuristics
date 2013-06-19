(ns metaheuristics.core)

;;
;; Common helpers
;;

(defn ^:export scramble
  "The Fisher–Yates shuffle (named after Ronald Fisher and Frank Yates),
  also known as the Knuth shuffle (after Donald Knuth), is an algorithm
  for generating a random permutation of a finite set—in plain terms,
  for randomly shuffling the set."
  [lst]
  (let [items (java.util.ArrayList. lst)]
    scramble (do
               (java.util.Collections/shuffle items)
               (seq items))))

(defn normal
  "Gaussian functions are widely used in statistics where they describe the
  normal distributions and in signal processing where they serve to define
  Gaussian filters."
  [mu sigma]
  (let [r (new java.util.Random)]
    (+ mu
       (* sigma
          (.nextGaussian r)))))
