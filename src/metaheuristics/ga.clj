(ns metaheuristics.ga
  "In a genetic algorithm, a population of candidate solutions to an
  optimization problem is evolved toward better solutions.

  A typical genetic algorithm requires: a genetic representation of the
  solution domain, a fitness function to evaluate the solution domain.

  Once the genetic representation and the fitness function are defined, a GA
  proceeds to initialize a population of solutions and then to improve it
  through repetitive application of the mutation, crossover, inversion and
  selection operators.

  Individual solutions are selected through a fitness-based process, where
  fitter solutions are typically more likely to be selected.

  Fitness (often denoted in population genetics models) is a central
  idea in evolutionary theory. "
  (:require [metaheuristics.core :refer [normal scramble]])
  (:use [clojure.set :only (intersection difference)]
        [clojure.contrib.math]
        [metaheuristics.testfunctions]))

;;
;; Helper functions
;;

(defn- garand-int
  "Generates a list of n random integers below with a maximum."
  [max n]
  (map (fn [_] (int (* max (rand)))) (range n)))

(defn- euclidean
  "Method alength works on Java arrays such as a String[] or Integer[].
  Expects these as input arguments for the function to calculate greatest
  common divisor (gcd)."
  [v1 v2]
  (Math/sqrt (areduce v1 i ret 0
	   (+ ret (Math/pow (- (aget v1 i)
                         (aget v2 i))
                      2)))))


;;
;; Structural
;;

(defstruct individual
  "Natural selection within an individual requires minimal the
  properties as outlined here."
  :tag :chromosome :steps :fitness)

(defstruct population
  "A population is made up of individuals."
  :poplist)

;;
;; Generation / initialization
;;

(defn- generate-chromosome
  "Generate a chromosome strand of n bits long."
  [n bits]
  (let [max (int (- (Math/pow 2 bits) 1))]
    (garand-int max n)))

;(doall (generate-chromosome 22 8))

(defn- init-individual
  "Intializes an individiual by giving it identity (hash-map with tag)
  who is made up of 22 random chromosomes. (??)"
  [n bits]
  (struct individual 0 (generate-chromosome n bits)
	  (map #(* 10 %) (take 22 (repeatedly rand))) 0))

;(init-individual 22 8)

(defn- init-population
  "Initialize a population by populating a list of individuals."
  [n dim bits]
  (let [poplist (map (fn [_] (init-individual dim bits)) (range n))]
    (struct population poplist)))

;;
;; Transformation
;;

(defn chromo-to-phenotype
  "Phenotypes result from the expression of an organism's genes as well as the
  influence of environmental factors and the interactions between the two."
  [chromosome]
  (let [sum (reduce + chromosome)]
    (double-array (for [gene chromosome] (/ gene sum)))))


;;
;; Strategies
;;

(defn- share
  "TODO DocString"
  [dist sigma alpha]
  (if (<= dist sigma)
    (- 1 (Math/pow (/ dist sigma) alpha))
    0))


(defn- fitness-sharing
  "http://en.wikipedia.org/wiki/Fitness_(biology)
  There are two commonly used measures of fitness; absolute fitness and
  relative fitness."
  [ind popu]
  (let [sigma 100 ;; 100 seems to be fine
        alpha 1
        dist-sum (reduce
                  +
                  (for [other (:poplist popu)]
                    (share (euclidean
                            (double-array (:chromosome ind))
                            (double-array (:chromosome other)))
                           sigma alpha)))]
    (assoc ind :fitness (/ (:fitness ind) dist-sum))))

;; (def foo (init-population 10 22 8))
;; (fitness-sharing (init-individual 22 8) foo)


;;
;; Evaluate process
;;

(defn- evaluate-individual
  "Evaluate an individual on basis of their fitness."
  [ind fitness]
  (let [w-pheno (chromo-to-phenotype (:chromosome ind))
        fitval  (fitness w-pheno)]
    (assoc ind :fitness fitval)))

(defn- evaluate-all
  "Evaluates performance of the population in terms of fitness.
  ((The average fitness of the whole population is the fitness of
  each genotype multiplied by its frequency: this is called mean fitness.))"
  [popu fitness fs?]
  (let [agentlist (for [ind (:poplist popu)
                        :when (= (:tag ind) 1)]
                    (agent ind))]
    (dorun (map #(send %1 evaluate-individual fitness) agentlist))
    (apply await agentlist)

    ;; fitness sharing
    (if fs? (do
	      (dorun (map #(send %1 fitness-sharing popu) agentlist))
	      (apply await agentlist)))

    (let [children (for [agent agentlist] @agent)
          parents  (for [p (:poplist popu)
                         :when (= (:tag p) 0)] p)]

      (assoc popu :poplist (concat parents children)))))


(defn- evaluate-all-firstrun
  "Evaluate all on the first run."
  [popu fitness]
  (let [agentlist (for [ind (:poplist popu)] (agent ind))]
    (dorun (map #(send %1 evaluate-individual fitness) agentlist))
    (apply await agentlist)
    (assoc popu :poplist (for [agent agentlist] @agent))))


;;
;; Genetic recombination (synapsis)
;;

(defn- gene-crossover
  "Genetic crossover is the exchange of genes."
  [gene1 gene2 len]
  (let [mask (int (- (Math/pow 2 len) 1))
        last-g1 (bit-and gene1 mask)
        last-g2 (bit-and gene2 mask)
        new-g1 (bit-or
                (bit-shift-left (bit-shift-right gene1 len) len) last-g2)
        new-g2 (bit-or
                (bit-shift-left (bit-shift-right gene2 len) len) last-g1)]
    (list new-g1 new-g2)))

;(gene-crossover 8 3 6)


(defn- crossover
  "Chromosomal crossover is the exchange of genetic material between homologous
  chromosomes that results in recombinant chromosomes."
  [chromo1 chromo2 bits]
  (let [split-pos (nth (range 1 bits) (rand-int (- bits 1)))
        split-length (- bits split-pos)
        result (map #(gene-crossover %1 %2 split-length) chromo1 chromo2)
        c1-new (map #(first %1) result)
        c2-new (map #(second %1) result)]
    (list c1-new c2-new)))

;(crossover '(255 12 18 238 210 199 88) '(4 8 15 16 23 42 108) 8)


(defn- mutation
  "Mutation can result in several different types of change in sequences.
  Mutations in genes can either have no effect, alter the product of a gene,
  or prevent the gene from functioning properly or completely.

  TODO: One study on genetic variations between different species of Drosophila
  suggests that if a mutation changes a protein produced by a gene, the result
  is likely to be harmful, with an estimated 70 percent of amino acid
  polymorphisms having damaging effects, and the remainder being either neutral
  or weakly beneficial.

  Due to the damaging effects that mutations can have on genes, organisms have
  mechanisms such as DNA repair to prevent or correct mutations."

  ;; Note: disable non-used percentage param
  [chromosome bits #_percentage]
  (for [gene chromosome]
    (loop [g gene b 0]
      (if (= b bits)
        g
        (let [new-g (if (< (rand) 0.25) (bit-flip g b) g)]
          (recur new-g (inc b)))))))

;(mutation '(34 21 57 56 11 10) 8)


(defn- do-offspring
  "TODO DocString"
  [acc parents-pair]
  (let [percentage 0.3
        p1c (:chromosome (first parents-pair))
        p2c (:chromosome (second parents-pair))
        children (crossover p1c p2c 8)
        child1 {:tag 1
                :chromosome (mutation (first children)
                                      8 percentage)
                :steps (list)
                :fitness 0}
        child2 {:tag 1
                :chromosome (mutation (second children)
                                      8 percentage)
                :steps (list)
                :fitness 0}
        old-pop (:poplist acc)]
    (assoc acc :poplist (conj old-pop child1 child2))))


(defn- adapted-crossover
  "TODO DocString"
  [p1 p2]
  (let [dim (count (:chromosome p1))
        dim2 (/ dim 2)
        p1ch (:chromosome p1)
        p1st (:steps p1)
        p2ch (:chromosome p2)
        p2st (:steps p2)
        childc (map #(if (= (int (* 2 (rand))) 0) %1 %2) p1ch p2ch)
        childs (map #(if (= (int (* 2 (rand))) 0) %1 %2) p1st p2st)]
    (list childc childs)))


(defn- adapted-mutation
  "TODO DocString"
  [chromosome steps]
  (let [dim 22
        bound 0.5
        tau (/ 1 (sqrt (* 2 (sqrt dim))))
        tauprime-rand (* (/ 1 (sqrt (* 2 dim))) (normal 0 1))
        new-steps-raw (for [s steps]
                        (* s (expt (Math/E)
                                   (+ tauprime-rand (* tau (normal 0 1))))))
        new-steps (map #(if (< % bound) bound %) new-steps-raw)
        new-pos (map #(+ %1 (* %2 (normal 0 1))) chromosome new-steps)]
    (list new-pos new-steps)))


(defn- do-offspring-adapted
  "TODO DocString"
  [acc parents-pair]
  (let [p1 (first parents-pair)
        p2 (second parents-pair)
        crossed (adapted-crossover p1 p2)
        mutated (adapted-mutation (first crossed) (second crossed))
        child {:tag 1
               :chromosome (first mutated)
               :steps (second mutated)
               :fitness 0}
        old-pop (:poplist acc)]
    (assoc acc :poplist (conj old-pop child))))


(defn- generate-offspring
  "TODO DocString"
  [acc parents adapted?]
  (if adapted?
    (reduce do-offspring-adapted acc parents)
    (reduce do-offspring acc parents)))


(defn- survivor-selection
  "TODO DocString"
  [popu ftype percent popsize]
  (let [sorted-all (sort-by :fitness ftype (:poplist popu))
        size-all (count sorted-all)
        size-top (* percent size-all)
        size-rest (- popsize size-top)
        splitted (split-at size-top sorted-all)
        top-percent-members (nth splitted 0)
        rest-offspring (take size-rest
                             (sort-by
                              :fitness ftype
                              (for [i (nth splitted 1)
                                    :when (if (= (:tag i) 1) true)] i)))]
    (assoc popu :poplist
      (concat top-percent-members rest-offspring))))


(defn- parent-selection
  "TODO DocString"
  [popu ftype percent popsize adapted?]
  (let [sorted (sort-by :fitness ftype (:poplist popu))
	size-parents (int (* percent popsize))
	parents (take size-parents sorted)

	ext-parents (if adapted?
		     (scramble (take (* 2 popsize) (cycle parents)))
		     (scramble (take popsize (cycle parents))))
	splitted (if adapted?
		   (split-at popsize ext-parents)
		   (split-at (/ popsize 2) ext-parents))]
    (map #(list %1 %2) (nth splitted 0) (nth splitted 1))))

(defn ga
  "Starts the GA algorithm.
  Parameter:
  fitness        - defines the fitness function used.
  The only argument is the position of the particle (a double-vector).
  ftype          - defines the type of optimization problem
                   (minimize (<) or maximize (>)).
  dim            - number of dimensions in solution.
  popsize        - size of the population.
  par-perc       - percentage of parents chosen for parent selection.
  surv-perc      - percentage of population chosen for next genertion.
  adpated?       - use mutation step size adaption (boolean).
  max-iterations - maximum number of iterations.
  "
  [fitness ftype dim popsize par-perc surv-perc adapted? max-iterations]
  (let [popu (init-population popsize 22 8)
        popu-evaluated (evaluate-all-firstrun popu fitness)]
    (loop [runs max-iterations popu popu-evaluated]
      (if DEBUG?
        (println "best:"
                 (:fitness (first (sort-by :fitness ftype (:poplist popu))))))
      (if (zero? runs)
        (first (sort-by :fitness ftype (:poplist popu)))
        (let [parents (parent-selection popu ftype par-perc popsize adapted?)
              popu-with-children (generate-offspring popu parents adapted?)
              popu-eva (evaluate-all popu-with-children fitness adapted?)
              popu-surv (survivor-selection popu-eva ftype surv-perc popsize)
              new-pop (assoc popu-surv :poplist
                        (for [i (:poplist popu-surv)] (assoc i :tag 0)))]
          (if DEBUG?
            (do (println "no. parent pairs:" (count parents))
              (println "no. parents + children:"
                       (count (:poplist popu-with-children)))
              (println "no. survivors: " (count (:poplist popu-surv)))))
          (recur (dec runs) new-pop))))))

(def DEBUG? false)

;; usage
;;(seq (chromo-to-phenotype (:chromosome
;;                           (ga griewank < 4 30 0.4 0.4 false 50))))


;;
;; Comments
;;

(comment "

  Problem:
  -------

  In genetic algorithms, the term of premature convergence means that a
  population for an optimization problem converged too early, resulting in
  being suboptimal. In this context, the parental solutions, through the aid of
  genetic operators, are not able to generate offsprings that are superior to
  their parents. Premature convergence can happen in case of loss of genetic
  variation (every individual in the population is identical, see convergence).

  Strategies for preventing premature convergence and to regain genetic
  variation can be:

  - a mating strategy called incest prevention,
  - uniform crossover,
  - favored replacement of similar individuals (preselection or crowding),
  - segmentation of individuals of similar fitness (fitness sharing),
  - increasing population size.

  The genetic variation can also be regained by mutation though this process is
  highly random!!

  Solution:
  --------

  Fitness sharing
  ")