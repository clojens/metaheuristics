(ns metaheuristics.iwo
  "Local search can be used on problems that can be formulated as finding a
  solution maximizing a criterion among a number of candidate solutions.
  Local search algorithms move from solution to solution in the space of
  candidate solutions (the search space) by applying local changes, until a
  solution deemed optimal is found or a time bound is elapsed."
  (:use [clojure.set :only (intersection difference)]
        [clojure.contrib.math]
        [metaheuristics.testfunctions]))



;; Invasive weed optimization
;; ==========================

;;
;; Structural
;;

;; FIXME Refactor into defrecord(2)

(defstruct plant
  :seedlist :position :pfit :tfit)

(defstruct population
  :plantlist :gbest   :gworst :minseed   :maxseed :sigmaInit :sigmaFinal)


;;
;; Initialization
;;

(defn- init-plant
  "Initialize fauna."
  [nDimensions max]
  (agent (struct plant
		 (list)
		 (double-array (for [i (range nDimensions)] (rand max)))
		 (double 0.0)
		 (double 0.0))))

(defn- init-population
  "Initialize population using seeds."
  [nPlants nDimensions nSeedMin nSeedMax sigmaInit sigmaFinal max]
  (let [plants (map (fn [_] (init-plant nDimensions max)) (range nPlants))
        gbest  (init-plant nDimensions max)
        gworst (init-plant nDimensions max)]
    (agent (struct population plants gbest gworst
                   (int nSeedMin)
                   (int nSeedMax)
                   sigmaInit
                   sigmaFinal))))

;;
;; Scenarios
;;

(defn- set-bestworst
  "Set a best-worst case scenario."
  [population ftype]
  (let [plants (for [p (:plantlist population)] @p)
        sorted (sort-by :pfit ftype plants)]
    (assoc population
      :gbest (agent (first sorted))
      :gworst (agent (last sorted)))))


;;
;; Evaluation
;;

(defn- eval-seed
  "TODO DocString"
  [seed fitness]
  (let [#^doubles pos (:position seed)
        fit (fitness pos)]
    (assoc seed :pfit fit)))

(defn- eval-plantseeds
  "TODO DocString"
  [plant fitness]
  (let [seeds (:seedlist @plant)]
    (doseq [s seeds] (send s eval-seed fitness))
    (apply await seeds)))

;;
;; Creation and (re)generation
;;

(defn- create-new-seed
  "Todo ds"
  [pos sigma]
  (let [dim (count pos)
        pos-offset (double-array
                    (map (fn [_] (normal 0 sigma)) (range dim)))
        newpos (amap pos i ret
                     (+ (aget pos i) (aget pos-offset i)))]
    (agent (struct plant (list) newpos 0.0))))

(defn- generate-seeds
  "Generates new seeds"
  [plant population maxIt modulation iteration]
  (let [pfit          (double (:pfit plant))
        #^doubles pos (:position plant)
        gbestfit      (double (:pfit @(:gbest population)))
        gworstfit     (double (:pfit @(:gworst population)))
        minSeed       (int (:minseed population))
        maxSeed       (int (:maxseed population))
        sigmaInit     (double (:sigmaInit population))
        sigmaFinal    (double (:sigmaFinal population))
        nSeeds (int (+
                     (* (/ (- pfit gworstfit) (- gbestfit gworstfit)) maxSeed)
                     (* (/ (- pfit gbestfit) (- gworstfit gbestfit)) minSeed)))
        sigma (+
               (* (/ (Math/pow (- maxIt iteration) modulation)
                     (Math/pow maxIt modulation))
                  (- sigmaInit sigmaFinal))
               sigmaFinal)
        newSeeds (map (fn [_] (create-new-seed pos sigma)) (range nSeeds))]
    (assoc plant :seedlist newSeeds)))

(defn- competition
  "Plant life competes over valuable real-estate. This would be for example
  for both ground (best minerals, soil) as for sunlight (grow higher, expand
  on top level to gain more light = photosynthesis) and underground roots for
  water."
  [population ftype limit style]
  (let [seeds     (for [p (:plantlist population) s (:seedlist @p)] @s)
        plants    (for [p (:plantlist population)] @p)
        ;; best up to limit
        newpop    (take limit (sort-by :pfit ftype (concat seeds plants)))
        newagents (for [np newpop] (agent np))]
    (assoc population
      :plantlist newagents
      :gbest (agent (first newpop))
      :gworst (agent (last newpop)))))


;; FIXME Replace debug with macro'd dbg and monadic spears/arrows/diamonds.

(defn- grow
  "Grow plant life."
  [fitness ftype population maxIt nPlantsMax modulation iteration]
  ;; generate new seeds for each plant-agent
  (if DEBUG? (println "creating new seeds..."))
  (dorun (map #(send % generate-seeds @population maxIt modulation iteration)
              (:plantlist @population)))
  (apply await (:plantlist @population))

  (if DEBUG?
    (println "no. of seeds:"
             (count
              (for [p (:plantlist @population)
                    seeds (:seedlist @p)]
                seeds))))

  ;; update seed fitness
  (if DEBUG? (println "updating seed fitness..."))
  (dorun (pmap #(eval-plantseeds % fitness) (:plantlist @population)))

  ;; competitive exclusion
  (if DEBUG? (println "competition..."))
  (send population competition ftype nPlantsMax 1)
  (await population)
  (if DEBUG? (println "best in generation" iteration":"
                      (:pfit @(:gbest @population)) "\n---")))

(defn iwo
  "Starts the IWO algorithm.
  Algorithm based on http://dx.doi.org/10.1016/j.ecoinf.2006.07.003

  Parameters:
  ----------
  fitness        - the fitness function to be used. only one parameter:
                   the position of the particle (double-array).
  ftype          - defines the type of optimization problem
                   (minimize (<) or maximize (>)).
  dim            - number of dimensions in solution.
  nplants        - number of initial plants.
  nplants-max    - maximum number of plants in population.
  seed-min       - minimum number of seeds generated per plant.
  seed-max       - maximum number of seeds generated per plant.
  sigma-init     - inital value for sigma (standard deviation).
  sigma-final    - final value for sigma (standard deviation).
  modulation     -  modulation index (usually 3).
  max-iterations - maximum number of iterations.
  max-feat       - maximum value for one feature."
  [fitness ftype dim nplants nplants-max seed-min seed-max
   sigma-init sigma-final modulation max-iterations max-feat]
  (let [population (init-population nplants dim seed-min seed-max
                                    sigma-init sigma-final max-feat)]
    ;; one time eval of initial plants
    (dorun (map #(send % eval-seed fitness) (:plantlist @population)))
    (apply await (:plantlist @population))
    (send population set-bestworst ftype)
    (await population)
    ;; start IWO
    (dorun
     (map
      #(grow fitness ftype population max-iterations nplants-max modulation %)
		(range max-iterations)))

    ;; return best solution
    @(:gbest @population)))

(def DEBUG? false)

;; testing
;; (seq (:position (iwo griewank < 100 10 30 0 5 10 0.1 3 100 10.0)))