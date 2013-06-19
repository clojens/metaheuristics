(ns metaheuristics.testfunctions
  "A multimodal test function. Mimimization problem. Widely used in
  global opimization problems. The input is one single double array with
  arbitrary length. Global minimum at x = (0,...0), f(x) = 0")

(defn rastrigin
  "In mathematical optimization, the Rastrigin function is a non-convex
  function used as a performance test problem for optimization algorithms.
  It is a typical example of non-linear multimodal function."
  [position]
  (let [sumvec (amap position i ret
                     (+ (- (Math/pow (aget position i) 2)
                           (* 10
                              (Math/cos (* 2
                                           (Math/PI)
                                           (aget position i)))))
                        10))]

    (areduce sumvec i ret 0
             (+ ret (aget sumvec i)))))

(defn griewank
  "A multimodal test function. Mimimization problem. Widely used in
  global opimization problems. The input is one single double array with
  arbitrary length. Global minimum at x = (0,...0), f(x) = 0"
  [position]
  (let [prod (areduce position i ret 0
                      (* ret (Math/cos (/ (aget position i)
                                          (Math/sqrt (+ i 1))))))
        sum (/ (areduce position i ret 0
                        (+ ret (Math/pow (aget position i) 2))) 4000)]
    (double
     (-
      (+ sum 1)
      prod))))
