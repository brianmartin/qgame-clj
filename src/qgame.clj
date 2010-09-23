;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; qgame.clj
;; 2010 Brian Martin (btm08@hampshire.edu)
;;
;; Translated from qgame.lisp (version 1.20031226):
;;   c) 1999-2004, Lee Spector (lspector@hampshire.edu)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(ns qgame
  (:require [clojure.contrib.math :as math]
            [clojure.set :as set])
  (:gen-class))

(use '(incanter core))
(use 'traceme)

(defn -main [& args]
  (println "hello"))

(defstruct quantum-system-struct
             :number-of-qubits     ;number of qubits in the system
             :amplitudes           ;an array of amplitudes
             :prior-probability    ;the prob. of having reached this
                                   ;system in the first place
             :oracle-count         ;number of oracle calls throughout
                                   ;the system's history
             :measurement-history  ;a list of measurements and their
                                   ;results in the history of this system
             :instruction-history
             :program
             :qubit-numbers
             :amplitude-address
             )

(defn quantum-system
  "Initializes and returns quantum-system-struct."
  [number-of-qubits pgm]
  (struct-map quantum-system-struct
    :number-of-qubits number-of-qubits
    :amplitudes (let [dim (math/expt 2 number-of-qubits)]
                  (vec (cons 1.0 (repeat (dec dim) 0.0))))
    :prior-probability 0
    :oracle-count 0
    :instruction-history []
    :program pgm
    :qubit-numbers (range number-of-qubits)
    :amplitude-address (vec (repeat number-of-qubits 0))
    ))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; quantum computer manipulation utilities

(defn set-address-components
  "returns qsys with :amplitude-address to refer to a particular amplitude, as
    indicated by the bits in count."
  [qsys cnt qubits]
  (with-local-vars [new-addr (qsys :amplitude-address)]
    (doseq [i (range (count qubits))]
      (var-set new-addr
        (assoc @new-addr (nth qubits i)
          (if (bit-test cnt i) 1 0))))
    (assoc qsys :amplitude-address @new-addr)))

(defn map-qubit-combinations
  "Calls function once for each of the 1/0 combinations of the provided
  qubits, with the right-most qubit varying the most
  NOTE: function must take and return a qsys"
  [qsys function qubits]
  (let [num-iterations (math/expt 2 (count qubits))
        qubits-reversed (reverse qubits)]
    (with-local-vars [qsys-mut qsys]
      (doseq [i (range num-iterations)]
        (var-set qsys-mut
                 (function
                   (set-address-components @qsys-mut i qubits-reversed))))
      @qsys-mut)))
      
(defn addressed-amplitude-index
  "Returns the index of the currently addressed amplitude."
  [qsys]
  (apply + (for [i (range (qsys :number-of-qubits))]
       (if (not (zero? (nth (qsys :amplitude-address) i)))
         (math/expt 2 i)
         0))))
  
(defn get-addressed-amplitude
  "Returns the amplitude currently addressed."
  [qsys]
  (nth (qsys :amplitudes) (addressed-amplitude-index qsys)))

(defn set-addressed-amplitude 
  "Sets the amplitude currently addressed to new-value"
  [qsys new-value]
  (assoc qsys :amplitudes 
    (assoc (qsys :amplitudes) (addressed-amplitude-index qsys) new-value)))

(defn extract-column
  "Returns a column from the amplitudes obtained by varying the listed qubits,
  with the right-most qubit varying the fastest."
  [qsys qubits-to-vary]
  (reverse
    (:col
      (map-qubit-combinations
        (assoc qsys :col [])
        (fn [qsys]
            (update-in qsys [:col] conj (get-addressed-amplitude qsys)))
        qubits-to-vary))))

(defn install-column 
  "Installs the given column in the amplitude positions obtained by varying the
  listed qubits, with the right-most qubit varying the fastest."
  [qsys column qubits-to-vary]
  (with-local-vars [c column]
    (map-qubit-combinations
      qsys
      (fn [q]
        (let [f (first @c)]
          (var-set c (rest @c))
          (set-addressed-amplitude q f)))
      qubits-to-vary)))

(defn apply-operator 
  "Applies the given matrix-form operator to the given qubits."
  [qsys operator qubits]
  (map-qubit-combinations
    qsys
    (fn [q] 
      (let [pre-column (extract-column q qubits)
            post-column (mmult operator pre-column)]
        (install-column q post-column qubits)))
    (seq (set/difference 
           (set (qsys :qubit-numbers))
           (set qubits)))))

(defn qc-output-probabilities
  "Returns a list of the probabilities for all combinations for the given
  qubits in binary order with the rightmost qubit varying fastest."
  [qsys qubits]
  (let [other-qubits (set/difference (set (qsys :qubit-numbers))
                                     (set qubits))]
    (with-local-vars [probabilities '()]
      (map-qubit-combinations
        qsys
        (fn [q]
          (var-set
            probabilities
            (with-local-vars [probability 0]
              (cons
                (do
                  (map-qubit-combinations
                    qsys
                    (fn [q]
                      (var-set probability (+ @probability (math/expt (abs (get-addressed-amplitude q)) 2)))
                      q)
                    other-qubits)
                  @probability)
                @probabilities)))
          q)
        qubits)
      (reverse @probabilities))))

(defn add-lists 
  "Helper for multi-qsys-output-probabilities.  Takes two lists and adds
  elements of the same index."
  [l1 l2]
  (if (empty? l1)
    nil
    (cons (+ (first l1) (first l2))
          (add-lists (rest l1) (rest l2)))))

(defn multi-qsys-output-probabilities
 "Returns a list of the probabilities for all combinations for the
  given qubits, in binary order with the rightmost qubit varying fastest.
  This function takes a LIST of quantum systems as input and sums the
  results across all systems." 
  [qsys-list qubits]
  (let [probabilities (map #(qc-output-probabilities % qubits) qsys-list)]
    (reduce add-lists probabilities)))


(defn expected-oracles
  "Returns the expected number of oracle calls for the given set of
  quantum systems."
  [qsys-list]
  (println qsys-list)
  (reduce + (map #(* (% :prior-probability) (% :oracle-count)) qsys-list)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; oracle gates

(defn binary-operator-matrix
  "Returns a matrix operator for a binary function with the
  given tt-right-column as the right column of the truth table."
  [tt-right-column]
  (let [column-length (count tt-right-column)
        operator-size (* 2 column-length)]
    (matrix
      (loop [i 0
             m (vec (repeat operator-size (vec (repeat operator-size 0))))]
      (let [offset (* i 2)
            new-m (if (zero? (nth tt-right-column i))
                    (-> m (assoc-in [offset offset] 1)
                          (assoc-in [(inc offset) (inc offset)] 1))
                    (-> m (assoc-in [offset (inc offset)] 1)
                          (assoc-in [(inc offset) offset] 1)))]
        (if (= i (dec column-length))
          new-m
          (recur (inc i) new-m)))))))

(defn oracle
  "Applies the oracle operator built from tt-right-column, which
  is the right column of the corresponding truth table."
  [qsys tt-right-column & qubits]
  (apply-operator
    (update-in qsys [:oracle-count] inc)
    (binary-operator-matrix tt-right-column)
    qubits))

(defn limited-oracle
  "If (oracle-count qsys) is less than max-calls then this applies 
  the oracle operator built from tt-right-column, which is the right 
  column of the corresponding truth table. Otherwise this does nothing."
  [qsys max-calls tt-right-column & qubits]
  (if (< (qsys :oracle-count) max-calls)
    (apply-operator
      (update-in qsys [:oracle-count] inc)
      (binary-operator-matrix tt-right-column)
      qubits)
    qsys))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; other quantum gates

(def s (/ 1 (sqrt 2)))

(defn qnot
  "Quantum NOT gate"
  [qsys q]
  (apply-operator qsys
                  (matrix [[0 1]
                           [1 0]])
                  (list q)))

(defn cnot
  "Quantum controlled NOT gate"
  [qsys q1 q2]
  (apply-operator qsys
                  (matrix [[1 0 0 0]
                           [0 1 0 0]
                           [0 0 0 1]
                           [0 0 1 0]])
                  (list q1 q2)))

(defn srn
  "Quantum Square-root-of-NOT gate"
  [qsys q]
  (apply-operator qsys
                  (matrix [[s (- s)]
                           [s    s ]])
                  (list q)))

(defn nand
  "Quantum NAND gate"
  [qsys q1 q2 q3]
  (apply-operator qsys
                  (binary-operator-matrix '(1 1 1 0))
                  (list q1 q2 q3)))

(defn hadamard
  "Quantum Hadamard gate"
  [qsys q]
  (apply-operator qsys
                  (matrix [[s    s ]
                           [s (- s)]])
                  (list q)))

(defn u-theta
  "Quantum U-theta (rotation) gate"
  [qsys q theta]
  (apply-operator qsys
                  (matrix [[   (cos theta)  (sin theta)]
                           [(- (sin theta)) (cos theta)]])
                  (list q)))

(defn cphase
  "Quantum conditional phase gate."
  [qsys q1 q2 alpha]
  (apply-operator qsys
                  (diag 1 1 1 (exp (* (sqrt -1.0) alpha)))
                  (list q1 q2)))

;TODO: how to deal with i?
(defn u2
  "Quantum U2 gate, implemented as:
       e^(i(-phi-psi+alpha))*cos(theta)  e^(i(-phi+psi+alpha))*sin(-theta)
       e^(i(phi-psi+alpha))*sin(theta)   e^(i(phi+psi+alpha))*cos(theta)"
  [qsys q phi theta psi alpha]
  )

(defn swap
  "A quantum gate that swaps the amplitudes for the two specified qubits."
  [qsys q1 q2]
  (apply-operator qsys
                  (matrix [[1 0 0 0]
                           [0 0 1 0]
                           [0 1 0 0]
                           [0 0 0 1]])
                  (list q1 q2)))

(defn printamps
  "For use in quantum programs; causes the amplitudes of the executing
     quantum system to be printed."
  [qsys]
  (println (qsys :amplitudes)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; utilities for measurement and branching

(defn end 
  "Marks the end of a measurement branch; this has no effect otherwise."
  [qsys] qsys)

(defn distance-to-next-unmatched-end
  "Returns 0 if there is no unmatched (end) in pgm; otherwise returns
  the number of instructions to the next unmatched (end) (counting the (end))."
  ([pgm] (distance-to-next-unmatched-end pgm 0 0 0))
  ([pgm num-measures num-ends distance-so-far]
   (if (empty? pgm)
     0
     (if (= (ffirst pgm) 'end)
       (if (zero? num-measures)
         (inc distance-so-far)
         (if (odd? num-ends)
           (distance-to-next-unmatched-end 
             (rest pgm)
             (dec num-measures)
             (dec num-ends)
             (inc distance-so-far))
           (distance-to-next-unmatched-end
             (rest pgm)
             num-measures
             (inc num-ends)
             (inc distance-so-far))))
       (if (= (ffirst pgm) 'measure)
           (distance-to-next-unmatched-end 
             (rest pgm)
             (dec num-measures)
             num-ends
             (inc distance-so-far))
           (distance-to-next-unmatched-end
             (rest pgm)
             num-measures
             num-ends
             (inc distance-so-far)))))))

(defn my-subseq [s x y] (drop x (drop-last (- (count s) y) s)))

(defn without-if-branch
  "Assuming that a MEASURE form has just been removed from the given
  program, returns the remainder of the program without the IF (measure-1)
  branch"
  [pgm]
  (let [distance-to-first-unmatched-end 
           (distance-to-next-unmatched-end pgm)
        distance-from-first-to-second-unmatched-end
           (distance-to-next-unmatched-end
               (drop distance-to-first-unmatched-end pgm))]
    (if (zero? distance-to-first-unmatched-end)
      ; it's all the if part
      nil
      ; there is some else part
      (if (zero? distance-from-first-to-second-unmatched-end)
        ; the else never ends
        (->> pgm (drop distance-to-first-unmatched-end))
        ; the else does end
        (list (my-subseq pgm distance-to-first-unmatched-end
                             (+ distance-to-first-unmatched-end
                                distance-from-first-to-second-unmatched-end
                                -1))
              (->> pgm (drop (+ distance-to-first-unmatched-end
                                distance-from-first-to-second-unmatched-end))))))))

(defn without-else-branch
  [pgm]
  (let [distance-to-first-unmatched-end 
           (distance-to-next-unmatched-end pgm)
        distance-from-first-to-second-unmatched-end
           (distance-to-next-unmatched-end
               (drop distance-to-first-unmatched-end pgm))]
    (if (zero? distance-to-first-unmatched-end)
      ; it's all the if part
      pgm
      ; there is some else part
      (if (zero? distance-from-first-to-second-unmatched-end)
        ; the else never ends
        (my-subseq pgm 0 (dec distance-to-first-unmatched-end))
        ; the else does end
        (list (my-subseq pgm 0 (dec distance-to-first-unmatched-end))
              (my-subseq pgm (+ distance-to-first-unmatched-end
                             distance-from-first-to-second-unmatched-end)))))))

(defn force-to
  "Collapses a quantum system to the provided measured-value for the
  provided qubit."
  [measured-value qubit qsys]
  (with-local-vars [qu (list qubit)]
    (map-qubit-combinations 
      qsys
      (fn [q]
        (let [pre-column (extract-column q qu)
              new-column (if (zero? measured-value)
                           (list (first pre-column) 0)
                           (list 0 (second pre-column)))]
          (install-column q new-column qu)))
      (remove #(= qu %) (qsys :qubit-numbers)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; top-level functions

;TODO: assoc measurement-history ???
(defn run-qsys
  "Takes a quantum system and returns the list of quantum systems that results
  from the execution of its program."
  [qsys]
  (if (or (nil? (qsys :program))
          (zero? (qsys :prior-probability)))
    (list qsys)
    (let [instruction (first (qsys :program))
          q (update-in qsys [:instruction-history] conj instruction)]
      (if (= (first instruction) 'halt)
        (list q)
        (if (= (first instruction) 'measure)
          ; it's a measurement so split state and return list of results
          (let [measurement-qubit (second instruction)
                probabilities (qc-output-probabilities q (list measurement-qubit))]
            (list
              ;; 1 branch
              (run-qsys
                (force-to 1 measurement-qubit
                  (-> q (assoc :prior-probability (second probabilities))
                        (update-in [:program] #(without-else-branch (rest %))))))
              ;; 0 branch
              (run-qsys
                (force-to 0 measurement-qubit
                  (-> q (assoc [:prior-probability] (first probabilities))
                        (update-in [:program] #(without-if-branch (rest %))))))))
          (let [resulting-sys (apply (first instruction) (cons qsys (rest instruction)))]
            (run-qsys (assoc resulting-sys :program (rest (resulting-sys :program))))))))))

(defn execute-quantum-program
  "Executes the provided quantum program with the specified number of qubits
  and the provided oracle truth table, returning a list of the resulting 
  quantum systems."
  ([pgm num-qubits] (execute-quantum-program pgm num-qubits nil))
  ([pgm num-qubits oracle-tt]
    (run-qsys (quantum-system num-qubits
                              (map #(replace {'ORACLE-TT oracle-tt} %) pgm)))))

(defn test-quantum-program
  "The top-level function to evaluate a quantum program relative to a list of
  a list of (oracle value) cases. Returns a list of: misses, max-error, average-error,
  max-expected-oracles, and average-expected-oracles."
  [pgm {:keys [num-qubits test-cases final-measurement-qubits threshold debug]
         :or {debug 2}}]
  (println pgm)
  (let [stats-all
        (let [num-cases (count test-cases)
              stats
                    (loop [cases test-cases
                           misses 0
                           max-error 0
                           total-error 0
                           max-expected-oracles 0
                           total-expected-oracles 0]
                      (let [resulting-systems (execute-quantum-program pgm num-qubits (ffirst cases))
                            raw-error (- 1.0 (nth (multi-qsys-output-probabilities resulting-systems final-measurement-qubits)
                                                  (second (first cases))))
                            expected-oracles (expected-oracles resulting-systems)]
                        (when (>= debug 2)
                          (println "case: " c "  error: " raw-error))
                        (if (empty? (rest cases))
                          {:misses misses :max-error max-error :total-error total-error :max-expected-oracles max-expected-oracles
                             :total-expected-oracles total-expected-oracles}
                          (recur
                            (rest cases)
                            (if (> raw-error threshold) (inc misses) misses)
                            (if (> raw-error max-error) raw-error max-error)
                            (+ total-error raw-error)
                            (if (> expected-oracles max-expected-oracles) expected-oracles max-expected-oracles)
                            (+ total-expected-oracles expected-oracles)))))]
            (-> stats (assoc :average-error  (/ (stats :total-error) num-cases))
                      (assoc :average-expected-oracles (/ (stats :total-expected-oracles) num-cases))))]
      (when (>= debug 1)
        (println stats-all))
    stats-all))

(defn test-herbs-grover
  []
  (test-quantum-program 
   '((hadamard 2)
     (hadamard 1)
     (u-theta 0 ,(/ pi 4))
     (oracle ORACLE-TT 2 1 0)
     (hadamard 2)
     (cnot 2 1)
     (hadamard 2)
     (u-theta 2 ,(/ pi 2))
     (u-theta 1 ,(/ pi 2))
     )
   {
     :num-qubits 3
     :test-cases '(((1 0 0 0) 0)
              ((0 1 0 0) 1)
              ((0 0 1 0) 2)
              ((0 0 0 1) 3))
     :final-measurement-qubits '(2 1)
     :threshold 0.48}))
