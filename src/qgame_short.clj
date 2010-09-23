(ns qgame
  (:require [clojure.contrib.math :as math]
            [clojure.set :as set])
  (:gen-class))

(use '(incanter core))

(defn -main [& args])

(defstruct quantum-system-struct
             :number-of-qubits     ;number of qubits in the system
             :amplitudes           ;an array of amplitudes
             :qubit-numbers
             :amplitude-address
             )

(defn quantum-system
  [number-of-qubits]
  (struct-map quantum-system-struct
    :number-of-qubits number-of-qubits
    :amplitudes (let [dim (math/expt 2 number-of-qubits)]
                  (vec (cons 1.0 (repeat (dec dim) 0.0))))
    :qubit-numbers (range number-of-qubits)
    :amplitude-address (vec (repeat number-of-qubits 0))
    ))

(defn set-address-components
  [qsys cnt qubits]
  (let [old-addr (qsys :amplitude-address)]
    (assoc qsys :amplitude-address 
     (eval
       (conj
         (for [i qubits]
           `(assoc ~i
              (if (bit-test ~cnt ~i) 1 0)))
         old-addr
         '->)))))

(defn map-qubit-combinations
  [qsys function qubits]
  (let [num-iterations (math/expt 2 (count qubits))
        qubits-reversed (reverse qubits)]
    (eval
      (conj 
        (interleave 
          (for [i (range num-iterations)] `(set-address-components ~i '~qubits-reversed))
          (repeat num-iterations function))
        'qsys 
        '->))))
    
(defn apply-qnot-operator 
  [qsys operator qubits]
  (map-qubit-combinations
    qsys
    (fn [q] 
      (let [pre-column [1 0 0 0]
            post-column (mmult [[0 1] [1 0]]  pre-column)]
        (qnot q post-column qubits)))
    (seq (set/difference 
           (set (qsys :qubit-numbers))
           (set qubits)))))
