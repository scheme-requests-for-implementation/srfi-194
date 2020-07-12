;
; zipf-test.scm
; Unit tests for the Zipf (zeta) distribution.
;
; Created by Linas Vepstas 10 July 2020
; Nominated for inclusion in srfi-194

; ------------------------------------------------------------------
; Debug utility for gnuplot graphing.
; You can use this to dump a vector to a tab-delimited file.
(define (vector-to-file vec filename)
  (define (write-vec)
    (generator-for-each
      (lambda (i) 
        (define index (+ i 1))
        (define val (vector-ref vec i))
        (display index)
        (display "  ")
        (display val)
        (newline))
      (make-iota-generator (vector-length vec))))
  (with-output-to-file filename write-vec))

; ------------------------------------------------------------------
; Simple test harness for exploring paramter space.
;
; Take REPS samples from the zeta distribution (ZGEN NVOCAB ESS QUE)
; Accumulate them into NVOCAB bins.
; Normalize to unit probability (i.e. divide by NVOCAB)
;
; The resulting distribution should uniformly converge to C/(k+q)^s
; for 1 <= k <= NVOCAB where C is a normalization constant.
;
; This compares the actual distribution to the expected convergent
; and reports an error if it is not within TOL of the convergent.
; i.e. it computes the Banach l_0 norm of (distribution-convergent)
;
(define (test-zipf ZGEN NVOCAB ESS QUE REPS TOL)

  ; Bin-counter containing accumulated histogram.
  (define bin-counts
    (let ((bin-counts (make-vector NVOCAB 0)))
     ; Accumulate samples into the histogram.
     (generator-for-each
       (lambda (SAMP)
         (define offset (- SAMP 1))
         (vector-set! bin-counts offset (+ 1 (vector-ref bin-counts offset))))
       (gtake (ZGEN NVOCAB ESS QUE) REPS))
     bin-counts))

  ; Verify the distribution is within tolerance.
  ; This is written out long-hand for easier debuggability.

  ; Frequency is normalized to be 0.0 to 1.0
  (define frequency (vector-map (lambda (n) (/ n REPS)) bin-counts))
  (define probility (vector-map (lambda (n) (inexact n)) frequency))

  ; Sequence 1..NVOCAB
  (define seq
    (generator->vector (make-iota-generator NVOCAB 1)))

  ; Sequence  1/(k+QUE)^ESS
  (define inv-pow (vector-map
                    (lambda (k) (expt (+ k QUE) (- (inexact ESS)))) seq))

  ; Hurwicz harmonic number sum_1..NVOCAB 1/(k+QUE)^ESS
  (define hnorm
    (vector-fold
      (lambda (sum cnt) (+ sum cnt)) 0 inv-pow))

  ; The expected distribution
  (define expect
    (vector-map (lambda (x) (/ x hnorm)) inv-pow))

  ; Convert to floating point.
  (define prexpect (vector-map inexact expect))

  ; The difference
  (define diff (vector-map (lambda (x y) (- x y)) probility prexpect))

  ; Maximum deviation from expected distribution (l_0 norm)
  (define l0-norm
    (vector-fold
      (lambda (sum x) 
        (if (< sum (abs x)) 
            (abs x) 
            sum)) 
      0 diff))

  ; The total counts in the bins should be equal to REPS
  (test-equal
    REPS
    (vector-fold
      (lambda (sum cnt) (+ sum cnt)) 0 bin-counts))
  
  ; Test for uniform convergence.
  (test-assert 
    (<= l0-norm TOL))

  ; Utility debug printing
  (vector-to-file probility "probility.dat")
  (vector-to-file prexpect "prexpect.dat")
  (vector-to-file diff "diff.dat")
  )

; Explore the parameter space.
; The error bounds have been selected to be sort-of-ish tight, in that
; the whole combined set of tests below will usually pass, failing only
; once out of every dozen(?) or two(?) invocations.  The failures are
; random but infrequent, exactly how they should be!
;
(define (zipf-test-group)

  ; Zoom into s->1
  (test-zipf make-zipf-generator 30 1.1     0 1000 4e-2)
  (test-zipf make-zipf-generator 30 1.01    0 1000 4e-2)
  (test-zipf make-zipf-generator 30 1.001   0 1000 4e-2)
  (test-zipf make-zipf-generator 30 1.0001  0 1000 4e-2)
  (test-zipf make-zipf-generator 30 1.00001 0 1000 4e-2)

  (test-zipf make-zipf-generator 30 (+ 1 1e-6)  0 1000 4e-2)
  (test-zipf make-zipf-generator 30 (+ 1 1e-8)  0 1000 4e-2)
  (test-zipf make-zipf-generator 30 (+ 1 1e-10) 0 1000 4e-2)
  (test-zipf make-zipf-generator 30 (+ 1 1e-12) 0 1000 4e-2)
  (test-zipf make-zipf-generator 30 (+ 1 1e-14) 0 1000 4e-2)
  (test-zipf make-zipf-generator 30 1           0 1000 4e-2)

  ; Verify improving uniform convergence
  (test-zipf make-zipf-generator 30 1  0 10000   9.5e-3)
  (test-zipf make-zipf-generator 30 1  0 100000  2.5e-3)
  (test-zipf make-zipf-generator 30 1  0 1000000 9.5e-4)

  ; Larger vocabulary
  (test-zipf make-zipf-generator 300 1.1     0 1000 4e-2)
  (test-zipf make-zipf-generator 300 1.01    0 1000 4e-2)
  (test-zipf make-zipf-generator 300 1.001   0 1000 4e-2)
  (test-zipf make-zipf-generator 300 1.0001  0 1000 4e-2)
  (test-zipf make-zipf-generator 300 1.00001 0 1000 4e-2)

  ; Larger vocabulary
  (test-zipf make-zipf-generator 3701 1.1     0 1000 4e-2)
  (test-zipf make-zipf-generator 3701 1.01    0 1000 4e-2)
  (test-zipf make-zipf-generator 3701 1.001   0 1000 4e-2)
  (test-zipf make-zipf-generator 3701 1.0001  0 1000 4e-2)
  (test-zipf make-zipf-generator 3701 1.00001 0 1000 4e-2)

  ; Huge vocabulary
  (test-zipf make-zipf-generator 43701 (+ 1 1e-6)  0 1000 2.5e-2)
  (test-zipf make-zipf-generator 43701 (+ 1 1e-7)  0 1000 2.5e-2)
  (test-zipf make-zipf-generator 43701 (+ 1 1e-9)  0 1000 2.5e-2)
  (test-zipf make-zipf-generator 43701 (+ 1 1e-12) 0 1000 2.5e-2)
  (test-zipf make-zipf-generator 43701 1           0 1000 2.5e-2)

  ; Large s, small range
  (test-zipf make-zipf-generator 5 1.1     0 1000 5e-2)
  (test-zipf make-zipf-generator 5 2.01    0 1000 5e-2)
  (test-zipf make-zipf-generator 5 4.731   0 1000 5e-2)
  (test-zipf make-zipf-generator 5 9.09001 0 1000 5e-2)
  (test-zipf make-zipf-generator 5 13.45   0 1000 5e-2)

  ; Large s, larger range
  (test-zipf make-zipf-generator 130 1.5     0 1000 4e-2)
  (test-zipf make-zipf-generator 130 2.03    0 1000 4e-2)
  (test-zipf make-zipf-generator 130 4.5     0 1000 4e-2)
  (test-zipf make-zipf-generator 130 6.66    0 1000 4e-2)

  ; Verify that accuracy improves with more samples.
  (test-zipf make-zipf-generator 129 1.1     0 10000 9.5e-3)
  (test-zipf make-zipf-generator 129 1.01    0 10000 9.5e-3)
  (test-zipf make-zipf-generator 129 1.001   0 10000 9.5e-3)
  (test-zipf make-zipf-generator 129 1.0001  0 10000 9.5e-3)
  (test-zipf make-zipf-generator 129 1.00001 0 10000 9.5e-3)

  ; Non-zero Hurwicz parameter
  (test-zipf make-zipf-generator 131 1.1     0.3    10000 9.5e-3)
  (test-zipf make-zipf-generator 131 1.1     1.3    10000 9.5e-3)
  (test-zipf make-zipf-generator 131 1.1     6.3    10000 9.5e-3)
  (test-zipf make-zipf-generator 131 1.1     20.23  10000 9.5e-3)

  ; Negative Hurwicz parameter. Must be greater than branch point at -0.5.
  (test-zipf make-zipf-generator 81 1.1     -0.1   1000 4e-2)
  (test-zipf make-zipf-generator 81 1.1     -0.3   1000 4e-2)
  (test-zipf make-zipf-generator 81 1.1     -0.4   1000 4e-2)
  (test-zipf make-zipf-generator 81 1.1     -0.499 1000 4e-2)

  ; A walk into a stranger corner of the parameter space.
  (test-zipf make-zipf-generator 131 1.1     41.483 10000 10.5e-3)
  (test-zipf make-zipf-generator 131 2.1     41.483 10000 10.5e-3)
  (test-zipf make-zipf-generator 131 6.1     41.483 10000 10.5e-3)
  (test-zipf make-zipf-generator 131 16.1    41.483 10000 10.5e-3)
  (test-zipf make-zipf-generator 131 46.1    41.483 10000 10.5e-3)
  (test-zipf make-zipf-generator 131 96.1    41.483 10000 10.5e-3)

  ; A still wilder corner of the parameter space.
  (test-zipf make-zipf-generator 131 1.1     1841.4 10000 9.9e-3)
  (test-zipf make-zipf-generator 131 1.1     1.75e6 10000 9.9e-3)
  (test-zipf make-zipf-generator 131 2.1     1.75e6 10000 9.9e-3)
  (test-zipf make-zipf-generator 131 12.1    1.75e6 10000 9.9e-3)
  (test-zipf make-zipf-generator 131 42.1    1.75e6 10000 9.9e-3)

  ; Lets try s less than 1
  (test-zipf make-zipf-generator 35 0.9     0 1000 4e-2)
  (test-zipf make-zipf-generator 35 0.99    0 1000 4e-2)
  (test-zipf make-zipf-generator 35 0.999   0 1000 4e-2)
  (test-zipf make-zipf-generator 35 0.9999  0 1000 4e-2)
  (test-zipf make-zipf-generator 35 0.99999 0 1000 4e-2)

  ; Attempt to force an overflow
  (test-zipf make-zipf-generator 437 (- 1 1e-6)  0 1000 3e-2)
  (test-zipf make-zipf-generator 437 (- 1 1e-7)  0 1000 3e-2)
  (test-zipf make-zipf-generator 437 (- 1 1e-9)  0 1000 3e-2)
  (test-zipf make-zipf-generator 437 (- 1 1e-12) 0 1000 3e-2)

  ; Almost flat distribution
  (test-zipf make-zipf-generator 36 0.8     0 1000 4e-2)
  (test-zipf make-zipf-generator 36 0.5     0 1000 4e-2)
  (test-zipf make-zipf-generator 36 0.1     0 1000 4e-2)

  ; A visit to crazy-town -- increasing, not decreasing exponent
  (test-zipf make-zipf-generator 36 0.0     0 1000 4e-2)
  (test-zipf make-zipf-generator 36 -0.1    0 1000 4e-2)
  (test-zipf make-zipf-generator 36 -1.0    0 1000 4e-2)
  (test-zipf make-zipf-generator 36 -3.0    0 1000 4e-2)

  ; More crazy with some Hurwicz on top.
  (test-zipf make-zipf-generator 16 0.0     0.5 1000 4e-2)
  (test-zipf make-zipf-generator 16 -0.2    2.5 1000 4e-2)
  (test-zipf make-zipf-generator 16 -1.3    10  1000 4e-2)
  (test-zipf make-zipf-generator 16 -2.9    100 1000 4e-2))
