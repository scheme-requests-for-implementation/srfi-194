; SPDX-FileCopyrightText: 2020 Arvydas Silanskas
; SPDX-FileCopyrightText: 2020 Linas Vepštas
; SPDX-License-Identifier: MIT
;
; zipf-test.scm
; Unit tests for the Zipf (zeta) distribution.
;
; Created by Linas Vepstas 10 July 2020
; Part of srfi-194

; ------------------------------------------------------------------
; Debug utility for gnuplot graphing.
; You can use this to dump a vector to a tab-delimited file.
(define (vector-to-file vec filename)
  (define (write-vec)
    (for-each
      (lambda (i)
        (define index (+ i 1))
        (define val (vector-ref vec i))
        (display index)
        (display "  ")
        (display val)
        (newline))
      (iota (vector-length vec))))
  (with-output-to-file filename write-vec))

; ------------------------------------------------------------------
; Test harness for exploring the Zipf (Riemann/Hurwicz zeta) distribution
; parameter space.
;
;     (test-zipf TEST-ID NVOCAB ESS QUE REPS TOL)
;
; * TEST-ID -- String ID, for debugging.
; * The next three parameters are presented to the generator as
;     (make-zipf-generator NVOCAB ESS QUE)
;   ++  NVOCAB -- Size of vocabulary to select from.
;   ++  ESS -- The Riemann zeta "s" exponent.
;   ++  QUE -- The Hurwicz zeta "q" offset.
; * REPS -- The number of samples to draw from the distribution.
; * TOL -- The test tolerance, governing the expected failure rate.
;
; The algorithm is roughly:
;   Take REPS samples (make-zipf-generator NVOCAB ESS QUE)
;   Accumulate them into NVOCAB histogram bins.
;   Normalize counts to unit probability (i.e. divide by NVOCAB)
;
; The resulting distribution should uniformly converge to C/(k+q)^s
; for 1 <= k <= NVOCAB where C is a normalization constant.
;
; This compares the actual distribution to the expected convergent
; and reports an error if it is not within TOL of the convergent.
; i.e. it computes the Banach l_0 norm of (distribution-convergent)
; TOL is to be given in units of standard deviations. So, for example,
; setting TOL to 6 gives a six-sigma bandpass, allowing the tests to
; usually pass.
;
(define (test-zipf TEST-ID NVOCAB ESS QUE REPS TOL)

  ; Default random number generator
  (define ZGEN make-zipf-generator)

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
  (define frequency (vector-map (lambda (i n) (/ n REPS)) bin-counts))
  (define probility (vector-map (lambda (i n) (inexact n)) frequency))

  ; Sequence 1..NVOCAB
  (define seq
    (cond-expand
     (gambit
      (list->vector (iota NVOCAB 1)))
     (else
      (vector-unfold (lambda (i x) (values x (+ x 1))) NVOCAB 1))))

  ; Sequence  1/(k+QUE)^ESS
  (define inv-pow
    (vector-map (lambda (i k) (expt (+ k QUE) (- (inexact ESS)))) seq))

  ; Hurwicz harmonic number sum_1..NVOCAB 1/(k+QUE)^ESS
  (define hnorm
    (vector-fold
      (lambda (i sum cnt) (+ sum cnt)) 0 inv-pow))

  ; The expected distribution
  (define expect
    (vector-map (lambda (i x) (/ x hnorm)) inv-pow))

  ; Convert to floating point.
  (define prexpect (vector-map (lambda (i x) (inexact x)) expect))

  ; The difference
  (define diff (vector-map (lambda (i x y) (- x y)) probility prexpect))

  ; Re-weight the tail by k^{s/2}. This seems give a normal error
  ; distribution. ... at least, for small q. Problems for large q
  ; and with undersampling; so we hack around that.
  (define err-dist
    (if (< 10 QUE) diff
        (vector-map (lambda (j i x) (* x (expt (+ i 1) (* 0.5 ESS))))
                    (list->vector (iota (vector-length diff)))
                    diff)))

  ; Normalize to unit root-mean-square.
  (define rms (/ 1 (sqrt (* 2 3.141592653 REPS))))
  (define norm-dist (vector-map (lambda (i x) (/ x rms)) err-dist))

  ; Maximum deviation from expected distribution (l_0 norm)
  (define l0-norm
    (vector-fold
      (lambda (i sum x) (if (< sum (abs x)) (abs x) sum)) 0 norm-dist))

  ; The mean.
  (define mean
    (/ (vector-fold (lambda (i sum x) (+ sum x)) 0 norm-dist)
       NVOCAB))

  (define root-mean-square
    (sqrt (/ (vector-fold (lambda (i sum x) (+ sum (* x x))) 0 norm-dist)
             NVOCAB)))

  ; The total counts in the bins should be equal to REPS
  (test-assert TEST-ID
    (equal? REPS
            (vector-fold
              (lambda (i sum cnt) (+ sum cnt)) 0 bin-counts)))

  ; Test for uniform convergence.
  (test-assert TEST-ID (<= l0-norm TOL))

  ; Should not random walk too far away.
  ; Could tighten this with a proper theory of the error distribution.
  (test-assert TEST-ID (< (abs mean) 3))
  ; I don't understand the error distribution ....
  ; (test-assert (and (< 0.4 root-mean-square) (< root-mean-square 1.5)))

  ; Utility debug printing
  ;(vector-to-file probility "probility.dat")
  ;(vector-to-file prexpect "prexpect.dat")
  ;(vector-to-file diff "diff.dat")
  #f)

; Explore the parameter space.
(define (zipf-test-group)
  ; (test-begin "srfi-194-zipf")

  ; The unit test computes something that is "almost" a standard
  ; deviation for the error distribution. Except, maybe not quite,
  ; I don't fully understand the theory. So most tests seem to come
  ; in fine in well-under a six-sigma deviation, but some of the wilder
  ; parameter choices misbehave, so six-sigma doesn't always work.
  ; Also, when the number of bins is large, its easy to under-sample;
  ; some bins end up empty and the std-dev is thrown off as a result.
  ; Thus, the tolerance bounds below are hand-adjusted.
  (define six-sigma 6.0)

  (define hack-que 3.0)

  ; Zoom into s->1
  (test-zipf "zoom-1"   30 1.1     0 1000 six-sigma)
  (test-zipf "zoom-2"   30 1.01    0 1000 six-sigma)
  (test-zipf "zoom-3"   30 1.001   0 1000 six-sigma)
  (test-zipf "zoom-4"   30 1.0001  0 1000 six-sigma)
  (test-zipf "zoom-5"   30 1.00001 0 1000 six-sigma)

  (test-zipf "zoom-6"   30 (+ 1 1e-6)  0 1000 six-sigma)
  (test-zipf "zoom-8"   30 (+ 1 1e-8)  0 1000 six-sigma)
  (test-zipf "zoom-10"  30 (+ 1 1e-10) 0 1000 six-sigma)
  (test-zipf "zoom-12"  30 (+ 1 1e-12) 0 1000 six-sigma)
  (test-zipf "zoom-14"  30 (+ 1 1e-14) 0 1000 six-sigma)
  (test-zipf "zoom-inf" 30 1           0 1000 six-sigma)

  ; Verify improving uniform convergence
  (test-zipf "uniform-1" 30 1  0 10000   six-sigma)
  (test-zipf "uniform-2" 30 1  0 100000  six-sigma)

  ; Larger vocabulary
  (test-zipf "mid-voc-1" 300 1.1     0 10000 six-sigma)
  (test-zipf "mid-voc-2" 300 1.01    0 10000 six-sigma)
  (test-zipf "mid-voc-3" 300 1.001   0 10000 six-sigma)
  (test-zipf "mid-voc-4" 300 1.0001  0 10000 six-sigma)
  (test-zipf "mid-voc-5" 300 1.00001 0 10000 six-sigma)

  ; Larger vocabulary. Take more samples....
  (test-zipf "large-voc-1" 3701 1.1     0 40000 six-sigma)
  (test-zipf "large-voc-2" 3701 1.01    0 40000 six-sigma)
  (test-zipf "large-voc-3" 3701 1.001   0 40000 six-sigma)
  (test-zipf "large-voc-4" 3701 1.0001  0 40000 six-sigma)
  (test-zipf "large-voc-5" 3701 1.00001 0 40000 six-sigma)

  ; Huge vocabulary; few samples. Many bins will be empty,
  ; causing the std-dev to get large.
  (test-zipf "huge-voc-1" 43701 (+ 1 1e-6)  0 60000 9.5)
  (test-zipf "huge-voc-2" 43701 (+ 1 1e-7)  0 60000 9.5)
  (test-zipf "huge-voc-3" 43701 (+ 1 1e-9)  0 60000 9.5)
  (test-zipf "huge-voc-4" 43701 (+ 1 1e-12) 0 60000 9.5)
  (test-zipf "huge-voc-5" 43701 1           0 60000 9.5)

  ; Large s, small range
  (test-zipf "big-s-lo-1" 5 1.1     0 1000 six-sigma)
  (test-zipf "big-s-lo-2" 5 2.01    0 1000 six-sigma)
  (test-zipf "big-s-lo-3" 5 4.731   0 1000 six-sigma)
  (test-zipf "big-s-lo-4" 5 9.09001 0 1000 six-sigma)
  (test-zipf "big-s-lo-5" 5 13.45   0 1000 8.0)

  ; Large s, larger range. Most histogram bins will be empty
  ; so allow much larger error margins.
  (test-zipf "bis-mid-1" 130 1.5     0 30000 six-sigma)
  (test-zipf "bis-mid-2" 130 2.03    0 30000 9.0)
  (test-zipf "bis-mid-3" 130 4.5     0 30000 16.0)
  (test-zipf "bis-mid-4" 130 6.66    0 30000 24.0)

  ; Verify that accuracy improves with more samples.
  (test-zipf "samp-bi-1" 129 1.1     0 10000 six-sigma)
  (test-zipf "samp-bi-2" 129 1.01    0 10000 six-sigma)
  (test-zipf "samp-bi-3" 129 1.001   0 10000 six-sigma)
  (test-zipf "samp-bi-4" 129 1.0001  0 10000 six-sigma)
  (test-zipf "samp-bi-5" 129 1.00001 0 10000 six-sigma)

  ; Non-zero Hurwicz parameter
  (test-zipf "hurw-1" 131 1.1     0.3    10000 six-sigma)
  (test-zipf "hurw-2" 131 1.1     1.3    10000 six-sigma)
  (test-zipf "hurw-3" 131 1.1     6.3    10000 six-sigma)
  (test-zipf "hurw-4" 131 1.1     20.23  10000 six-sigma)

  ; Negative Hurwicz parameter. Must be greater than branch point at -0.5.
  (test-zipf "hneg-1" 81 1.1     -0.1   1000 six-sigma)
  (test-zipf "hneg-2" 81 1.1     -0.3   1000 six-sigma)
  (test-zipf "hneg-3" 81 1.1     -0.4   1000 six-sigma)
  (test-zipf "hneg-4" 81 1.1     -0.499 1000 six-sigma)

  ; A walk into a stranger corner of the parameter space.
  (test-zipf "big-h-1" 131 1.1     41.483 10000 hack-que)
  (test-zipf "big-h-2" 131 2.1     41.483 10000 hack-que)
  (test-zipf "big-h-3" 131 6.1     41.483 10000 hack-que)
  (test-zipf "big-h-4" 131 16.1    41.483 10000 hack-que)
  (test-zipf "big-h-5" 131 46.1    41.483 10000 hack-que)
  (test-zipf "big-h-6" 131 96.1    41.483 10000 hack-que)

  ; A still wilder corner of the parameter space.
  (test-zipf "huhu-1" 131 1.1     1841.4 10000 hack-que)
  (test-zipf "huhu-2" 131 1.1     1.75e6 10000 hack-que)
  (test-zipf "huhu-3" 131 2.1     1.75e6 10000 hack-que)
  (test-zipf "huhu-4" 131 12.1    1.75e6 10000 hack-que)
  (test-zipf "huhu-5" 131 42.1    1.75e6 10000 hack-que)

  ; Lets try s less than 1
  (test-zipf "small-s-1" 35 0.9     0 1000 six-sigma)
  (test-zipf "small-s-2" 35 0.99    0 1000 six-sigma)
  (test-zipf "small-s-3" 35 0.999   0 1000 six-sigma)
  (test-zipf "small-s-4" 35 0.9999  0 1000 six-sigma)
  (test-zipf "small-s-5" 35 0.99999 0 1000 six-sigma)

  ; Attempt to force an overflow
  (test-zipf "ovfl-1" 437 (- 1 1e-6)  0 1000 six-sigma)
  (test-zipf "ovfl-2" 437 (- 1 1e-7)  0 1000 six-sigma)
  (test-zipf "ovfl-3" 437 (- 1 1e-9)  0 1000 six-sigma)
  (test-zipf "ovfl-4" 437 (- 1 1e-12) 0 1000 six-sigma)

  ; Almost flat distribution
  (test-zipf "flat-1" 36 0.8     0 1000 six-sigma)
  (test-zipf "flat-2" 36 0.5     0 1000 six-sigma)
  (test-zipf "flat-3" 36 0.1     0 1000 six-sigma)

  ; A visit to crazy-town -- increasing, not decreasing exponent
  (test-zipf "neg-s-1" 36 0.0     0 1000 six-sigma)
  (test-zipf "neg-s-2" 36 -0.1    0 1000 six-sigma)
  (test-zipf "neg-s-3" 36 -1.0    0 1000 six-sigma)
  (test-zipf "neg-s-4" 36 -3.0    0 1000 six-sigma)

  ; More crazy with some Hurwicz on top.
  (test-zipf "neg-shu-1" 16 0.0     0.5 1000 six-sigma)
  (test-zipf "neg-shu-2" 16 -0.2    2.5 1000 six-sigma)
  (test-zipf "neg-shu-3" 16 -1.3    10  1000 six-sigma)
  (test-zipf "neg-shu-4" 16 -2.9    100 1000 six-sigma)

  ; (test-end "srfi-194-zipf")
  )
