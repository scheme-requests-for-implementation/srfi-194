
; Guile-specific test infrastructure
(use-modules (srfi srfi-64))

; Debug utility for gnuplot
(define (vector-to-file vec filename)
	(let ((outport (open-file filename "w")))
   	(vector-for-each
			(lambda (i x) (format outport "~A	~A\n" (+ i 1) x))
			vec)
      (close outport)))

; Take REPS samples from the zeta distribution (ZGEN NVOCAB ESS QUE)
;
; and whaaaatever
;
(define (test-zipf ZGEN NVOCAB ESS QUE REPS TOL)

	; Bin-counter containing accumulated histogram.
	(define bin-counts (make-vector NVOCAB 0))

	; Accumulate samples into the histogram.
	(generator-for-each
		(lambda (SAMP)
			(define offset (- SAMP 1))
			(vector-set! bin-counts offset (+ 1 (vector-ref bin-counts offset))))
		(gtake (ZGEN NVOCAB ESS QUE) REPS))

	; The total counts in the bins should be equal to REPS
	(test-assert
		(equal? REPS
			(vector-fold
				(lambda (i sum cnt) (+ sum cnt)) 0 bin-counts)))

	; Verify the distribution is within tolerance.
	; This is written out long-hand for easier debuggability.

	; Frequency is normalized to be 0.0 to 1.0
	(define frequency (vector-map (lambda (i n) (/ n REPS)) bin-counts))
	(define probility (vector-map (lambda (i n) (exact->inexact n)) frequency))

	; Sequence 1..NVOCAB
	(define seq
		(vector-unfold (lambda (i x) (values x (+ x 1))) NVOCAB 1))

	; Sequence  1/(k+QUE)^ESS
	(define inv-pow (vector-map (lambda (i k) (expt (+ k QUE) (- ESS))) seq))

	; Hurwicz harmonic number sum_1..NVOCAB 1/(k+QUE)^ESS
	(define hnorm
		(vector-fold
			(lambda (i sum cnt) (+ sum cnt)) 0 inv-pow))
	; (format #t "Norm = ~A\n" (exact->inexact hnorm))

	; The expected distribution
	(define expect
		(vector-map (lambda (i x) (/ x hnorm)) inv-pow))

	; Convert to floating point.
	(define prexpect (vector-map (lambda (i x) (exact->inexact x)) expect))

	; The difference
	(define diff (vector-map (lambda (i x y) (- x y)) probility prexpect))

	; Maximum deviation from expected distribution (l_0 norm)
	(define l0-norm
		(vector-fold
			(lambda (i sum x) (if (< sum (abs x)) (abs x) sum)) 0 diff))

	; Test for uniform convergence.
	(test-assert (<= l0-norm TOL))

	(format #t "N=~D s=~9,6f q=~9,6f SAMP=~D norm=~9,6f tol=~9,6f ~A\n"
		NVOCAB ESS QUE REPS l0-norm TOL (if (<= l0-norm TOL) "PASS" "FAIL"))

	; Utility debug printing
	;(vector-to-file probility "foo.dat")
	;(vector-to-file prexpect "bar.dat")
	;(vector-to-file diff "baz.dat")

	#f
)

(test-begin "srfi-194-zipf")
(test-zipf make-zipf-generator/zri 30 1.1     0 1000 4e-2)
(test-zipf make-zipf-generator/zri 30 1.01    0 1000 4e-2)
(test-zipf make-zipf-generator/zri 30 1.001   0 1000 4e-2)
(test-zipf make-zipf-generator/zri 30 1.0001  0 1000 4e-2)
(test-zipf make-zipf-generator/zri 30 1.00001 0 1000 4e-2)

(test-zipf make-zipf-generator/zri 300 1.1     0 1000 4e-2)
(test-zipf make-zipf-generator/zri 300 1.01    0 1000 4e-2)
(test-zipf make-zipf-generator/zri 300 1.001   0 1000 4e-2)
(test-zipf make-zipf-generator/zri 300 1.0001  0 1000 4e-2)
(test-zipf make-zipf-generator/zri 300 1.00001 0 1000 4e-2)

(test-zipf make-zipf-generator/zri 3701 1.1     0 1000 4e-2)
(test-zipf make-zipf-generator/zri 3701 1.01    0 1000 4e-2)
(test-zipf make-zipf-generator/zri 3701 1.001   0 1000 4e-2)
(test-zipf make-zipf-generator/zri 3701 1.0001  0 1000 4e-2)
(test-zipf make-zipf-generator/zri 3701 1.00001 0 1000 4e-2)

(test-zipf make-zipf-generator/zri 43701 (+ 1 1e-6)  0 1000 2e-2)
(test-zipf make-zipf-generator/zri 43701 (+ 1 1e-7)  0 1000 2e-2)
(test-zipf make-zipf-generator/zri 43701 (+ 1 1e-9)  0 1000 2e-2)
(test-zipf make-zipf-generator/zri 43701 (+ 1 1e-12) 0 1000 2e-2)

(test-zipf make-zipf-generator/zri 5 1.1     0 1000 4e-2)
(test-zipf make-zipf-generator/zri 5 2.01    0 1000 4e-2)
(test-zipf make-zipf-generator/zri 5 4.731   0 1000 4e-2)
(test-zipf make-zipf-generator/zri 5 9.09001 0 1000 4e-2)
(test-zipf make-zipf-generator/zri 5 13.45   0 1000 4e-2)

(test-zipf make-zipf-generator/zri 130 1.5     0 1000 4e-2)
(test-zipf make-zipf-generator/zri 130 2.03    0 1000 4e-2)
(test-zipf make-zipf-generator/zri 130 4.5     0 1000 4e-2)
(test-zipf make-zipf-generator/zri 130 6.66    0 1000 4e-2)

(test-zipf make-zipf-generator/zri 129 1.1     0 10000 9.5e-3)
(test-zipf make-zipf-generator/zri 129 1.01    0 10000 9.5e-3)
(test-zipf make-zipf-generator/zri 129 1.001   0 10000 9.5e-3)
(test-zipf make-zipf-generator/zri 129 1.0001  0 10000 9.5e-3)
(test-zipf make-zipf-generator/zri 129 1.00001 0 10000 9.5e-3)

(test-zipf make-zipf-generator/zri 131 1.1     0.3    10000 9.5e-3)
(test-zipf make-zipf-generator/zri 131 1.1     1.3    10000 9.5e-3)
(test-zipf make-zipf-generator/zri 131 1.1     6.3    10000 9.5e-3)
(test-zipf make-zipf-generator/zri 131 1.1     20.23  10000 9.5e-3)
(test-zipf make-zipf-generator/zri 131 1.1     41.483 10000 9.5e-3)

(test-end "srfi-194-zipf")
