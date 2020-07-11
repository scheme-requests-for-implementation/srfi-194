
; Guile-specific test infrastructure
(use-modules (srfi srfi-64))


; Debug wrapper
(define (xmake-zipf-generator n s)
	(define foo (make-zipf-generator n s))
	(lambda () (foo)))

;(define (xmake-zipf-generator n s)
;	(define foo (make-zipf-generator/zri n s))
;	(lambda () (foo)))

; Debug utility for gnuplot
(define (vector-to-file vec filename)
	(let ((outport (open-file filename "w")))
   	(vector-for-each
			(lambda (i x) (format outport "~A	~A\n" (+ i 1) x))
			vec)
      (close outport)))

; Take REPS samples from the zipf distribution (NVOCAB ESS)
;
; and whaaaatever
;
(define (test-zipf NVOCAB ESS REPS TOL)

	; Bin-counter containing accumulated histogram.
	(define bin-counts (make-vector NVOCAB 0))

	; Accumulate samples into the histogram.
	(generator-for-each
		(lambda (SAMP)
			(define offset (- SAMP 1))
			(vector-set! bin-counts offset (+ 1 (vector-ref bin-counts offset))))
		(gtake (xmake-zipf-generator NVOCAB ESS) REPS))

	; The total counts in the bins should be equal to REPS
	(test-assert
		(equal? REPS
			(vector-fold
				(lambda (i sum cnt) (+ sum cnt)) 0 bin-counts)))

	; Verify the distribution is within tolerance.

	; Frequency is normalized to be 0.0 to 1.0
	(define frequency (vector-map (lambda (i n) (/ n REPS)) bin-counts))
	(define probility (vector-map (lambda (i n) (exact->inexact n)) frequency))

	; Sequence 1..NVOCAB
	(define seq
		(vector-unfold (lambda (i x) (values x (+ x 1))) NVOCAB 1))

	; Sequence  1/n^ESS
	(define inv-pow (vector-map (lambda (i n) (expt n (- ESS))) seq))

	; Harmonic number sum_1..NVOCAB 1/n^ESS
	(define hnorm
		(vector-fold
			(lambda (i sum cnt) (+ sum cnt)) 0 inv-pow))

	; The expected distribution
	(define expect
		(vector-map (lambda (i x) (/ x hnorm)) inv-pow))

	(define prexpect (vector-map (lambda (i x) (exact->inexact x)) expect))

	(format #t "Norm = ~A\n" (exact->inexact hnorm))
	(vector-to-file probility "foo.dat")

	#f
)

(test-begin "srfi-194-zipf")
(test-group "Test zipf 30 1"
	(test-assert #f)
)
(test-end "srfi-194-zipf")
