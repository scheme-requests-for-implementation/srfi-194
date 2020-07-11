;
; Create a Zipf random distribution.
;
; Created by Linas Vepstas 10 July 2020
; Nominated for inclusion in srfi-194
;
; Not optimized for speed!
;
; Implementation from ZRI algorithm presented in the paper:
; "Rejection-inversion to generate variates from monotone discrete
; distributions", Wolfgang Hörmann and Gerhard Derflinger
; ACM TOMACS 6.3 (1996): 169-184
;
; Hörmann and Derflinger use "q" everywhere, when they really mean "s".
; Thier "q" is not the standard q-series deformation. Its just "s".
; The notation in the code below has been changed to reflect
; conventional usage.
;
; ------------------------------------------------------------

; The Hurwicz zeta distribution 1 / (k+q)^s for 1 <= k <= n integer
; The Zipf distribution is recovered by setting q=0.
;
; The exponent `s` must be a real number not equal to 1.
; Accuracy is diminished for |1-s|< 1e-6. The accuracy is roughly
; equal to 1e-15 / |1-s| where 1e-15 == 64-bit double-precision ULP.
;
; Example usage:
;    (define zgen (make-zipf-generator 50 1.01 0))
;    (generator->list zgen 10)
;
(define (make-zipf-generator/zri n s q)

	; The hat function h(x) = 1 / (x+q)^s
	(define (hat x s)
		(expt (+ x q) (- s)))

	; The integral of hat(x)
	; H(x) = (x+q)^{1-s} / (1-s)
	; Note that H(x) is always negative.
	(define (big-h x s)
		(define 1ms (- 1 s))
		(/ (expt (+ q x) 1ms) 1ms))

	; The inverse function of H(x)
	(define (big-h-inv x s)
		(define 1ms (- 1 s))
		(define oms (/ 1 1ms))
		(- (expt (* x 1ms) oms) q))

	; Lower and upper bounds for the uniform random generator.
	; Note that both are negative for all values of s.
	(define big-h-half (big-h 0.5 s))
	(define big-h-n (big-h (+ n 0.5) s))

	; Rejection cut
	(define cut (- 1 (big-h-inv (- (big-h 1.5 s) (expt (+ 1 q) (- s))) s)))

	; Uniform distribution
	(define dist (make-random-real-generator big-h-half big-h-n))

	; Attempt to hit the dartboard. Return #f if we fail,
	; otherwise return an integer between 1 and n.
	(define (try)
		(define u (dist))
		(define x (big-h-inv u s))
		(define kflt (floor (+ x 0.5)))
		(define k (inexact->exact kflt))
		(if (or
			(<= (- k x) cut)
			(>= u (- (big-h (+ k 0.5) s) (hat k s)))) k #f))

	; Did we hit the dartboard? If not, try again.
	(define (loop-until)
		(define k (try))
		(if k k (loop-until)))

	; Return the generator.
	loop-until
)

; The Hurwicz zeta distribution 1 / (k+q)^s for 1 <= k <= n integer
; The Zipf distribution is recovered by setting q=0.
;
; The exponent `s` must be a real number close to 1.
; Accuracy is diminished for |1-s|> 2e-4. The accuracy is roughly
; equal to 0.05 * |1-s|^4 due to exp(1-s) being expanded to 4 terms.
;
; This handles the special case of s==1 perfectly.
(define (make-zipf-generator/one n s q)

	; The hat function h(x) = 1 / (x+q)^s
	; Written for s->1 i.e. 1/(x+q)(x+q)^{s-1}
	(define (hat x)
		(define xpq (+ x q))
		(/ (expt xpq (- 1 s)) xpq))

	; Expansion of exp(1-s)/(1-s) for s->1
	(define 1ms (- 1 s))
	(define (trm n u) (+ 1 (/ (* 1ms u) n)))
	(define exn (trm 2 (trm 3 (trm 4 1))))

	; The integral of hat(x)
	; H(x) = log (x+q)
	(define (big-h x)
		(* (log (+ q x)) exn))

	; The inverse function of H(x)
	(define (big-h-inv x)
		(- (exp (/ x exn)) q))

	; Lower and upper bounds for the uniform random generator.
	; Note that both are negative for all values of s.
	(define big-h-half (big-h 0.5))
	(define big-h-n (big-h (+ n 0.5)))

	; Rejection cut
	(define cut (- 1 (big-h-inv (- (big-h 1.5) (/ 1 (+ 1 q))))))

	; Uniform distribution
	(define dist (make-random-real-generator big-h-half big-h-n))

	; Attempt to hit the dartboard. Return #f if we fail,
	; otherwise return an integer between 1 and n.
	(define (try)
		(define u (dist))
		(define x (big-h-inv u))
		(define kflt (floor (+ x 0.5)))
		(define k (inexact->exact kflt))
		(if (or
			(<= (- k x) cut)
			(>= u (- (big-h (+ k 0.5)) (hat k)))) k #f))

	; Did we hit the dartboard? If not, try again.
	(define (loop-until)
		(define k (try))
		(if k k (loop-until)))

	; Return the generator.
	loop-until
)

*unspecified*
