;
; Create a Zipf random distribution.
;
; Created by Linas Vepstas 10 July 2020
; Nominated for inclusion in srfi-194
; Not yet unit-tested, may be buggy!
;
; Implementation from ZRI algorithm presented in the paper:
; "Rejection-inversion to generate variates from monotone discrete
; distributions", Wolfgang Hörmann and Gerhard Derflinger
; ACM TOMACS 6.3 (1996): 169-184
;
; Example usage:
;    (define zgen (make-zipf-generator 50 1.01))
;    (generator->list zgen 10)
;
;
; ------------------------------------------------------------
; Begin with private functions, not to be exported.
; Public API at bottom.

; The Hurwicz zeta offfset, called "v" in the original paper.
; Can be any value greater than zero.
(define vee 1)

; The hat function h(x) = 1 / (x+v) ^ q)
(define (hat x q)
	(expt (+ x vee) (- q))
)

; The integral of hat(x)
; H(x) = (v+x)^{1-q} / (1-q)
; Note that H(x) is always negative.
(define (big-h x q)
	(define 1mq (- 1 q))
	(/ (expt (+ vee x) 1mq) 1mq)
)

; The inverse function of H(x)
(define (big-h-inv x q)
	(define 1mq (- 1 q))
	(define omq (/ 1 1mq))
	(- (expt (* x 1mq) omq) vee)
)

; ------------------------------------------------------------

; The public API
; q must be greater than 1
(define (make-zipf-generator/zri n q)
	(define big-h-half (big-h 0.5 q))
	(define big-h-n (big-h (+ n 0.5) q))

	(define cut (- 1 (big-h-inv (- (big-h 1.5 q) (expt (+ 1 vee) (- q))) q)))

	(define dist (make-random-real-generator big-h-half big-h-n))

	; Attempt to hit the dartboard. Return #f if we fail,
	; otherwise return an integer between 1 and n.
	(define (try)
		(define u (dist))
		(define x (big-h-inv u q))
		(define k (floor (+ x 0.5)))
		(if (or
			(<= (- k x) cut)
			(>= u (- (big-h (+ k 0.5) q) (hat k q)))) k #f))

	; Did we hit the dartboard? If not, try again.
	(define (loop-until)
		(define k (try))
		(if k k (loop-until)))

	; Return the generator.
	loop-until
)

*unspecified*
