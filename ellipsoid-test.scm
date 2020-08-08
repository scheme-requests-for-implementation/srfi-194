; 
; ellipsoid-test.scm
;
; Verify that the diistribution of points on the surface of an
; ellipsoid is uniform.
;
; Test proceeds by taking 2-D slices through the ellipsoid, and
; verifying uniformity on that slice. Thus, the core test is for
; ellipses.
;

; Sort a list of 2D vectors of floats into clock-wise order.
; Assumes that `pts` is a list of 2D vectors of floats.
(define (clockwise pts)
	(sort pts (lambda (a b)
		(if (and (< 0 (vector-ref a 1)) (< 0 (vector-ref b 1)))
			(< (vector-ref b 0) (vector-ref a 0))
			(if (and (< (vector-ref a 1) 0) (< (vector-ref b 1) 0))
				(< (vector-ref a 0) (vector-ref b 0))
				(< (vector-ref b 1) (vector-ref a 1)))))))

; Verfiy that the routine above is not broken.
; Returns #t if it is OK.
(define (test-clockwise)
	(define  clock (list
		'#(1 1e-3) '#(0.8 0.2) '#(0.2 0.8)
		'#(0 1) '#(-0.2 0.8) '#(-0.8 0.2) '#(-1 1e-3)
		'#(-1 -1e-3) '#(-0.8 -0.2) '#(-0.2 -0.8)
		'#(0 -1) '#(0.2 -0.8) '#(0.8 -0.2) '#(1 -1e-3)))

	(equal? (clockwise clock) clock))

; Vector subtraction
; Example usage: (vector-diff '#( 2 3) '#(0.5 0.7))
(define (vector-sub a b)
   (vector-map (lambda (idx ea eb) (- ea eb)) a b))

; Newton differences - compute the difference between neighboring
; points. Assumes `pts` is a list of vectors.  Should be called with
; `rv` set to the null list. (tail-recursive helper) 
(define (delta pts rv)
   (if (null? (cdr pts)) (reverse! rv)
      (delta (cdr pts) (cons (vector-diff (car pts) (cadr pts)) rv))))

; Assumes that `points` is a list of 2D vectors of floats.
(define (verify-ellipse points)
	; Place in sorted order.
	(define ordered-points (clockwise points))

	; Difference between neghboring points.
	(define diffs (delta ordered-points '()))

	; Compute the distances between neighboring points
	(define dists (map l2-norm diffs))

