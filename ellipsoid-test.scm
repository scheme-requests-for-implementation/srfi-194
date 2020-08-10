;
; ellipsoid-test.scm
;
; Verify that the distribution of points on the surface of an
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

; Compute sum of a list of numbers
(define (sum lst) (fold (lambda (x sum) (+ sum x)) 0 lst))

; -----------------------------------------------------------
; Stuff for the complete elliptic integral
(define pi 3.14159265358979)

; factorial
(define (fact n rv)
   (if (zero? n) rv (fact (- n 1) (* n rv))))

; Double factorial
; https://en.wikipedia.org/wiki/Double_factorial
(define (double-fact n rv)
   (if (<= n 0) rv (double-fact (- n 2) (* n rv))))

; Complete elliptic integral per wikipedia, see the Ivorty& Bessel
; expansion. Here `a` and `b` are the axes.
; https://en.wikipedia.org/wiki/Ellipse
(define (complete-elliptic a b)
   (define rh (/ (- a b) (+ a b)))
   (define h (* rh rh))

	(define precision 1e-10)

   (define (ivory term n twon hn fact-n dfact-n sum)
      (if (< term precision) (+ sum term)
         ; (format #t "yo n= ~A term=~A 2^n=~A h^n=~A n!=~A n!!=~A sum=~A\n"
         ; n term twon hn fact-n dfact-n sum)
         (ivory
            (/ (* dfact-n dfact-n hn) (* twon twon fact-n fact-n))
            (+ n 1)
            (* 2 twon)
            (* h hn)
            (* (+ n 1) fact-n)
            (* (- (* 2 n) 1) dfact-n)
            (+ term sum))))

   (* pi (+ a b) (+ 1 (/ h 4)
      (ivory (/ (* h h) 64) 3 8 (* h h h) 6 3 0.0))))

; -----------------------------------------------------------

; Assumes that `points` is a list of 2D vectors of floats.
(define (verify-ellipse points)
	; Place in sorted order.
	(define ordered-points (clockwise points))

	; Difference between neghboring points.
	(define diffs (delta ordered-points '()))

	; Compute the distances between neighboring points
	(define dists (map l2-norm diffs))

	; Sum of the intervals
	(define perimeter (sum dists))

	; Find major and minor axes
	(define major
		(fold (lambda (MAJ x)
			(if (< MAJ x) x MAJ))
			0
			(map l2-norm points)))

	(define minor
		(fold (lambda (MIN x)
			(if (< MIN x) MIN x))
			1.0e308
			(map l2-norm points)))

	; The expected perimiter
	(define perim-exact (complete-elliptic major minor))

	; The normalized difference of measured and expected perimeters
	; Should almost always be less than ten, often less than two.
	(define error
		(abs (* (/ (- perimeter perim-exact) perim-exact) (length points))))

	(format #t "Number of points: ~A\n" (length points))
	(format #t "Measured perimeter: ~A\n" perimeter)
	(format #t "Major and minor axes: ~A ~A\n" major minor)
	(format #t "Expected perimeter: ~A\n" perim-exact)
	(format #t "Relative error: ~A\n" error)
)
