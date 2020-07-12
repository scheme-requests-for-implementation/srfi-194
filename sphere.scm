;
; sphere.scm
; Uniform distributions on a sphere, and a ball.
; Submitted for inclusion in srfi-194
;
; Algorithm based on BoxMeuller as described in
; http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
;


; make-sphere-generator N - return a generator of points uniformly
; distributed on an N-dimensional sphere.
; This implements the BoxMeuller algorithm, that is, of normalizing
; N+1 Gaussian random variables.
(define (make-sphere-generator N)

  (define gaussg (make-normal-generator))

  ; Banach l2-norm aka root-mean-square distance.
  (define (l2-norm VEC)
    (sqrt (vector-fold 
            (lambda (sum x) 
              (+ sum (* x x))) 
            0 VEC)))

  (define np1 (+ N 1))

  ; Create a vector of N+1 values.
  (lambda ()
    (define vect (generator->vector gaussg np1))
    (define norm (/ 1.0 (l2-norm vect)))
    (vector-map (lambda (x) 
                  (* x norm)) 
                vect)))

; make-ball-generator N - return a generator of points uniformly
; distributed inside an N-dimensional ball.
; This implements the Harman-Lacko-Voelker Dropped Coordinate method.
(define (make-ball-generator N)
  (define sphereg (make-sphere-generator (+ N 1)))
  ; Create a vector of N+2 values, and drop the last two.
  ; (The sphere already added one, so we only add one more)
  (lambda () (take (sphereg) N)))

