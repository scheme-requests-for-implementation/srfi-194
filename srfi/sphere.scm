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
(define (make-sphere-generator arg)
  (cond
    ((integer? arg) (make-sphere-generator* (make-vector (+ 1 arg) 1.0)))
    ((vector? arg) (make-sphere-generator* arg))
    (else (error "expected argument to either be a number (dimension), or vector (axis length for the dimensions)"))))

(define (make-sphere-generator* dim-sizes)
  (define gaussg-vec
    (vector-map
      (lambda (size)
        (make-normal-generator 0.0 size))
      dim-sizes))
  ; Banach l2-norm aka root-mean-square distance.
  (define (l2-norm VEC)
    (sqrt (vector-fold
            (lambda (sum x l)
              (+ sum (/ (* x x)
                        (* l l)
                        )))
            0
            VEC
            dim-sizes)))

  (lambda ()
    (define vect
      (vector-map
        (lambda (gaussg)
          (gaussg))
        gaussg-vec))
    (define norm (/ 1.0 (l2-norm vect)))
    (vector-map (lambda (x)
                  (* x norm))
                vect)))

; make-ball-generator N - return a generator of points uniformly
; distributed inside an N-dimensional ball.
; This implements the Harman-Lacko-Voelker Dropped Coordinate method.
(define (make-ball-generator arg)
  (define dim-sizes
    (cond
      ((integer? arg) (make-vector (+ 2 arg) 1.0))
      ((vector? arg) (vector-append arg (vector 1.0 1.0)))
      (else (error "expected argument to either be a number (dimension), or vector (axis length for the dimensions)"))))
  (define N (- (vector-length dim-sizes) 2))
  (define sphereg (make-sphere-generator (+ N 2)))
  ; Create a vector of N+2 values, and drop the last two.
  ; (The sphere already added one, so we only add one more)
  (lambda ()
    (vector-map
      (lambda (el dim-size _) (* el dim-size))
      (sphereg)
      dim-sizes
      (make-vector N #f))))

