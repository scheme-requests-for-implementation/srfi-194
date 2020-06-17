;;
;; Parameters & syntax
;;
(define current-random-source (make-parameter default-random-source))

(define-syntax with-random-source
  (syntax-rules ()
    ((_ random-source proc arg ...)
     (begin
       (unless (random-source? random-source)
         (error "expected random source"))
       (parameterize ((current-random-source random-source))
                     (proc arg ...))))))

;;
;; Primitive randoms
;;

(define (make-random-integer-generator low-bound up-bound)
     (when (not (integer? low-bound))
       (error "expected integer"))
     (when (not (integer? up-bound))
       (error "expected integer"))
     (let ((rand-int-proc (random-source-make-integers (current-random-source)))
           (range (- up-bound low-bound)))
       (lambda ()
         (+ low-bound (rand-int-proc range)))))

(define (make-random-u1-generator)
  (make-random-integer-generator 0 2))
(define (make-random-u8-generator) 
  (make-random-integer-generator 0 256))
(define (make-random-s8-generator) 
  (make-random-integer-generator -128 128))
(define (make-random-u16-generator) 
  (make-random-integer-generator 0 65536))
(define (make-random-s16-generator) 
  (make-random-integer-generator -32768 32768))
(define (make-random-u32-generator) 
  (make-random-integer-generator 0 (expt 2 32)))
(define (make-random-s32-generator) 
  (make-random-integer-generator (- (expt 2 31)) (expt 2 31)))
(define (make-random-u64-generator) 
  (make-random-integer-generator 0 (expt 2 64)))
(define (make-random-s64-generator) 
  (make-random-integer-generator (- (expt 2 63)) (expt 2 63)))

(define (make-random-real-generator low-bound up-bound)
  (when (not (number? low-bound))
    (error "expected number"))
  (when (not (number? up-bound))
    (error "expected number"))
  (let ((rand-real-proc (random-source-make-reals (current-random-source)))
        (range (- up-bound low-bound)))
    (lambda ()
      (+ low-bound (* range (rand-real-proc))))))

(define (make-random-complex-generator real-lower-bound imag-lower-bound
                                       real-upper-bound imag-upper-bound)
  (let ((real-gen (make-random-real-generator real-lower-bound real-upper-bound))
        (imag-gen (make-random-real-generator imag-lower-bound imag-upper-bound)))
    (lambda ()
      (make-rectangular (real-gen) (imag-gen)))))

(define (make-random-boolean-generator)
  (define u1 (make-random-u1-generator))
  (lambda ()
    (zero? (u1))))

(define (make-random-char-generator str)
  (when (not (string? str))
    (error "expected string"))
  (let* ((int-gen (make-random-integer-generator 0 (string-length str))))
   (lambda ()
     (string-ref str (int-gen)))))

(define (make-random-string-generator k str)
  (let ((char-gen (make-random-char-generator str))
        (int-gen (make-random-integer-generator 0 k)))
    (lambda ()
      (generator->string char-gen (int-gen)))))

;;
;; Non-uniform distributions
;;

(define PI (* 4 (atan 1.0)))

(define (make-bernoulli-generator p)
  (when (or (< p 0) (> p 1))
    (error "expected 0 <= p <= 1"))
  (let ((rand-real-proc (random-source-make-reals (current-random-source)))) 
   (lambda ()
     (if (<= (rand-real-proc) p)
         0
         1))))

(define (make-categorical-generator pvec)
  (define prob-sum
    (vector-fold 
      (lambda (sum p)
        (unless (and (number? p)
                     (> p 0))
          (error "parameter must be a vector of positive numbers"))
        (+ sum p))
      0
      pvec))
  (unless (= prob-sum 1)
    (error "sum of given probabilities must be equal to 1"))
  (let ((real-gen (make-random-real-generator 0 1)))
   (lambda ()
     (define roll (real-gen))
     (let it ((sum 0)
              (i 0))
       (if (< roll (+ sum (vector-ref pvec i)))
           i
           (it (+ sum (vector-ref pvec i))
               (+ i 1)))))))

(define make-normal-generator
  (case-lambda
    (()
     (make-normal-generator 0.0 1.0))
    ((mean)
     (make-normal-generator mean 1.0))
    ((mean deviation)
     (let ((rand-real-proc (random-source-make-reals (current-random-source)))
           (state #f))
      (lambda ()
        ;;Box-Muller
        (if state
            (let ((result state))
             (set! state #f)
             result)
            (let ((r (sqrt (* -2 (log (rand-real-proc)))))
                  (theta (* 2 PI (rand-real-proc))))
              (set! state (+ mean (* deviation r (cos theta))))
              (+ mean (* deviation r (sin theta))))))))))

(define (make-exponential-generator mean)
  (let ((rand-real-proc (random-source-make-reals (current-random-source))))
   (lambda ()
     (- (* mean (log (rand-real-proc)))))))

(define (make-geometric-generator p)
  (let ((c (/ (log (- 1.0 p))))
        (rand-real-proc (random-source-make-reals (current-random-source))))
    (lambda ()
      (ceiling (* c (log (rand-real-proc)))))))

;; Draw from poisson distribution with mean L, variance L.
;; For small L, we use Knuth's method.  For larger L, we use rejection
;; method by Atkinson, The Computer Generation of Poisson Random Variables,
;; J. of the Royal Statistical Society Series C (Applied Statistics), 28(1),
;; pp29-35, 1979.  The code here is a port by John D Cook's C++ implementation
;; (http://www.johndcook.com/stand_alone_code.html )
(define (make-poisson-generator L)
  (let ((rand-real-proc (random-source-make-reals (current-random-source))))
   (if (< L 36)
       (make-poisson/small rand-real-proc L)
       (make-poisson/large rand-real-proc L))))

;private
(define (make-poisson/small rand-real-proc L)
  (lambda ()
    (do ((exp-L (exp (- L)))
         (k 0 (+ k 1))
         (p 1.0 (* p (rand-real-proc))))
        ((<= p exp-L) (- k 1)))))

;private
(define (make-poisson/large rand-real-proc L)
  (let* ((c (- 0.767 (/ 3.36 L)))
         (beta (/ PI (sqrt (* 3 L))))
         (alpha (* beta L))
         (k (- (log c) L (log beta))))
    (define (loop)
      (let* ((u (rand-real-proc))
             (x (/ (- alpha (log (/ (- 1.0 u) u))) beta))
             (n (exact (floor (+ x 0.5)))))
        (if (< n 0)
            (loop)
            (let* ((v (rand-real-proc))
                   (y (- alpha (* beta x)))
                   (t (+ 1.0 (exp y)))
                   (lhs (+ y (log (/ v (* t t)))))
                   (rhs (+ k (* n (log L)) (- (log-of-fact n)))))
              (if (<= lhs rhs)
                  n
                  (loop))))))
    loop))

;private
;log(n!) table for n 1 to 256. Vector, where nth index corresponds to log((n+1)!)
;Computed on first invocation of `log-of-fact`
(define log-fact-table #f)

;private
;computes log-fact-table
;log(n!) = log((n-1)!) + log(n)
(define (make-log-fact-table!)
   (define table (make-vector 256))
   (vector-set! table 0 0)
   (do ((i 1 (+ i 1)))
       ((> i 255) #t)
       (vector-set! table i (+ (vector-ref table (- i 1))
                               (log (+ i 1)))))
   (set! log-fact-table table))

;private
;returns log(n!)
;adapted from https://www.johndcook.com/blog/2010/08/16/how-to-compute-log-factorial/
(define (log-of-fact n)
  (when (not log-fact-table)
    (make-log-fact-table!))
  (cond
    ((<= n 1) 0)
    ((<= n 256) (vector-ref log-fact-table (- n 1)))
    (else (let ((x (+ n 1)))
           (+ (* (- x 0.5)
                 (log x))
              (- x)
              (* 0.5
                 (log (* 2 PI)))
              (/ 1.0 (* x 12.0)))))))



(define (gsampling . generators-lst)
  (let ((gen-vec (list->vector generators-lst))
        (rand-int-proc (random-source-make-integers (current-random-source))))

       ;remove exhausted generator at index
       (define (remove-gen index)
         (define new-vec (make-vector (- (vector-length gen-vec) 1)))
         ;when removing anything but first, copy all elements before index
         (when (> index 0)
           (vector-copy! new-vec 0 gen-vec 0 index))
         ;when removing anything but last, copy all elements after index
         (when (< index (- (vector-length gen-vec) 1))
           (vector-copy! new-vec index gen-vec (+ 1 index)))
         (set! gen-vec new-vec))

       ;randomly pick generator. If it's exhausted remove it, and pick again
       ;returns value (or eof, if all generators are exhausted)
       (define (pick)
         (let* ((index (rand-int-proc (vector-length gen-vec)))
                (gen (vector-ref gen-vec index))
                (value (gen)))
           (if (eof-object? value)
               (begin
                 (remove-gen index)
                 (if (= (vector-length gen-vec) 0)
                     (eof-object)
                     (pick)))
               value)))

       (lambda ()
         (if (= 0 (vector-length gen-vec))
             (eof-object)
             (pick)))))


