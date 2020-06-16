(import
  (srfi-194)
  (scheme base)
  (scheme inexact)
  (srfi 27)
  (srfi 64))

(cond-expand
  ((library (srfi 158)) (import (srfi 158)))
  ((library (srfi 121)) (import (srfi 121))))

;; for testing custom rand source,
;; since it's implementation specific
(cond-expand
  (gauche (import 
            (gauche base)
            (math mt-random))))

(define (assert-number-generator/all-in-range gen from to)
  (test-assert
    (generator-every
      (lambda (num)
        (and (>= num from)
             (< num to)))
      (gtake gen 1000))))

(define (assert-number-generator gen from to)
  (define range (- to from))
  (define lower-quarter (+ from (* 0.25 range)))
  (define upper-quarter (- to (* 0.25 range)))
  (assert-number-generator/all-in-range gen from to)
  (test-assert
    (generator-any
      (lambda (num)
        (and (>= num from)
             (< num lower-quarter)))
      gen))
  (test-assert
    (generator-any
      (lambda (num)
        (and (>= num lower-quarter)
             (< num upper-quarter)))
      gen))

  (test-assert
    (generator-any
      (lambda (num)
        (and (>= num upper-quarter)
             (< num to)))
      gen)))

(define (assert-int-generator gen byte-size signed?)
  (define from (if signed?
                   (- (expt 2 (- byte-size 1)))
                   0))
  (define to (if signed?
                 (expt 2 (- byte-size 1))
                 (expt 2 byte-size)))
  (assert-number-generator gen from to))

(test-begin "srfi-194")

(test-group "Test with-random-source syntax"
            ;;create and consume generators that are made with different source
            ;;with various order, and check that order doesn't change the outcome
            (define (test-multiple-sources gen1-maker gen1-expect 
                                           gen2-maker gen2-expect)
              
              ;;create gen1, consume gen1, create gen2, consume gen2
              (let ((gen1 (gen1-maker)))
               (test-equal (generator->list gen1) gen1-expect)
               (let ((gen2 (gen2-maker)))
                (test-equal (generator->list gen2) gen2-expect)))
              
              ;;create gen1, create gen2, consume gen1, consume gen2
              (let ((gen1 (gen1-maker))
                    (gen2 (gen2-maker)))
                (test-equal (generator->list gen1) gen1-expect)
                (test-equal (generator->list gen2) gen2-expect)))
            
            (define multiple-sources-testcase
              (cond-expand
                (gauche (list (lambda ()
                                (gtake (with-random-source
                                       (make <mersenne-twister> :seed 0)
                                       make-random-integer-generator 0 10) 
                                     5))
                              '(5 5 7 8 6)
                              (lambda ()
                                (gtake (with-random-source
                                       (make <mersenne-twister> :seed 1)
                                       make-random-integer-generator 0 10) 
                                     5))
                              '(4 9 7 9 0)))
                (else #f)))
            
            (when multiple-sources-testcase
              (apply test-multiple-sources multiple-sources-testcase))
            
            (with-random-source default-random-source
                                make-random-integer-generator 0 10))

(test-group "Test random int"
            (assert-number-generator
              (make-random-integer-generator 1 100)
              1 100)

            (for-each
              (lambda (testcase)
                (define make-gen (car testcase))
                (define byte-size (cadr testcase))
                (define signed? (caddr testcase))
                (assert-int-generator (make-gen) byte-size signed?))
              (list
                (list make-random-u8-generator 8 #f)
                (list make-random-s8-generator 8 #t)
                (list make-random-u16-generator 16 #f)
                (list make-random-s16-generator 16 #t)
                (list make-random-u32-generator 32 #f)
                (list make-random-s32-generator 32 #t)
                (list make-random-u64-generator 64 #f)
                (list make-random-s64-generator 64 #t)))
            
            ;;test u1 separately, since it will fail quarter checks due to small range
            (assert-number-generator/all-in-range (make-random-u1-generator) 0 2)
            (test-assert
              (generator-any
                (lambda (v) (= v 0))
                (gtake (make-random-u1-generator) 100)))
            (test-assert
              (generator-any
                (lambda (v) (= v 1))
                (gtake (make-random-u1-generator) 100))))

(test-group "Test random real"
            (assert-number-generator
              (make-random-real-generator 1.0 5.0)
              1.0 5.0)

            (test-assert
              (generator-any
                (lambda (v)
                  (not (= v (floor v))))
                (make-random-real-generator 1.0 5.0))))

(test-group "Test random bool"
            (test-assert
              (generator-every
                (lambda (v)
                  (or (eq? v #t)
                      (eq? v #f)))
                (gtake (make-random-boolean-generator) 10000)))

            (test-assert
              (generator-any
                (lambda (v)
                  (eq? #t v))
                (make-random-boolean-generator)))

            (test-assert
              (generator-any
                (lambda (v)
                  (eq? #f v))
                (make-random-boolean-generator))))

(test-group "Test random char"
            (test-assert
              (generator-every
                (lambda (v)
                  (or (equal? v #\a)
                      (equal? v #\b)))
                (gtake (make-random-char-generator "ab")
                       10000)))

            (test-assert
              (generator-any
                (lambda (v)
                  (equal? v #\a))
                (make-random-char-generator "ab")))

            (test-assert
              (generator-any
                (lambda (v)
                  (equal? v #\b))
                (make-random-char-generator "ab"))))

(test-group "Test random string"
            (test-assert
              (generator-every
                (lambda (str)
                  (and (< (string-length str) 5)
                       (every (lambda (c)
                                   (or (equal? c #\a)
                                         (equal? c #\b)))
                                 (string->list str))))
                (gtake (make-random-string-generator 5 "ab")
                       10000)))

            (test-assert
              (generator-any
                (lambda (str)
                  (equal? "abb" str))
                (make-random-string-generator 4 "ab"))))

(test-group "Test Bernoulli"
            (define g (make-bernoulli-generator 0.7))
            (define expect 7000)
            (define actual (generator-count
                                (lambda (i) (= i 0))
                                (gtake g 10000)))
            (define ratio (inexact (/ actual expect)))
            (test-assert (> ratio 0.9))
            (test-assert (< ratio 1.1)))

(test-group "Test poisson"
            ;;TODO import from somewhere?
            (define (fact k)
              (cond
                ((<= k 1) 1)
                (else (* k (fact (- k 1))))))
            (define (expected-fraction L k)
              (/ (* (exact (expt L k)) (exact (exp (- L))))
                 (fact k)))

            (define (test-poisson L poisson-gen test-points)
             (generator-every
               (lambda (k)
                 (define expect (expected-fraction L k))
                 (define actual (/ (generator-count
                                     (lambda (i) (= i k))
                                     (gtake poisson-gen 10000))
                                   10000))
                 (define ratio (inexact (/ actual expect)))
                 (test-assert (> ratio 0.8))
                 (test-assert (< ratio 1.2)))
               (list->generator test-points)))

            (test-poisson 2 (make-poisson-generator 2) '(1 2 3))
            (test-poisson 40 (make-poisson-generator 40) '(30 40 50))
            (test-poisson 280 (make-poisson-generator 280) '(260 280 300)))

(test-group "Test normal"
            (define frac-at-1dev 0.34134)
            (define frac-at-2dev 0.47725)
            (define frac-at-3dev 0.49865)

            (define (test-normal-at-point gen count-from count-to expected-fraction)
              (define actual (/ (generator-count
                                  (lambda (n)
                                    (and (>= n count-from)
                                         (< n count-to)))
                                  (gtake gen 10000))
                                10000.0))
              (test-assert (and (> actual (* 0.9 expected-fraction))
                                (< actual (* 1.1 expected-fraction)))))

            (define (test-normal gen mean deviation)
              (test-normal-at-point gen mean (+ mean deviation) frac-at-1dev)
              (test-normal-at-point gen mean (+ mean (* 2 deviation)) frac-at-2dev)
              (test-normal-at-point gen mean (+ mean (* 3 deviation)) frac-at-3dev))

            (test-normal (make-normal-generator) 0.0 1.0)
            (test-normal (make-normal-generator 1.0) 1.0 1.0)
            (test-normal (make-normal-generator 1.0 2.0) 1.0 2.0))

(test-group "Test exponential"
            (define (expected-fraction mean x)
              (- 1 (exp (* (- (/ 1.0 mean)) x))))

            (define (test-exp-at-point gen count-to expected)
              (define actual (/ (generator-count
                                  (lambda (n)
                                    (< n count-to))
                                  (gtake gen 10000))
                                10000.0))
              (test-assert (> actual (* 0.9 expected)))
              (test-assert (< actual (* 1.1 expected))))

            (define (test-exp gen mean)
              (test-exp-at-point gen 1 (expected-fraction mean 1))
              (test-exp-at-point gen 2 (expected-fraction mean 2))
              (test-exp-at-point gen 3 (expected-fraction mean 3)))

            (test-exp (make-exponential-generator 1) 1)
            (test-exp (make-exponential-generator 1.5) 1.5))

(test-group "Test geometric"
            (define (expected-fraction p x)
              (* (expt (- 1 p) (- x 1)) p))

            (define (test-geom-at-point gen p x)
              (define expected (expected-fraction p x))
              (define actual (/ (generator-count
                                  (lambda (n)
                                    (= n x))
                                  (gtake gen 100000))
                                100000))
              (define ratio (/ actual expected))
              (test-assert (> ratio 0.9))
              (test-assert (< ratio 1.1)))

            (define (test-geom gen p)
              (test-geom-at-point gen p 1)
              (test-geom-at-point gen p 3)
              (test-geom-at-point gen p 5))

            (test-geom (make-geometric-generator 0.5) 0.5))

(test-group "Test uniform sampling"
            (test-equal
              '()
              (generator->list (gsampling)))
            (test-equal
              '()
              (generator->list (gsampling (generator) (generator))))
            (test-equal
              '(1 1 1)
              (generator->list (gsampling (generator) (generator 1 1 1))))
            (test-assert
              (generator-any
                (lambda (el)
                  (= el 1))
                (gsampling (circular-generator 1) (circular-generator 2))))
            (test-assert
              (generator-any
                (lambda (el)
                  (= el 2))
                (gsampling (circular-generator 1) (circular-generator 2)))))

(test-group "Test weighted sampling"
            (test-equal
              '()
              (generator->list (gweighted-sampling)))
            (test-equal
              '()
              (generator->list (gweighted-sampling 1 (generator) 2 (generator))))
            (test-equal
              '(1 1 1)
              (generator->list (gweighted-sampling 1 (generator) 2 (generator 1 1 1))))
            (test-assert
              (generator-any
                (lambda (el)
                  (= el 1))
                (gweighted-sampling 1 (circular-generator 1) 2 (circular-generator 2))))
            (test-assert
              (generator-any
                (lambda (el)
                  (= el 2))
                (gweighted-sampling 1 (circular-generator 1) 2 (circular-generator 2))))
            (let ((expected (/ 1 3))
                  (actual (/ (generator-count
                               (lambda (el) (= el 1))
                               (gtake (gweighted-sampling 1 (circular-generator 1) 2 (circular-generator 2))
                                      10000))
                             10000.0)))
              (test-assert (> actual (* 0.9 expected)))
              (test-assert (< actual (* 1.1 expected)))))

(test-end "srfi-194")

