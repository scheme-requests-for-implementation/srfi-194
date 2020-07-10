(define-library (srfi-194)
  (import (scheme base)
          (scheme case-lambda)
          (scheme inexact)
          (scheme complex)
          (srfi 133) ;; vector-fold
          (srfi 27))

  (cond-expand
    ((library (srfi 158)) (import (srfi 158)))
    ((library (srfi 121)) (import (srfi 121))))

  (export
    
    clamp-real-number
    
    current-random-source
    with-random-source
    
    make-random-integer-generator
    make-random-u1-generator
    make-random-u8-generator make-random-s8-generator
    make-random-u16-generator make-random-s16-generator
    make-random-u32-generator make-random-s32-generator
    make-random-u64-generator make-random-s64-generator
    make-random-boolean-generator
    make-random-char-generator
    make-random-string-generator
    make-random-real-generator
    make-random-complex-generator

    make-bernoulli-generator 
    make-categorical-generator
    make-normal-generator
    make-exponential-generator
    make-geometric-generator
    make-poisson-generator

    gsampling)

  (include "srfi-194-impl.scm"))
