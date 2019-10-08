(import :gerbil/gambit
        :gerbil/gambit/random
        :std/generic
        :std/misc/list
        :std/misc/repr)

(export #t)

(defstruct vec (x y z)
  constructor: :init!)

(defmethod {:pr vec}
  (lambda (self port options)
    (for-each (cut display <> port)
              ["(vec " (vec-x self) " " (vec-y self) " " (vec-z self) ")"])))

(defmethod {:init! vec}
  (case-lambda
    ((self obj)
     (match obj
       ((vec x y z)
        (set! (vec-x self) x)
        (set! (vec-y self) y)
        (set! (vec-z self) z))
       ([a b]
        (set! (vec-x self) a)
        (set! (vec-y self) b)
        (set! (vec-z self) 0))
       ([a b c]
        (set! (vec-x self) a)
        (set! (vec-y self) b)
        (set! (vec-z self) c))
       ((vector a b)
        (set! (vec-x self) a)
        (set! (vec-y self) b)
        (set! (vec-z self) 0))
       ((vector a b c)
        (set! (vec-x self) a)
        (set! (vec-y self) b)
        (set! (vec-z self) c))))
    ((self x y)
     (set! (vec-x self) x)
     (set! (vec-y self) y)
     (set! (vec-z self) 0))
    ((self x y z)
     (set! (vec-x self) x)
     (set! (vec-y self) y)
     (set! (vec-z self) z))))

; TODO: Provide inexact versions
; TODO: Guard against inf & div by 0

; {add! self vec-rep} -> self!
(defmethod {add! vec}
  (case-lambda
    ((self obj)
     (match obj
       ((vec x y z)
        (set! (vec-x self) (+ (vec-x self) x))
        (set! (vec-y self) (+ (vec-y self) y))
        (set! (vec-z self) (+ (vec-z self) z))
        self)
       ([a b]
        (set! (vec-x self) (+ (vec-x self) a))
        (set! (vec-y self) (+ (vec-y self) b))
        self)
       ([a b c]
        (set! (vec-x self) (+ (vec-x self) a))
        (set! (vec-y self) (+ (vec-y self) b))
        (set! (vec-z self) (+ (vec-z self) c))
        self)
       ((vector a b)
        (set! (vec-x self) (+ (vec-x self) a))
        (set! (vec-y self) (+ (vec-y self) b))
        self)
       ((vector a b c)
        (set! (vec-x self) (+ (vec-x self) a))
        (set! (vec-y self) (+ (vec-y self) b))
        (set! (vec-z self) (+ (vec-z self) c))
        self)))
    ((self x y)
     (set! (vec-x self) (+ (vec-x self) x))
     (set! (vec-y self) (+ (vec-y self) y))
     self)
    ((self x y z)
     (set! (vec-x self) (+ (vec-x self) x))
     (set! (vec-y self) (+ (vec-y self) y))
     (set! (vec-z self) (+ (vec-z self) z))
     self)))

; {sub! self vec-rep} -> self!
(defmethod {sub! vec}
  (case-lambda
    ((self obj)
     (match obj
       ((vec x y z)
        (set! (vec-x self) (- (vec-x self) x))
        (set! (vec-y self) (- (vec-y self) y))
        (set! (vec-z self) (- (vec-z self) z))
        self)
       ([a b]
        (set! (vec-x self) (- (vec-x self) a))
        (set! (vec-y self) (- (vec-y self) b))
        self)
       ([a b c]
        (set! (vec-x self) (- (vec-x self) a))
        (set! (vec-y self) (- (vec-y self) b))
        (set! (vec-z self) (- (vec-z self) c))
        self)
       ((vector a b)
        (set! (vec-x self) (- (vec-x self) a))
        (set! (vec-y self) (- (vec-y self) b))
        self)
       ((vector a b c)
        (set! (vec-x self) (- (vec-x self) a))
        (set! (vec-y self) (- (vec-y self) b))
        (set! (vec-z self) (- (vec-z self) c))
        self)))
    ((self x y)
     (set! (vec-x self) (- (vec-x self) x))
     (set! (vec-y self) (- (vec-y self) y))
     self)
    ((self x y z)
     (set! (vec-x self) (- (vec-x self) x))
     (set! (vec-y self) (- (vec-y self) y))
     (set! (vec-z self) (- (vec-z self) z))
     self)))

; {mul! self scalar} -> self!
(defmethod {mul! vec}
  (lambda (self s)
    (set! (vec-x self) (* (vec-x self) s))
    (set! (vec-y self) (* (vec-y self) s))
    (set! (vec-z self) (* (vec-z self) s))
    self))

; {div! self scalar} -> self!
(defmethod {div! vec}
  (lambda (self s)
    (when (not (= s 0))
      (set! (vec-x self) (/ (vec-x self) s))
      (set! (vec-y self) (/ (vec-y self) s))
      (set! (vec-z self) (/ (vec-z self) s)))
    self))

; {mag-sq self} -> squared-magnitude
(defmethod {mag-sq vec}
  (lambda (self)
    (let ((x (vec-x self))
          (y (vec-y self))
          (z (vec-z self)))
      (+ (* x x)
         (* y y)
         (* z z)))))

; {mag self} -> magnitude
(defmethod {mag vec}
  (lambda (self)
    (sqrt {mag-sq self})))

; {dot self vec-rep} -> dot-product
(defmethod {dot vec}
  (case-lambda
    ((self obj)
     (def (dot x y (z 0))
       (+ (* (vec-x self) x)
          (* (vec-y self) y)
          (* (vec-z self) z)))
     (match obj
       ((vec x y z)
        (dot x y z))
       ([a b]
        (dot a b))
       ([a b c]
        (dot a b c))
       ((vector a b)
        (dot a b))
       ((vector a b c)
        (dot a b c))))
    ((self x y)
     (+ (* (vec-x self) x)
        (* (vec-y self) y)
        (* (vec-z self) 0)))
    ((self x y z)
     (+ (* (vec-x self) x)
        (* (vec-y self) y)
        (* (vec-z self) z)))))

; {cross self vec-rep} -> cross-product
(defmethod {cross vec}
    (case-lambda
     ((self obj)
      (def (cross x y (z 0))
        (vec (- (* (vec-y self) z) (* (vec-z self) y))
             (- (* (vec-z self) x) (* (vec-x self) z))
             (- (* (vec-x self) y) (* (vec-y self) x))))
      (match obj
        ((vec x y z)
         (cross x y z))
        ([a b]
         (cross a b))
        ([a b c]
         (cross a b c))
        ((vector a b)
         (cross a b))
        ((vector a b c)
         (cross a b c))))
     ((self x y)
      (vec (- (* (vec-y self) 0) (* (vec-z self) y))
           (- (* (vec-z self) x) (* (vec-x self) 0))
           (- (* (vec-x self) y) (* (vec-y self) x))))
     ((self x y z)
      (vec (- (* (vec-y self) z) (* (vec-z self) y))
           (- (* (vec-z self) x) (* (vec-x self) z))
           (- (* (vec-x self) y) (* (vec-y self) x))))))

; {dist self vec-rep} -> distance
(defmethod {dist vec}
  (case-lambda
    ((self obj)
     (def (dist v)
       {mag {sub! v self}})
     (match obj
       ((vec x y z)
        (dist (vec x y z)))
       ([a b]
        (dist (vec a b)))
       ([a b c]
        (dist (vec a b c)))
       ((vector a b)
        (dist (vec a b)))
       ((vector a b c)
        (dist (vec a b c)))))
    ((self x y)
     {mag {sub! (vec x y) self}})
    ((self x y z)
     {mag {sub! (vec x y z) self}})))

; {normalize! self} -> self!
(defmethod {normalize! vec}
  (lambda (self)
    (let (len {mag self})
      (when (not (= len 0))
        ; Uses mul! to avoid dup div by 0 check
        {mul! self (/ 1 len)}))
    self))

; {limit! self max-mag} -> self!
(defmethod {limit! vec}
  (lambda (self lim-mag)
    (let (squared-mag {mag-sq self})
      (when (> squared-mag (* lim-mag lim-mag))
        {mul! {div! self (sqrt squared-mag)} lim-mag}))
    self))

; {mag! self new-mag} -> self!
(defmethod {mag! vec}
  (lambda (self scalar)
    {mul! {normalize! self} scalar}))

; {heading self} -> xy-angle-of-rotation
(defmethod {heading vec}
  (lambda (self)
    (atan (vec-y self) (vec-x self))))

; {rotate self angle} -> self!
(defmethod {rotate vec}
  (lambda (self a)
    (let ((new-heading (+ {heading self} a))
          (mag {mag self}))
      (set! (vec-x self) (* (cos new-heading) mag))
      (set! (vec-y self) (* (sin new-heading) mag)))
    self))

; {angle-between self vec-rep} -> angle
(defmethod {angle-between vec}
  (case-lambda
    ((self obj)
     (def (a-between v)
       (let (dot-mag-mag (/ {dot self v} (* {mag self} {mag v})))
         (acos (min 1 (max -1 dot-mag-mag)))))
     (match obj
       ((? vec?)
        (a-between obj))
       ([a b]
        {angle-between self a b})
       ([a b c]
        {angle-between self a b c})
       ((vector a b)
        {angle-between self a b})
       ((vector a b c)
        {angle-between self a b c})))
    ((self x y)
     {angle-between self (vec x y)})
    ((self x y z)
     {angle-between self (vec x y z)})))

; {lerp! self vec-rep} -> self!
(defmethod {lerp! vec}
  (case-lambda
    ((self obj amnt)
     (def (lerp! x y z amount)
       (set! (vec-x self) (+ (vec-x self) (* (- x (vec-x self)) amount)))
       (set! (vec-y self) (+ (vec-y self) (* (- y (vec-y self)) amount)))
       (set! (vec-z self) (+ (vec-z self) (* (- z (vec-z self)) amount))))
     (match obj
       ((vec a b c)
        (lerp! a b c amnt)
        self)
       ([a b]
        (lerp! a b 0 amnt)
        self)
       ([a b c]
        (lerp! a b c amnt)
        self)
       ((vector a b)
        (lerp! a b 0 amnt)
        self)
       ((vector a b c)
        (lerp! a b c amnt)
        self)))
    ((self x y amnt)
     (set! (vec-x self) (+ (vec-x self) (* (- x (vec-x self)) amnt)))
     (set! (vec-y self) (+ (vec-y self) (* (- y (vec-y self)) amnt)))
     (set! (vec-z self) (+ (vec-z self) (* (- 0 (vec-z self)) amnt)))
     self)
    ((self x y z amnt)
     (set! (vec-x self) (+ (vec-x self) (* (- x (vec-x self)) amnt)))
     (set! (vec-y self) (+ (vec-y self) (* (- y (vec-y self)) amnt)))
     (set! (vec-z self) (+ (vec-z self) (* (- z (vec-z self)) amnt)))
     self)))

; ->vector
; ->list
; ->values

;;;;;; Redo below

(def (vec-add a b)
  (vec (+ (vec-x a) (vec-x b))
       (+ (vec-y a) (vec-y b))
       (+ (vec-z a) (vec-z b))))

(def (vec-sub a b)
  (vec (- (vec-x a) (vec-x b))
       (- (vec-y a) (vec-y b))
       (- (vec-z a) (vec-z b))))

(def (vec-mul v s)
  (vec (* (vec-x v) s)
       (* (vec-y v) s)
       (* (vec-z v) s)))

(def (vec-div v s)
  (vec (/ (vec-x v) s)
       (/ (vec-y v) s)
       (/ (vec-z v) s)))

; (vec-equal)


; vec-from-angle
; vec-from-angles
; vec-random-2d
; vec-random-3d
; vec-add
; vec-sub
; vec-mul
; vec-div
