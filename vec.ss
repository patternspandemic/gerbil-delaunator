; TODO: Provide inexact versions, vec->inexact, vec->exact
; TODO: Guard against inf & div by 0

; Possible additions:
; - inexact equal-accuracy
; - proportional predicate
; - symmetric vec
; - component projections
; - max/min component
; - ratios
; - abs vec
; - inverse vec
; - 

(import :gerbil/gambit/random
        :std/misc/repr)

(export #t)

(defstruct vec (x y z)
  constructor: :init!
  equal: #t
  print: #t)

(defmethod {:pr vec}
  (lambda (self port options)
    (for-each (cut display <> port)
              ["(vec " (vec-x self) " " (vec-y self) " " (vec-z self) ")"])))


;;; vec methods


; {:init! self vec-rep} -> #!void
(defmethod {:init! vec}
  (case-lambda
    ((self obj)
     (match obj
       ((vec x y z) {:init! self x y z})
       ([a b] {:init! self a b 0})
       ([a b c] {:init! self a b c})
       ((vector a b) {:init! self a b 0})
       ((vector a b c) {:init! self a b c})))
    ((self x y) {:init! self x y 0})
    ((self x y z)
     (set! (vec-x self) x)
     (set! (vec-y self) y)
     (set! (vec-z self) z))))

; {set! self vec-rep} -> self!
(defmethod {set! vec}
  (case-lambda
    ((self obj)
     (match obj
       ((vec x y z) {set! self x y z})
       ([a b] {set! self a b 0})
       ([a b c] {set! self a b c})
       ((vector a b) {set! self a b 0})
       ((vector a b c) {set! self a b c})))
    ((self x y) {set! self x y 0})
    ((self x y z)
     (set! (vec-x self) x)
     (set! (vec-y self) y)
     (set! (vec-z self) z)
     self)))

; {add! self vec-rep} -> self!
(defmethod {add! vec}
  (case-lambda
    ((self obj)
     (match obj
       ((vec x y z) {add! self x y z})
       ([a b] {add! self a b})
       ([a b c] {add! self a b c})
       ((vector a b) {add! self a b})
       ((vector a b c) {add! self a b c})))
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
       ((vec x y z) {sub! self x y z})
       ([a b] {sub! self a b})
       ([a b c] {sub! self a b c})
       ((vector a b) {sub! self a b})
       ((vector a b c) {sub! self a b c})))
    ((self x y)
     (set! (vec-x self) (- (vec-x self) x))
     (set! (vec-y self) (- (vec-y self) y))
     self)
    ((self x y z)
     (set! (vec-x self) (- (vec-x self) x))
     (set! (vec-y self) (- (vec-y self) y))
     (set! (vec-z self) (- (vec-z self) z))
     self)))

; {mult! self scalar} -> self!
(defmethod {mult! vec}
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
     (match obj
       ((vec x y z) {dot self x y z})
       ([a b] {dot self a b})
       ([a b c] {dot self a b c})
       ((vector a b) {dot self a b})
       ((vector a b c) {dot self a b c})))
    ((self x y)
     (+ (* (vec-x self) x)
        (* (vec-y self) y)))
    ((self x y z)
     (+ (* (vec-x self) x)
        (* (vec-y self) y)
        (* (vec-z self) z)))))

; {cross self vec-rep} -> cross-product
(defmethod {cross vec}
  (case-lambda
    ((self obj)
     (match obj
       ((vec x y z) {cross self x y z})
       ([a b] {cross self a b})
       ([a b c] {cross self a b c})
       ((vector a b) {cross self a b})
       ((vector a b c) {cross self a b c})))
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
     (match obj
       ((vec x y z) {dist self x y z})
       ([a b] {dist self a b})
       ([a b c] {dist self a b c})
       ((vector a b) {dist self a b})
       ((vector a b c) {dist self a b c})))
    ((self x y) {mag {sub! (vec x y) self}})
    ((self x y z) {mag {sub! (vec x y z) self}})))

; {normalize! self} -> self!
(defmethod {normalize! vec}
  (lambda (self)
    (let (len {mag self})
      (when (not (= len 0))
        ; Uses mult! to avoid dup div by 0 check
        {mult! self (/ 1 len)}))
    self))

; {limit! self max-mag} -> self!
(defmethod {limit! vec}
  (lambda (self lim-mag)
    (let (squared-mag {mag-sq self})
      (when (> squared-mag (* lim-mag lim-mag))
        {mult! {div! self (sqrt squared-mag)} lim-mag}))
    self))

; {mag! self new-mag} -> self!
(defmethod {mag! vec}
  (lambda (self scalar)
    {mult! {normalize! self} scalar}))

; {heading self} -> xy-angle-of-rotation
(defmethod {heading vec}
  (lambda (self)
    (atan (vec-y self) (vec-x self))))

; {rotate! self angle} -> self!
(defmethod {rotate! vec}
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
       ((? vec?) (a-between obj))
       ([a b] {angle-between self a b})
       ([a b c] {angle-between self a b c})
       ((vector a b) {angle-between self a b})
       ((vector a b c) {angle-between self a b c})))
    ((self x y) {angle-between self (vec x y)})
    ((self x y z) {angle-between self (vec x y z)})))

; {lerp! self vec-rep} -> self!
(defmethod {lerp! vec}
  (case-lambda
    ((self obj amnt)
     (match obj
       ((vec x y z) {lerp! self x y z amnt})
       ([a b] {lerp! self a b 0 amnt})
       ([a b c] {lerp! self a b c amnt})
       ((vector a b) {lerp! self a b 0 amnt})
       ((vector a b c) {lerp! self a b c amnt})))
    ((self x y amnt) {lerp! self x y 0 amnt})
    ((self x y z amnt)
     (set! (vec-x self) (+ (vec-x self) (* (- x (vec-x self)) amnt)))
     (set! (vec-y self) (+ (vec-y self) (* (- y (vec-y self)) amnt)))
     (set! (vec-z self) (+ (vec-z self) (* (- z (vec-z self)) amnt)))
     self)))


;;; vec procs


; (vec-add v1 v2 [target = #f]) -> vec
(def (vec-add v1 v2 (target #f))
  (if (not target)
    (set! target (vec v1))
    {set! target v1})
  {add! target v2})

; (vec-sub v1 v2 [target = #f]) -> vec
(def (vec-sub v1 v2 (target #f))
  (if (not target)
    (set! target (vec v1))
    {set! target v1})
  {sub! target v2})

; (vec-mult v s [target = #f]) -> vec
(def (vec-mult v s (target #f))
  (if (not target)
    (set! target (vec v))
    {set! target v})
  {mult! target s})

; (vec-div v s [target = #f]) -> vec
(def (vec-div v s (target #f))
  (if (not target)
    (set! target (vec v))
    {set! target v})
  {div! target s})

; (vec-lerp v1 v2 amnt [target = #f]) -> vec
(def (vec-lerp v1 v2 amnt (target #f))
  (if (not target)
    (set! target (vec v1))
    {set! target v1})
  {lerp! target v2 amnt})

; (vec-from-angle angle [mag = 1]) -> vec
(def (vec-from-angle angle (mag 1))
  (vec (* (cos angle) mag)
       (* (sin angle) mag)
       0))

; (vec-from-angles theta phi [mag = 1]) -> vec
(def (vec-from-angles theta phi (mag 1))
  (let ((cos-phi (cos phi))
        (sin-phi (sin phi))
        (cos-theta (cos theta))
        (sin-theta (sin theta)))
    (vec (* mag sin-theta sin-phi)
         (* -1 mag cos-theta)
         (* mag sin-theta cos-phi))))
  
; (vec-random-2d) -> vec
(def (vec-random-2d)
  (vec-from-angle (* (random-real) 3.141592653589793 2)))

; (vec-random-3d) -> vec
(def (vec-random-3d)
  (let* ((angle (* (random-real) 3.141592653589793 2))
         (vz (1- (* (random-real) 2)))
         (vz-base (sqrt (- 1 (* vz vz))))
         (vx (* vz-base (cos angle)))
         (vy (* vz-base (sin angle))))
    (vec vx vy vz)))
