; Up to date with mapbox/Delaunator at 35223e3182bc1c87d7d6d03ef56b0e1dd8e3eef9 Aug. 30, 2019
; Some of these loops can may be able to be replaced with while/until.

(import :std/iter
        :std/logger
        :std/srfi/113
        :std/srfi/128
        :std/sugar
        :gerbil/gambit/bits
        :gerbil/gambit/exact
        :gerbil/gambit/hvectors)

(export #t)
; TODO: Probably remove slot accessors and provide iterators for them.
; (export
;   Delaunator-coords
;   Delaunator-triangles
;   Delaunator-halfedges
;   Delaunator-hull
;   ;bounds
;   delaunator/from
;   halfedge-ids-of-triangle
;   triangle-id-of-edge
;   next-halfedge
;   prev-halfedge
;   ;point-ids-of-triangle
;   ;triangle-ids-adjacent-to-triangle
;   ;iter-triangle-edges
;   ;iter-triangles
;   ;triangle-center
;   ;iter-voronoi-edges
;   ;edges-around-point
;   ;iter-voronoi-regions
;   )


(def EPSILON (expt 2 -52))
(def EDGE_STACK (make-u32vector 512))


(def (circumcenter ax ay bx by cx cy)
  (let* ((dx (- bx ax))
         (dy (- by ay))
         (ex (- cx ax))
         (ey (- cy ay))
         (bl (+ (* dx dx) (* dy dy)))
         (cl (+ (* ex ex) (* ey ey)))
         (d  (/ 0.5 (- (* dx ey) (* dy ex))))
         (x (+ ax (* (- (* ey bl) (* dy cl)) d)))
         (y (+ ay (* (- (* dx cl) (* ex bl)) d))))
    (values x y)))

(def (circumradius ax ay bx by cx cy)
  (let* ((dx (- bx ax))
         (dy (- by ay))
         (ex (- cx ax))
         (ey (- cy ay))
         (bl (+ (* dx dx) (* dy dy)))
         (cl (+ (* ex ex) (* ey ey)))
         (d  (/ 0.5 (- (* dx ey) (* dy ex))))
         (x (* (- (* ey bl) (* dy cl)) d))
         (y (* (- (* dx cl) (* ex bl)) d)))
    (+ (* x x) (* y y))))

(def (dist ax ay bx by)
  (let ((dx (- ax bx))
        (dy (- ay by)))
    (+ (* dx dx) (* dy dy))))

(def (orient-if-sure px py rx ry qx qy)
  (let ((l (* (- ry py) (- qx px)))
        (r (* (- rx px) (- qy py))))
    (if (>= (abs (- l r)) (* 3.3306690738754716e-16 (abs (+ l r))))
        (- l r)
        #f))) ;0)))

(def (orient rx ry qx qy px py)
  (< (or (orient-if-sure px py rx ry qx qy)
         (orient-if-sure rx ry qx qy px py)
         (orient-if-sure qx qy px py rx ry)
         0)
      0))

(def (quicksort ids dists left right)

  (def (swap vec i j)
    (let (tmp (u32vector-ref vec i))
      (u32vector-set! vec i (u32vector-ref vec j))
      (u32vector-set! vec j tmp)))

  (if (<= (- right left) 20)

      ; then
      (do-while ((i (1+ left) (1+ i)))
                ((<= i right))
        (let* ((temp (u32vector-ref ids i))
               (temp-dist (f64vector-ref dists temp))
               (j (1- i)))
          (while (and (>= j left)
                      (> (f64vector-ref dists (u32vector-ref ids j)) temp-dist))
            (u32vector-set! ids (1+ j) (u32vector-ref ids j))
            (set! j (1- j)))
          (u32vector-set! ids (1+ j) temp)))

      ; else
      (let ((median (arithmetic-shift (+ left right) -1))
            (i (1+ left))
            (j right))
        (swap ids median i)
        (if (> (f64vector-ref dists (u32vector-ref ids left))
               (f64vector-ref dists (u32vector-ref ids right)))
            (swap ids left right))
        (if (> (f64vector-ref dists (u32vector-ref ids i))
               (f64vector-ref dists (u32vector-ref ids right)))
            (swap ids i right))
        (if (> (f64vector-ref dists (u32vector-ref ids left))
               (f64vector-ref dists (u32vector-ref ids i)))
            (swap ids left i))

        (let* ((temp (u32vector-ref ids i))
               (temp-dist (f64vector-ref dists temp)))

          (call/cc (lambda (esc-break)
            (let lp ()
              (let lpi ()
                (set! i (1+ i))
                (if (< (f64vector-ref dists (u32vector-ref ids i)) temp-dist)
                    (lpi)))
              (let lpi ()
                (set! j (1- j))
                (if (> (f64vector-ref dists (u32vector-ref ids j)) temp-dist)
                    (lpi)))
              (if (< j i) (esc-break))
              (swap ids i j)
              (lp))
          )) ; end breaking call/cc

          (u32vector-set! ids (1+ left) (u32vector-ref ids j))
          (u32vector-set! ids j temp)

          (if (>= (- right (1+ i)) (- j left))
              ; then
              (begin (quicksort ids dists i right)
                     (quicksort ids dists left (1- j)))
              ; else
              (begin (quicksort ids dists left (1- j))
                     (quicksort ids dists i right)))))))

(def (pseudo-angle dx dy)
  (let* ((p (/ dx (+ (abs dx) (abs dy))))
         (n (if (> dy 0) (- 3 p) (1+ p))))
    (/ n 4)))

(def (in-circle ax ay bx by cx cy px py)
  (let* ((dx (- ax px))
         (dy (- ay py))
         (ex (- bx px))
         (ey (- by py))
         (fx (- cx px))
         (fy (- cy py))
         (ap (+ (* dx dx) (* dy dy)))
         (bp (+ (* ex ex) (* ey ey)))
         (cp (+ (* fx fx) (* fy fy))))
    (< (+ (- (* dx (- (* ey cp) (* bp fy)))
             (* dy (- (* ex cp) (* bp fx))))
             (* ap (- (* ex fy) (* ey fx))))
        0)))


(defclass Delaunator
  (coords           ; initial flattened coords
   min-X            ; computed min x of coords
   min-Y            ; computed min y of coords
   max-X            ; computed max x of coords
   max-Y            ; computed max y of coords
   triangles        ; trimmed version of _triangles
   halfedges        ; trimmed version of _halfedges
   hull             ; halfedges of the final convex hull
   triangles-length ; length to trim to
   _triangles       ; triangles indexed by halfedge
   _halfedges       ; halfedges indexed by their complement
   _hash-size       ; computed size for _hull-hash
   _hull-start      ; temporary hull start
   _hull-prev       ; edge to prev edge
   _hull-next       ; edge to next edge
   _hull-tri        ; edge to adjacent triangle
   _hull-hash       ; angular edge hash
   _ids             ; helper for sorting points
   _dists           ; helper for sorting points
   _cx              ; temporary circumcenter x
   _cy              ; temporary circumcenter y
   _point-to-leftmost-halfedge-index) ; for fast lookup of leftmost incoming halfedge
  constructor: :init!)


(defmethod {:init! Delaunator}
  (lambda (self coordinates)
    ; (unless (f64vector? coordinates)
    ;         (error "Expected coordinates to be f64vector"))
    (let* ((n (arithmetic-shift (f64vector-length coordinates) -1)) ; num of coords
           (max-triangles (max (- (* 2 n) 5) 0))
           (hash-size (exact (ceiling (sqrt n)))))
      (set! (@ self coords) coordinates)
      (set! (@ self min-X) +inf.0)
      (set! (@ self min-Y) +inf.0)
      (set! (@ self max-X) -inf.0)
      (set! (@ self max-Y) -inf.0)
      (set! (@ self triangles-length) 0)
      ; Vectors that will store the triangulation graph
      (set! (@ self _triangles) (make-u32vector (* max-triangles 3)))
      (set! (@ self _halfedges) (make-s32vector (* max-triangles 3)))
      ; Temporary vectors for tracking the edges of the advancing convex hull
      (set! (@ self _hash-size) hash-size)
      (set! (@ self _hull-prev) (make-u32vector n))
      (set! (@ self _hull-next) (make-u32vector n))
      (set! (@ self _hull-tri) (make-u32vector n))
      (set! (@ self _hull-hash) (make-s32vector hash-size -1))
      ; Temporary vectors for sorting points
      (set! (@ self _ids) (make-u32vector n))
      (set! (@ self _dists) (make-f64vector n))
      ; Index points to their leftmost incoming halfedge
      ; TODO: Tune hash table options...
      (set! (@ self _point-to-leftmost-halfedge-index) (make-hash-table)))
    {update self}))


; TODO: update is not written to support a change in the number of coords.
;       Support it or guard against it.
(defmethod {update Delaunator}
  (lambda (self)
    (call/cc
      (lambda (escape)

        (let ((coords    (@ self coords))
              (hull-prev (@ self _hull-prev))
              (hull-next (@ self _hull-next))
              (hull-tri  (@ self _hull-tri))
              (hull-hash (@ self _hull-hash))
              (n (arithmetic-shift (f64vector-length (@ self coords)) -1)))

          ; Find min/max bounds of input coords
          (do-while ((i 0 (1+ i)))
                    ((< i n))
            (let ((x (f64vector-ref coords (* 2 i)))
                  (y (f64vector-ref coords (1+ (* 2 i)))))
              (if (< x (@ self min-X)) (set! (@ self min-X) x))
              (if (< y (@ self min-Y)) (set! (@ self min-Y) y))
              (if (> x (@ self max-X)) (set! (@ self max-X) x))
              (if (> y (@ self max-Y)) (set! (@ self max-Y) y))
              (u32vector-set! (@ self _ids) i i))) ; Also init _ids

          (let ((cx (/ (+ (@ self min-X) (@ self max-X)) 2))
                (cy (/ (+ (@ self min-Y) (@ self max-Y)) 2))
                (min-dist +inf.0)
                (i0 0)
                (i1 0)
                (i2 0))

            ; Pick a seed point close to the center
            (do-while ((i 0 (1+ i)))
                      ((< i n))
              (let (d (dist cx cy (f64vector-ref coords (* 2 i)) (f64vector-ref coords (1+ (* 2 i)))))
                (when (< d min-dist)
                  (set! i0 i)
                  (set! min-dist d))))

            (let ((i0x (f64vector-ref coords (* 2 i0)))
                  (i0y (f64vector-ref coords (1+ (* 2 i0)))))
              (set! min-dist +inf.0) ; reset min-dist measure

              ; Find the point closest to the seed i0
              (do-while ((i 0 (1+ i)))
                        ((< i n))
                (unless (= i i0)
                  (let (d (dist i0x i0y (f64vector-ref coords (* 2 i)) (f64vector-ref coords (1+ (* 2 i)))))
                    (when (and (< d min-dist)
                              (> d 0))
                      (set! i1 i)
                      (set! min-dist d)))))

              (let ((i1x (f64vector-ref coords (* 2 i1)))
                    (i1y (f64vector-ref coords (1+ (* 2 i1))))
                    (min-radius +inf.0))

                ; Find the 3rd point which forms the smallest circumcircle with the first two
                (do-while ((i 0 (1+ i)))
                          ((< i n))
                  (unless (or (= i i0) (= i i1))
                    (let (r (circumradius i0x i0y i1x i1y (f64vector-ref coords (* 2 i)) (f64vector-ref coords (1+ (* 2 i)))))
                      (when (< r min-radius)
                        (set! i2 i)
                        (set! min-radius r)))))

                (let ((i2x (f64vector-ref coords (* 2 i2)))
                      (i2y (f64vector-ref coords (1+ (* 2 i2)))))

                  (when (= min-radius +inf.0)
                    ; order collinear points by dx (or dy if all are identical)
                    ; and return the list as a hull
                    (do-while ((i 0 (1+ i)))
                              ((< i n))
                      (let* ((c0 (f64vector-ref coords 0))
                             (c1 (f64vector-ref coords 1))
                             (dx (- (f64vector-ref coords (* 2 i)) c0))
                             (dy (- (f64vector-ref coords (1+ (* 2 i))) c1)))
                        (f64vector-set! (@ self _dists) i (if (zero? dx) dy dx))))
                    (quicksort (@ self _ids) (@ self _dists) 0 (1- n))

                    (let ((temp-hull (make-u32vector n))
                          (j 0)
                          (d0 -inf.0))

                      (do-while ((i 0 (1+ i)))
                                ((< i n))
                        (let* ((id (u32vector-ref (@ self _ids) i))
                               (ref-dists-id (f64vector-ref (@ self _dists) id)))
                          (when (> ref-dists-id d0)
                            (u32vector-set! temp-hull j id)
                            (set! j (1+ j))
                            (set! d0 ref-dists-id))))
                      (set! (@ self hull) (subu32vector temp-hull 0 j))
                      (set! (@ self triangles) (make-u32vector 0))
                      (set! (@ self halfedges) (make-s32vector 0))
                      (escape))) ; End special case early escape (2 coords?)

                  ; Swap the order of the seed points for counter-clockwise orientation
                  (when (orient i0x i0y i1x i1y i2x i2y)
                    (let ((i i1)
                          (x i1x)
                          (y i1y))
                      (set! i1 i2)
                      (set! i1x i2x)
                      (set! i1y i2y)
                      (set! i2 i)
                      (set! i2x x)
                      (set! i2y y)))

                  (let-values (((center-x center-y) (circumcenter i0x i0y i1x i1y i2x i2y)))
                    (set! (@ self _cx) center-x)
                    (set! (@ self _cy) center-y)
                    (do-while ((i 0 (1+ i)))
                              ((< i n))
                      (f64vector-set! (@ self _dists) i (dist (f64vector-ref coords (* 2 i)) (f64vector-ref coords (1+ (* 2 i))) center-x center-y))))

                  ; Sort the points by distance from the seed triangle circumcenter
                  (quicksort (@ self _ids) (@ self _dists) 0 (1- n))

                  ; Set up the seed triangle as the starting hull
                  (set! (@ self _hull-start) i0)

                  (u32vector-set! hull-next i0 i1)
                  (u32vector-set! hull-prev i2 i1)
                  (u32vector-set! hull-next i1 i2)
                  (u32vector-set! hull-prev i0 i2)
                  (u32vector-set! hull-next i2 i0)
                  (u32vector-set! hull-prev i1 i0)

                  (u32vector-set! hull-tri i0 0)
                  (u32vector-set! hull-tri i1 1)
                  (u32vector-set! hull-tri i2 2)

                  (s32vector-fill! hull-hash -1)
                  (s32vector-set! hull-hash {_hash-key self i0x i0y} i0)
                  (s32vector-set! hull-hash {_hash-key self i1x i1y} i1)
                  (s32vector-set! hull-hash {_hash-key self i2x i2y} i2)

                  (set! (@ self triangles-length) 0)
                  {_add-triangle self i0 i1 i2 -1 -1 -1}

                  (let ((hull-size 3)
                        (xp +inf.0)
                        (yp +inf.0))

                    ; The main update loop...
                    (let lp ((k 0))
                      (when (< k (u32vector-length (@ self _ids)))
                        (call/cc (lambda (esc-continue)

                          (let* ((i (u32vector-ref (@ self _ids) k))
                                 (x (f64vector-ref coords (* 2 i)))
                                 (y (f64vector-ref coords (1+ (* 2 i)))))

                            ; Skip near-duplicate points
                            (if (and (> k 0)
                                    (<= (abs (- x xp)) EPSILON)
                                    (<= (abs (- y yp)) EPSILON))
                              (esc-continue))

                            (set! xp x)
                            (set! yp y)

                            ; Skip seed triangle points
                            (if (or (= i i0) (= i i1) (= i i2))
                              (esc-continue))

                            ; Find a visible edge on the convex hull using edge hash
                            (let ((start 0)
                                  (key {_hash-key self x y}))
                              (let lpi ((j 0))
                                (when (< j (@ self _hash-size))
                                  (set! start (s32vector-ref hull-hash (modulo (+ key j) (@ self _hash-size))))
                                  (unless (and (not (= start -1))
                                               (not (= start (u32vector-ref hull-next start))))
                                    (lpi (1+ j)))))

                              (set! start (u32vector-ref hull-prev start))

                              (let* ((e start)
                                     (q (u32vector-ref hull-next e)))
                                (call/cc (lambda (esc-break)
                                  (let lpi ()
                                    (let (orient-result (orient x y
                                                                (f64vector-ref coords (* 2 e)) (f64vector-ref coords (1+ (* 2 e)))
                                                                (f64vector-ref coords (* 2 q)) (f64vector-ref coords (1+ (* 2 q)))))
                                      (when (not orient-result)
                                        (set! e q)
                                        (when (= e start)
                                          (set! e -1)
                                          (esc-break))
                                        (set! q (u32vector-ref hull-next e))
                                        (lpi))))))

                                ; Likely a near-duplicate point, skip it
                                (if (= e -1) (esc-continue))

                                ; Add the first triangle from the point
                                (let (t {_add-triangle self e i (u32vector-ref hull-next e) -1 -1 (u32vector-ref hull-tri e)})

                                  ; Recursively flip triangles from the point until they satisfy the Delaunay condition
                                  (u32vector-set! hull-tri i {_legalize self (+ t 2)})
                                  (u32vector-set! hull-tri e t) ; Keep track of boundary triangles on the hull
                                  (set! hull-size (1+ hull-size))

                                  ; Walk forward through the hull, adding more triangles and flipping recursively
                                  (let (n (u32vector-ref hull-next e))
                                    (set! q (u32vector-ref hull-next n))
                                    (let lpi ()
                                      (let (orient-result (orient x y
                                                            (f64vector-ref coords (* 2 n)) (f64vector-ref coords (1+ (* 2 n)))
                                                            (f64vector-ref coords (* 2 q)) (f64vector-ref coords (1+ (* 2 q)))))
                                        (when orient-result
                                          (set! t {_add-triangle self n i q (u32vector-ref hull-tri i) -1 (u32vector-ref hull-tri n)})
                                          (u32vector-set! hull-tri i {_legalize self (+ t 2)})
                                          (u32vector-set! hull-next n n) ; Mark as removed
                                          (set! hull-size (1- hull-size))
                                          (set! n q)
                                          (set! q (u32vector-ref hull-next n))
                                          (lpi))))

                                    ; Walk backward from the other side, adding more triangles and flipping
                                    (when (= e start)
                                      (set! q (u32vector-ref hull-prev e))
                                      (let lpi ()
                                        (let (orient-result (orient x y
                                                              (f64vector-ref coords (* 2 q)) (f64vector-ref coords (1+ (* 2 q)))
                                                              (f64vector-ref coords (* 2 e)) (f64vector-ref coords (1+ (* 2 e)))))
                                          (when orient-result
                                            (set! t {_add-triangle self q i e -1 (u32vector-ref hull-tri e) (u32vector-ref hull-tri q)})
                                            {_legalize self (+ t 2)}
                                            (u32vector-set! hull-tri q t)
                                            (u32vector-set! hull-next e e) ; Mark as removed
                                            (set! hull-size (1- hull-size))
                                            (set! e q)
                                            (set! q (u32vector-ref hull-prev e))
                                            (lpi)))))

                                    ; Update the hull indices
                                    (set! (@ self _hull-start) e)
                                    (u32vector-set! hull-prev i e)
                                    (u32vector-set! hull-next e i)
                                    (u32vector-set! hull-prev n i)
                                    (u32vector-set! hull-next i n)

                                    ; Save the two new edges in the hash table
                                    (s32vector-set! hull-hash {_hash-key self x y} i)
                                    (s32vector-set! hull-hash {_hash-key self (f64vector-ref coords (* 2 e)) (f64vector-ref coords (1+ (* 2 e)))} e)

                        ))))))) ; end continuing call/cc
                        (lp (1+ k))
                    )) ; end main update loop

                    (let ((triangles-length (@ self triangles-length)))

                      (set! (@ self hull) (make-u32vector hull-size))
                      (do-while ((i 0 (1+ i))
                                 (e (@ self _hull-start) (u32vector-ref hull-next e)))
                                ((< i hull-size))
                        (u32vector-set! (@ self hull) i e))

                      ; Build the index of point id to leftmost incoming halfedge,
                      ; useful for retrieving adhoc voronoi regions.
                      ; TODO: Clear to old size for more efficiency?
                      (hash-clear! (@ self _point-to-leftmost-halfedge-index)) ; reset index from prev update
                      (do-while ((e 0 (1+ e)))
                                ((< e (@ self triangles-length)))
                        (let ((endpoint (u32vector-ref (@ self _triangles) (next-halfedge e)))
                              (index (@ self _point-to-leftmost-halfedge-index)))
                          (if (or (not (hash-key? index endpoint))
                                  (= (s32vector-ref (@ self _halfedges) endpoint) -1))
                            (hash-put! index endpoint e))))

                      ; Trim typed triangle mesh vectors
                      (set! (@ self triangles) (subu32vector (@ self _triangles) 0 triangles-length))
                      (set! (@ self halfedges) (subs32vector (@ self _halfedges) 0 triangles-length)))

    )))))))) ; end escaping call/cc
)) ; end update method


(defmethod {_hash-key Delaunator}
  (lambda (self x y)
    (let ((cx (@ self _cx))
          (cy (@ self _cy))
          (hs (@ self _hash-size)))
      (modulo
        (exact (floor (* (pseudo-angle (- x cx) (- y cy)) hs)))
        hs))))


(defmethod {_legalize Delaunator}
  (lambda (self a)
    (let ((triangles (@ self _triangles))
          (halfedges (@ self _halfedges))
          (coords (@ self coords))
          (i 0)
          (ar 0))

      (call/cc
        (lambda (esc-break)

          (let lp ()
            (let ((b (s32vector-ref halfedges a))
                  (a0 (- a (modulo a 3))))

              (set! ar (+ a0 (modulo (+ a 2) 3)))

              (when (= b -1) ; convex hull edge
                (if (= i 0) (esc-break))
                (set! i (1- i))
                (set! a (u32vector-ref EDGE_STACK i))
                (lp))

              (let* ((b0 (- b (modulo b 3)))
                     (al (+ a0 (modulo (+ a 1) 3)))
                     (bl (+ b0 (modulo (+ b 2) 3)))
                     (p0 (u32vector-ref triangles ar))
                     (pr (u32vector-ref triangles a))
                     (pl (u32vector-ref triangles al))
                     (p1 (u32vector-ref triangles bl))
                     (illegal (in-circle
                                (f64vector-ref coords (* 2 p0)) (f64vector-ref coords (1+ (* 2 p0)))    ; ax ay
                                (f64vector-ref coords (* 2 pr)) (f64vector-ref coords (1+ (* 2 pr)))    ; bx by
                                (f64vector-ref coords (* 2 pl)) (f64vector-ref coords (1+ (* 2 pl)))    ; cx cy
                                (f64vector-ref coords (* 2 p1)) (f64vector-ref coords (1+ (* 2 p1)))))) ; px py

                (if illegal
                  ; then
                  (begin
                    (u32vector-set! triangles a p1)
                    (u32vector-set! triangles b p0)

                    (let ((hbl (s32vector-ref halfedges bl))
                          (br (+ b0 (modulo (1+ b) 3))))

                      (when (= hbl -1) ; edge swapped on other side of hull (rare) -fix halfedge reference
                        (call/cc
                          (lambda (done)
                            (let lpi ((e (@ self _hull-start)))
                              (when (= (u32vector-ref (@ self _hull-tri) e) bl)
                                (u32vector-set! (@ self _hull-tri) e a)
                                (done)
                              (set! e (u32vector-ref (@ self _hull-prev) e)))
                              (unless (= e (@ self _hull-start))
                                (lpi e))))))

                      {_link self a hbl}
                      {_link self b (s32vector-ref halfedges ar)}
                      {_link self ar bl}

                      (when (< i (u32vector-length EDGE_STACK))
                        ; "Don't worry about hitting the cap: it can only happen on extremely degenerate input."
                        (u32vector-set! EDGE_STACK i br)
                        (set! i (1+ i)))))
                  ; else
                  (begin
                    (if (= i 0) (esc-break))
                    (set! i (1- i))
                    (set! a (u32vector-ref EDGE_STACK i))))
                (lp)))) ; end let lp

      )) ; end breaking call/cc

      ar ; return ar!
))) ; end _legalize method


(defmethod {_link Delaunator}
  (lambda (self a b)
    (let (halfedges (@ self _halfedges))
      (s32vector-set! halfedges a b)
      (unless (= b -1) (s32vector-set! halfedges b a)))))


(defmethod {_add-triangle Delaunator}
  (lambda (self i0 i1 i2 a b c)
    (let ((t (@ self triangles-length))
          (triangles (@ self _triangles)))
      (u32vector-set! triangles t i0)
      (u32vector-set! triangles (+ t 1) i1)
      (u32vector-set! triangles (+ t 2) i2)
      {_link self t a}
      {_link self (+ t 1) b}
      {_link self (+ t 2) c}
      (set! (@ self triangles-length) (+ t 3))
      t))) ; return position of added triangle


; Delaunator -> min-X min-Y max-X max-Y
; The min-X min-Y max-X max-Y values describing the rectangular bounds of the triangulation.
(defmethod {bounds Delaunator}
  (lambda (self)
    (values (@ self min-X) (@ self min-Y)
            (@ self max-X) (@ self max-Y))))


; TODO: support taking vector or list?
(def (delaunator/from points ; vector of pairs
                      get/x: (get/x car)
                      get/y: (get/y cadr))
  (let* ((n (vector-length points))
         (coords (make-f64vector (* n 2))))
    (for ((p points)
          (k (in-range n)))
      (f64vector-set! coords (* 2 k) (inexact (get/x p)))
      (f64vector-set! coords (1+ (* 2 k)) (inexact (get/y p))))
    (make-Delaunator coords)))


;;; Helper procs below are not written for performance :/

; triangle-id -> List-of-halfedge-ids
; The halfedge ids of a triangle with id `t`.
(def (halfedge-ids-of-triangle t)
  (list (* 3 t) (1+ (* 3 t)) (+ 2 (* 3 t))))

; halfedge-id -> triangle-id
; The id of the triangle for which halfedge with id `e` is a part. (also the id of the 1st point?)
(def (triangle-id-of-edge e)
  (floor (/ e 3)))

; halfedge-id -> halfedge-id
; The id of the next halfedge of the triangle for which halfedge with id `e` is a part.
(def (next-halfedge e)
  (if (= (modulo e 3) 2)
      (- e 2)
      (1+ e)))

; halfedge-id -> halfedge-id
; The id of the previous halfedge of the triangle for which halfedge with id `e` is a part.
(def (prev-halfedge e)
  (if (= (modulo e 3) 0)
      (+ e 2)
      (1- e)))

; Delaunator triangle-id -> List-of-point-ids
; The point ids composing the triangle with id `t`.
(defmethod {point-ids-of-triangle Delaunator}
  (lambda (self t)
    (let (pids (map (lambda (e) (u32vector-ref (@ self triangles) e))
                    (halfedge-ids-of-triangle t)))
      (values (car pids) (cadr pids) (caddr pids)))))

; Delaunator triangle-id -> List-of-triangle-ids
; The triangle ids adjacent to triangle with id `t`.
(defmethod {triangle-ids-adjacent-to-triangle Delaunator}
  (lambda (self t)
    (filter-map (lambda (e)
      (let (opposite (s32vector-ref (@ self halfedges) e))
        (if (>= opposite 0)
            (triangle-id-of-edge opposite)
            #f)))
      (halfedge-ids-of-triangle t))))


; Delaunator -> (#!void -> point-id f64vector)
; Provides an iterator yielding values for each point of the triangulation.
; The values yielded are the id of the point, and a f64vector describing the
; point's location.
(defmethod {iter-points Delaunator}
  (lambda (self)
    (lambda ()
      (let (coords (@ self coords))
        (do-while ((p 0 (1+ p)))
                  ((< p (/ (f64vector-length coords) 2)))
          (yield p (subf64vector coords (* 2 p) (+ (* 2 p) 2))))))))


; TODO: Attempt to replace yielded p with the id of the same point in coords
; Delaunator -> (#!void -> point-id f64vector)
; Provides an iterator yielding values for each point of the triangulation's
; hull. The values yielded are the id of the point, and a f64vector describing
; the point's location.
(defmethod {iter-hull-points Delaunator}
  (lambda (self)
    (lambda ()
      (let ((coords (@ self coords))
            (hull (@ self hull)))
        (do-while ((p 0 (1+ p)))
                  ((< p (u32vector-length hull)))
          (yield p (subf64vector coords (* 2 (u32vector-ref hull p)) (+ (* 2 (u32vector-ref hull p)) 2))))))))

; TODO: Attempt to replace yielded e with the halfedge-id of the edge p -> q, but how to find? May require additional index.
; Delaunator -> (#!void -> ???-id f64vector f64vector)
; Provides an iterator yielding values for each edge of the triangulation's hull.
; The values yielded are the id of the halfedge running from p to q, a f64vector
; describing the point the edge starts at, and a f64vector describing the point
; the edge ends at.
(defmethod {iter-hull-edges Delaunator}
  (lambda (self)
    (lambda ()
      (let* ((coords (@ self coords))
             (hull (@ self hull))
             (hull-next (@ self _hull-next)))
            ;  (index (@ self _point-to-leftmost-halfedge-index))) ; BROKEN?
        (do-while ((e 0 (1+ e)))
                  ((< e (u32vector-length hull)))
          (let* ((pid (u32vector-ref hull e))
                 (qid (u32vector-ref hull-next pid))
                 ;(halfedge-id (hash-get index qid)) ; BROKEN
                 (p (subf64vector coords (* 2 pid) (+ (* 2 pid) 2)))
                 (q (subf64vector coords (* 2 qid) (+ (* 2 qid) 2))))
            ; TODO: Confirm triangles[halfedge-id] is same as pid
            (yield e p q)))))))


; Delaunator -> (#!void -> halfedge-id f64vector f64vector)
; Provides an iterator yielding values for each edge of the triangulation.
; The values yielded are the id of the halfedge chosen for the edge, a
; f64vector describing the point the edge starts at, and a f64vector
; describing the point the edge ends at.
(defmethod {iter-triangle-edges Delaunator}
  (lambda (self)
    (lambda ()
      (let ((coords (@ self coords))
            (triangles (@ self triangles)))
        (do-while ((e 0 (1+ e)))
                  ((< e (@ self triangles-length)))
          (when (> e (s32vector-ref (@ self halfedges) e))
            (let* ((pid (u32vector-ref triangles e))
                   (qid (u32vector-ref triangles (next-halfedge e)))
                   (p (subf64vector coords (* 2 pid) (+ (* 2 pid) 2)))
                   (q (subf64vector coords (* 2 qid) (+ (* 2 qid) 2))))
              (yield e p q))))))))


; Delaunator -> (#!void -> triangle-id f64vector f64vector f64vector)
; Provides an iterator yielding values for each triangle of the triangulation.
; The values yielded are the id of the triangle, and three f64vector, each 
; describing a point of the triangle.
(defmethod {iter-triangles Delaunator}
  (lambda (self)
    (lambda ()
      (do-while ((t 0 (1+ t)))
                ((< t (/ (@ self triangles-length) 3)))
        (let-values (((pid1 pid2 pid3) {point-ids-of-triangle self t})
                     ((coords) (@ self coords)))
          (let ((p1 (subf64vector coords (* 2 pid1) (+ (* 2 pid1) 2)))
                (p2 (subf64vector coords (* 2 pid2) (+ (* 2 pid2) 2)))
                (p3 (subf64vector coords (* 2 pid3) (+ (* 2 pid3) 2))))
            (yield t p1 p2 p3)))))))


; Delaunator triangle-id -> f64vector
; The circumcenter of triangle with id `t`.
(defmethod {triangle-center Delaunator}
  (lambda (self t)
    (let-values (((pid1 pid2 pid3) {point-ids-of-triangle self t})
                 ((coords) (@ self coords)))
      (let* ((p1 (subf64vector coords (* 2 pid1) (+ (* 2 pid1) 2)))
             (p2 (subf64vector coords (* 2 pid2) (+ (* 2 pid2) 2)))
             (p3 (subf64vector coords (* 2 pid3) (+ (* 2 pid3) 2)))
             ((values cx cy) (circumcenter (f64vector-ref p1 0) (f64vector-ref p1 1)
                                           (f64vector-ref p2 0) (f64vector-ref p2 1)
                                           (f64vector-ref p3 0) (f64vector-ref p3 1))))
        (f64vector cx cy)))))


; Delaunator -> (#!void -> halfedge-id f64vector f64vector)
; Provides an iterator yielding values for each bisecting voronoi edge of the
; graph dual to the triangulation. The values yielded are the id of the
; halfedge chosen for the bisected edge, a f64vector describing the
; circumcenter of the triangle for which that halfedge is a part, and a
; f64vector describing the circumcenter of the adjacent triangle for which that
; halfedge's compliment is a part. By default, no voronoi edges bisecting the
; hull of the triangulation are provided. For such bisectors to be included, a
; bounding region must be provided to clip these edges against.
; TODO: Support voronoi edges of the hull.
(defmethod {iter-voronoi-edges Delaunator}
  (lambda (self) ;(lambda (self bounds: (bounds #f)) ...)
    (lambda ()
      (do-while ((e 0 (1+ e)))
                ((< e (@ self triangles-length)))
        (when (< e (s32vector-ref (@ self halfedges) e)) ; excludes halfedges on hull
          (let* ((p {triangle-center self (triangle-id-of-edge e)})
                 (q {triangle-center self (triangle-id-of-edge (s32vector-ref (@ self halfedges) e))}))
            (yield e p q)))))))


; Delaunator halfedge-id -> List-of-halfedge-ids
; The ids of all halfedges pointing to the point that `start` points to.
; `start` should not point to a point on the hull.
(defmethod {edge-ids-around-point Delaunator}
  (lambda (self start)
    (let ((incoming start)
          (result []))
      (let lp ()
        (set! result (cons incoming result))
        (let (outgoing (next-halfedge incoming))
            (set! incoming (s32vector-ref (@ self halfedges) outgoing)))
        (unless (or (= incoming -1)
                    (= incoming start))
          (lp)))
      result)))


; Delaunator -> (#!void -> point-id f64vector-of-f64vector)
; Provides an iterator yielding values for each region of the voronoi diagram.
; The values yielded are the id of the point to which the region belongs, and
; a list of f64vectors, describing the region's polygon.
; TODO: Support voronoi regions of the hull.
(defmethod {iter-voronoi-regions Delaunator}
  (lambda (self)  ;(lambda (self bounds: (bounds #f)) ...)
    (let* ((seen (set (make-default-comparator)))
           (hull-ids (@ self hull))
           (hull-length (u32vector-length hull-ids)))
      ; Do not yield regions of hull points by default
      (for ((i (in-range hull-length)))
        (set-adjoin! seen (u32vector-ref hull-ids i)))
      (lambda ()
        (do-while ((e 0 (1+ e)))
                  ((< e (@ self triangles-length)))
          (let ((p (u32vector-ref (@ self triangles) (next-halfedge e))))
            (unless (set-contains? seen p)
              (set-adjoin! seen p)
              (let* ((edge-ids {edge-ids-around-point self e})
                     (triangle-ids (map triangle-id-of-edge edge-ids))
                     (vertices (map (lambda (tid) {triangle-center self tid}) triangle-ids)))
                (yield p vertices)))))))))

(defmethod {voronoi-region Delaunator}
  (lambda (self point-id)
    (let* ((incoming-halfedge (hash-get (@ self _point-to-leftmost-halfedge-index) point-id))
           (p (u32vector-ref (@ self triangles) (next-halfedge incoming-halfedge)))
           (edge-ids {edge-ids-around-point self incoming-halfedge})
           (triangle-ids (map triangle-id-of-edge edge-ids))
           (vertices (map (lambda (tid) {triangle-center self tid}) triangle-ids)))
      (values p vertices))))


; {iter-triangle-centers} ?