;;; :patternspandemic/delaunator unit tests

(import :std/test
        :std/sugar
        "delaunator")

(export delaunator-test)

(def delaunator-test
  (test-suite "Delaunator"

    (def (failing-check thing)
      (check thing => #f))

    (test-case "This should fail."
      (failing-check 7))

  )
)
