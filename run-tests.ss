#!/usr/bin/env gxi

(import :std/test
        "delaunator-test")

(def tests
  [delaunator-test])

(apply run-tests! tests)
(test-report-summary!)

(case (test-result)
  ((OK) (exit 0))
  (else (exit 1)))
