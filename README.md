# bestfit-ujk

Exactly verified Lyapunov certificates for the average-case waste of **Best Fit** bin packing under the
discrete uniform distribution **U{j,k}** (bin capacity k, item sizes i.i.d. uniform on {1..j}): for which
(j,k) is the expected waste O(1), and for which is it Θ(n)?

Results (all backed by certificates in `certs/`):

* the classification is complete for every k ≤ 25, including the seven pairs left open by
  Coffman–Johnson–Shor–Weber (STOC 1993);
* for every j ≤ 12 and every k: linear waste iff 8 ≤ j and j+3 ≤ k ≤ K_j, with K_8..K_12 = 11, 13, 17, 20, 24;
* certified intervals for the waste rate of the linear cases with k ≤ 14; in 26 linear cases, all linear
  waste is in bins with exactly one unit free.

`results/SUMMARY.md` lists every pair with its certificate; `notes/method.md` has the theorems and proofs.

## Use

```
make            # build (C++17 and C; no other dependencies)
make test       # ~1 min: chain moves vs explicit Best Fit, a few certificates checked end to end
make verify     # verify every certificate with two independent exact verifiers; regenerate results/
make verify-py  # additionally a third (Python) verifier on certificates with <= 300k patterns
make simcheck   # simulate every certified pair and compare with its certificate (sanity check)
make figures    # figures/*.png (needs matplotlib)
```

Verify one certificate by hand:
```
C=certs/9_13_linear_vi_M1-2-1-2-2-2_R.cert.gz; M=1,2,1,2,2,2
./bfcert verify 9 13 $M $C reach
./bfcert states 9 13 $M lp/x && ./verify2 9 13 $M lp/x $C reach
```
Certificate format and the claim each mode certifies: `certs/README.md`.

## Layout

| path | what |
|---|---|
| `src/bfchain.hpp` | gap-count Markov chain, invariants, truncation patterns, abstract reachability |
| `src/bfcert.cpp` | exact verifier (`verify`), psi table (`states`), Theorem D containment check (`ycover`), self-test (`test`) |
| `src/verify2.cpp` | second exact verifier, independent of `bfchain.hpp` |
| `scripts/verify_indep.py` | third exact verifier (Python) |
| `scripts/verify_all.py` | verifies all certificates; writes `results/verification.csv` and derives `classification.csv`, `structure.csv`, `rates.csv` |
| `scripts/e2e_check.py` | checks a certificate on the states of a literal Best Fit run (no abstraction) |
| `scripts/make_tables.py`, `sim_check.py`, `make_figures.py` | `results/SUMMARY.md`; simulation cross-check; figures |
| `src/bfvi.cpp`, `scripts/kelley.py`, `run_vi_case.sh`, `classify.py`, `structure_all.py`, `rates_all.py` | certificate search (value iteration + cutting planes; needs numpy, highspy) |
| `src/bfsim.c`, `src/bfstats.c`, `scripts/boundary_sweep.py` | simulation (search guidance and the empirical boundary; never part of a proof) |

## License

MIT (see `LICENSE`).
