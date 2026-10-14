# Certificates

Each `*.cert.gz` file is a gzip-compressed Lyapunov certificate for Best Fit under U{j,k} (all three
verifiers read `.gz` directly). File name: `{j}_{k}_{mode}_vi_M{M_1-...-M_d0}_R.cert.gz` (`_R` =
patterns restricted to the abstract reachable closure, i.e. verify with the `reach` flag; `M` =
truncation levels of the small gaps 1..d0, where d0 = ⌊(k−1)/2⌋ or, when j < ⌊(k−1)/2⌋, possibly d0 = j).

Format (text):
```
mode <linear|bounded|bounded2|ylo|yhi>
D <common denominator>
c <k-1 integers>      numerators of c_1..c_{k-1}
psi
<integers ...>        numerators of psi over the compressed truncated-state index
```
The compressed index of psi is the enumeration order of `bfcert states` (`lp/<tag>/psi_states.txt` lists,
for each index, the truncated small-gap levels and the big gap); it depends only on (k, M).

| mode | claim when VERDICT = CERTIFIED | theorem |
|---|---|---|
| linear | min drift over all reachable patterns > 0 ⇒ E W_n ≥ γ n − O(1) | A |
| bounded | c_g > 0 (small g), max drift over saturated patterns < 0 ⇒ sup E e^{ηW_n} < ∞ | B |
| bounded2 | as bounded with gap 1 excluded (c_1 = 0) ⇒ waste = gap-1 bins + O(1) | C |
| ylo / yhi | c_1 = 1; min / max drift = lower / upper bound on the gap-1 creation rate (+ `bfcert ycover`) | D |

Verify all certificates and regenerate the tables: `make verify`. A single one:
```
./bfcert verify j k M_1,...,M_d0 certs/<file>.cert.gz reach          # exact C++ (__int128)
./bfcert states j k M_1,...,M_d0 lp/<tag>                              # writes lp/<tag>/psi_states.txt
./verify2       j k M_1,...,M_d0 lp/<tag> certs/<file>.cert.gz reach   # independent C++ (src/verify2.cpp)
python3 scripts/verify_indep.py j k M_1,...,M_d0 lp/<tag> certs/<file>.cert.gz reach   # independent
./bfcert ycover j k 0,M_2,...,M_d0                                     # Theorem D containment (ylo/yhi)
python3 scripts/e2e_check.py certs/<file>.cert.gz                      # literal Best Fit run
```
