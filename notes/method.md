# Best Fit under U{j,k}: certificates for bounded and linear expected waste

## Problem
Bins have capacity k, items are i.i.d. uniform on {1,...,j} (1 <= j <= k-1), packed online by Best Fit
(BF). Let W_n be the total empty space in the open bins after n items (waste). For j = k-1 the
expected waste is Theta(sqrt n); for j <= k-2 the optimal packing has O(1) expected waste, and the
question (Coffman, Johnson, Shor, Weber, STOC 1993) is for which (j,k) BF has E[W_n] = O(1) and for
which E[W_n] = Theta(n).

## Markov chain
State x = (x_1, ..., x_{k-1}), x_g = number of open bins with gap g. An item of size s goes to the
bin with the smallest gap g* >= s with x_{g*} > 0 (x_{g*}--, x_{g*-s}++ if g* > s); otherwise a new bin
is opened (x_{k-s}++). The waste is W(x) = sum_g g x_g.

**Invariant (CJSW93).** Any two open bins have gaps g, h with g + h < k. (A bin is opened with item s
only if every other gap is < s; afterwards gaps only shrink.) Hence with d0 = floor((k-1)/2) the
"small" gaps 1..d0 have unbounded counts while at most one bin has a "big" gap G > d0, and it
excludes small gaps h >= k - G.

**Lemma 0 (second invariant).** At most one open bin has a gap larger than j.
*Proof.* Let A be older than B. B was opened by an item s <= j that fitted into no open bin, so A's gap
was < s <= j at that moment, and gaps never increase. Hence of any two open bins the older has gap < j;
only the newest bin can have gap > j. []

Consequently, when j < (k-1)/2 we may take the "small" gaps to be 1..j (d0 = j) and treat every gap > j
as the single "big" bin; the code accepts either d0 = (k-1)/2 or d0 = j (selected by the number of
truncation levels). Both encodings over-approximate the reachable set, so either yields valid proofs.

## Certificate class
Fix truncation levels M_g >= 0 (small gaps g). Let m(x) = (min(x_g, M_g))_g together with the big
gap (0 if none), and

    V(x) = sum_{g <= d0} c_g x_g + psi(m(x)).

**Lemma 1 (finite description).** The one-step drift DV(x) = E[V(X_{t+1}) - V(X_t) | X_t = x] depends
on x only through its pattern p(x) = (min(x_g, M_g + 1))_g and the big gap.
*Proof.* The BF move for item s depends only on which counts are positive (a 1-limited chain), and
min(x_g + Delta_g, M_g) = min(p_g + Delta_g, M_g) for Delta_g in {-1, 0, 1}. []

So DV ranges over finitely many values, one per pattern; the patterns of reachable states are
contained in the closure of the empty pattern under abstract moves (a saturated level M_g + 1 that is
decremented may stay at M_g + 1 or drop to M_g), which we compute by BFS ("reach").

**Theorem A (linear waste).** If DV(x) >= gamma > 0 for every pattern in the reachable closure, then
for all n,   E[W_n] >= (gamma n - (max psi - psi(m(0)))) / C,   C = max_g max(c_g, 0)/g,
and liminf W_n / n >= gamma / C almost surely.
*Proof.* E V(X_n) = V(0) + sum_{t<n} E DV(X_t) >= V(0) + gamma n, and V(x) - V(0) <= C W(x) + max psi -
psi(m(0)). The a.s. statement follows from the SLLN for martingales with bounded increments applied
to V(X_n) - sum_{t<n} DV(X_t). []

**Theorem B (bounded waste).** If c_g > 0 for all small g and DV(x) <= -gamma < 0 for every reachable
pattern with at least one saturated coordinate (p_g = M_g + 1), then sup_n E[exp(eta W_n)] < infinity
for some eta > 0; in particular sup_n E[W_n] < infinity.
*Proof.* V has bounded increments. Non-saturated patterns are single states, so the exceptional
set U = {x : x_g <= M_g for all g} is finite; with a = max_U V the drift is <= -gamma on {V > a}.
Hajek (1982, Thm 2.3) gives sup_n E exp(eta V(X_n)) < infinity, and W(x) <= (k-1)(1 + sum_g x_g) <=
const (1 + V(x) - min psi) because c_g > 0. []

**Theorem C (structure of linear waste).** Suppose c_1 = 0, c_g > 0 for 2 <= g <= d0, and
DV(x) <= -gamma < 0 for every reachable pattern in which some gap g >= 2 is saturated ("bounded2"
certificate; psi may still depend on min(x_1, M_1)). Then sup_n E[exp(eta sum_{g>=2} x_g(n))] <
infinity; hence W_n = x_1(n) + O_P(1) and E[W_n] = E[x_1(n)] + O(1): all linear waste sits in bins
with gap exactly 1.
*Proof.* On patterns with no saturated gap g >= 2, V <= a := sum_{g>=2} c_g M_g + max psi (x_1
enters V only through psi). So the drift is <= -gamma on {V > a}; Hajek's theorem applies as in B. []

**Theorem D (waste rate).** Let Y be the decoupled chain in which gap-1 bins form an infinite reservoir
(size-1 items are no-ops; the coordinate y_1 counts gap-1 bins *created*), and a(y) the probability
that the next item creates a gap-1 bin. Fix truncation levels with M_1 = 0 and suppose
(i) X has a linear (A) and a bounded2 (C) certificate;
(ii) h_lo, h_hi are Y-certificates with c_1 = 1 (c_g, g >= 2, of any sign; psi does not see y_1) with
    rho_lo <= D h_lo(y)   and   D h_hi(y) <= rho_hi     on every pattern of the reachable closure of Y
(D h already contains the term c_1 a(y) = a(y));
(iii) every pattern of the reachable closure of X with p_1 = 1 (i.e. x_1 >= 1) lies in the reachable
closure of Y (`bfcert ycover`, checked for all rate certificates).
Then almost surely, and also for E[W_n]/n,
    rho_lo - 1/j  <=  liminf W_n/n  <=  limsup W_n/n  <=  rho_hi - 1/j.
*Proof.* By (i), liminf W_n/n > 0 a.s. (A, martingale SLLN) and S_n := sum_{g>=2} x_g(n) = o(n) a.s.
(C: sup_n E e^{eta S_n} < infinity, Markov's inequality and Borel-Cantelli). Since
W_n <= x_1(n) + (k-1)(1 + S_n), x_1(n) -> infinity, so a.s. x_1(n) = 0 only finitely often. Whenever
x_1 >= 1, the gaps >= 2 make exactly the move of Y (size-1 items fill a gap-1 bin; larger items never
use one), so with N_n = number of gap-1 bins created up to n and
H_n = N_n + sum_{g>=2} c_g x_g(n) + psi(.), the conditional increment E[H_{n+1} - H_n | F_n] equals
D h(y) evaluated at the Y-pattern of X_n, which by (iii) is in the closure where (ii) holds. H has
bounded increments, so H_n - sum_{t<n} E[dH_t | F_t] = o(n) a.s., and the finitely many steps with
x_1 = 0 change the sum by O(1). With H_n - N_n = O(1 + S_n) = o(n) this gives
rho_lo <= liminf N_n/n <= limsup N_n/n <= rho_hi a.s. Gap-1 bins only disappear when a size-1 item
fills one, which happens for every size-1 item except finitely many, so x_1(n) = N_n - n/j + o(n)
a.s. (SLLN). Finally W_n = x_1(n) + O(1 + S_n), and 0 <= W_n/n <= k-1 lets Fatou's lemma transfer the
bounds to E[W_n]/n. []
(No stationarity argument is needed, so the "ybounded" certificates computed earlier are not used.)

## Computing certificates: an adversarial MDP
For fixed c the best (psi, gamma) is a linear program with one constraint per pattern p:

    (1/j) sum_s [ c . Delta_s(p) + psi(q'_s(p)) ] - psi(q(p))  >= gamma        (linear case)

where q(p) = min(p, M) and q'_s(p) = min(p + Delta_s(p), M). Grouping constraints by q, this is the
Bellman inequality of an average-reward MDP on truncated states q in which an *adversary* chooses which
pattern p in pre(q) (i.e. which coordinates at level M_g are "really" larger) is realised. Hence

    gamma*(c) = min over adversary policies pi of  c . d_pi,

d_pi being the long-run mean increment vector of the gap counts under pi; gamma*(c) is concave in c,
and a subgradient is d_pi for the optimal pi. We maximise gamma*(c) (c_1 = 1, box) by Kelley's
cutting-plane method; each evaluation is a relative value iteration (bfvi), which also returns psi.
Any psi with min_q (T psi - psi)(q) > 0 is already a certificate (no convergence needed). The bounded
case is the max-adversary version restricted to saturated patterns (minimise over c >= 1).

**Proposition E (duality and monotonicity).** Fix M and c. The largest gamma for which some psi
satisfies the linear-mode inequalities equals the optimal (minimal) long-run average reward
min_pi c.d_pi of the finite adversarial MDP above (average-reward LP duality; the MDP is finite because
only reachable patterns enter). Hence a certificate with truncation M exists iff max_c min_pi c.d_pi > 0
(linear) or min_{c>0} max_pi c.d_pi < 0 (bounded). If M <= M' componentwise, every certificate at
truncation M is one at M' (psi'(q) = psi(min(q, M))), so the attainable margins are monotone in M.

The LP/MDP duality also explains failures: an adversary can "forget" increments of a saturated
coordinate (declare x_g = M_g when it is M_g + 1), which produces spurious growing drift vectors;
larger truncation levels on the fluctuating coordinates remove them.

## Verification
Certificates are (c, psi) with a common denominator D (integers). `bfcert verify` recomputes every
pattern's drift in exact integer arithmetic (__int128). `verify2` (C++) and `scripts/verify_indep.py`
(Python) are independent implementations (explicit Best Fit on concrete representative states, their
own reachability BFS, exact integers), sharing only the psi index table psi_states.txt. Sharing it is
harmless: the table is just a lookup that defines *some* function V of the state, and any V
satisfying the checked inequalities is a valid certificate; a missing or duplicated entry makes the
verifiers abort.
`scripts/e2e_check.py` is a further sanity check that does not use the pattern abstraction at all: it
runs literal Best Fit on a list of bins and, at every visited state, checks that the state's pattern
lies in the verified closure and that the exact one-step drift of V has the certified sign.
