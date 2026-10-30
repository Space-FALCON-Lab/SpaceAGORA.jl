"""Self-test of the oracle DP and policies on synthetic drifts with closed-form optima, plus brute force.
Run:  python3 -I oracle_selftest.py     (no data needed)

Closed forms (band [-w, w] after shrinking, constant lever L unless stated, burn opportunity at every pass):
 1. constant drift d per pass over K passes starting inside the band:  dV* = L * max(0, x_K - w)
    (the offset must absorb everything past the upper edge); minimum burn count = ceil((x_K - w) / (2w + d)).
 2. sinusoid of amplitude A > w with N extrema (maxima/minima alternate) in the record:
    dV* = L * (A - w) * (2 N - 1)   (first extremum needs A - w, each later reversal 2 (A - w)).
 3. step of height S > w at a pass k0, lever L1 before and L2 after (pre-positioning allowed):
    dV* = min(L1, L2) pre-positions up to w;  L1 < L2:  L1 * min(w, S - w) + L2 * max(0, S - 2w);  L1 >= L2: L2 * (S - w);
    S <= w: 0.
 4. brute-force enumeration over all grid paths for random small instances (with burn penalty and min burn dV).
"""
import importlib.util, itertools, math, pathlib, sys
import numpy as np

here = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("oracle_lib", here / "oracle_lib.py")
lib = importlib.util.module_from_spec(spec); sys.modules["oracle_lib"] = lib; spec.loader.exec_module(lib)

fails = 0
def check(name, got, want, tol=1e-9):
    global fails
    ok = abs(got - want) <= tol * max(1.0, abs(want))
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got {got:.9g}  want {want:.9g}")
    fails += (not ok)

def run_all(x, w, L, res):
    n = len(x); lo = np.full(n, -w, float); hi = np.full(n, w, float); Lv = np.full(n, L) if np.isscalar(L) else np.asarray(L)
    d = lib.dp_schedule(lo - x, hi - x, Lv, res)
    return d, lib.lazy(x, lo, hi, Lv), lib.near_passive(x, lo, hi, Lv), lib.reset_opposite(x, lo, hi, Lv), lo, hi, Lv

L = 7.0e-5
# 1. constant drift
d_, K, w = 300.0, 100, 4000.0
x = d_ * np.arange(1, K + 1)
dp, lz, *_ = run_all(x, w, L, 100.0)
want = L * (x[-1] - w)
check("1 const drift: DP", dp["dv"], want); check("1 const drift: lazy", lz[0], want)
lo = np.full(K, -w); hi = np.full(K, w); Lv = np.full(K, L)
need = math.ceil((x[-1] - w) / (2 * w + d_))
big = lib.dp_schedule(lo - x, hi - x, Lv, 100.0, lam=10.0)
check("1 const drift: min burn count at large lambda", big["burns"], need)
check("1 const drift: dV at min burns", big["dv"], want)
check("1 const drift: greedy hold burns", lib.hold(x, lo, hi, Lv, "nearest")[1], need)
# 2. sinusoid
A, w, P, K = 10000.0, 6000.0, 40, 131
k = np.arange(1, K + 1); x = np.round(A * np.sin(2 * np.pi * k / P))
N = int(sum(1 for j in range(1, K) if abs(x[j - 1]) == A))   # extrema present at integer sample points
dp, lz, *_ = run_all(x, w, L, 1.0)
want = L * (A - w) * (2 * N - 1)
check(f"2 sinusoid (N={N} extrema): DP", dp["dv"], want); check("2 sinusoid: lazy", lz[0], want)
# 3. step
for L1, L2, S, w in ((1.0e-4, 2.0e-4, 9000.0, 6000.0), (1.0e-4, 2.0e-4, 15000.0, 6000.0), (2.0e-4, 1.0e-4, 15000.0, 6000.0), (1e-4, 2e-4, 5000.0, 6000.0)):
    K, k0 = 60, 30
    x = np.where(np.arange(K) < k0, 0.0, S)
    Lv = np.where(np.arange(K) < k0, L1, L2)
    dp, lz, *_ = run_all(x, w, Lv, 100.0)
    want = 0.0 if S <= w else (L1 * min(w, S - w) + L2 * max(0.0, S - 2 * w) if L1 < L2 else L2 * (S - w))
    check(f"3 step S={S:.0f} L1={L1:g} L2={L2:g}: DP", dp["dv"], want)
    lazy_want = L2 * max(0.0, S - w)
    check(f"3 step: lazy (reactive, no pre-position)", lz[0], lazy_want)
# 4. brute force
rng = np.random.default_rng(12345)
for trial in range(40):
    n = 5; res = 1.0
    Ak = rng.integers(-4, 2, n).astype(float); Bk = Ak + rng.integers(1, 5, n)
    Lv = rng.uniform(0.5, 2.0, n); lam = rng.choice([0.0, 0.7, 3.0]); dvmin = rng.choice([0.0, 1.0, 2.5])
    d = lib.dp_schedule(Ak, Bk, Lv, res, lam=lam, dvmin=dvmin)
    grid = range(-8, 9); best = math.inf
    for path in itertools.product(*[[g for g in grid if Ak[j] <= g <= Bk[j]] for j in range(n)]):
        prev = 0; c = 0.0; ok = True
        for j, g in enumerate(path):
            dd = g - prev
            if dd != 0:
                cost = Lv[j] * abs(dd)
                if dvmin and cost < dvmin - 1e-12: ok = False; break
                c += cost + lam
            prev = g
        if ok: best = min(best, c)
    got = math.inf if d is None else d["obj"]
    if not (math.isinf(best) and d is None):
        check(f"4 brute force trial {trial} lam={lam} dvmin={dvmin}", got, best, 1e-9)
# 5. exact (candidate-offset) mode: equals the grid mode on integer data, <= fine grid on continuous data, and matches an LP (scipy, if present)
for trial in range(30):
    n = int(rng.integers(4, 9)); Ak = rng.uniform(-5, 1, n); Bk = Ak + rng.uniform(0.5, 4, n); Lv = rng.uniform(0.5, 2.0, n)
    lam = float(rng.choice([0.0, 0.7]))
    ex = lib.dp_schedule(Ak, Bk, Lv, None, lam=lam); gr = lib.dp_schedule(Ak, Bk, Lv, 1e-3, lam=lam)
    if ex is None: continue
    ok = ex["obj"] <= gr["obj"] + 1e-9 and gr["obj"] - ex["obj"] <= 2e-3 * (Lv.sum() + 1) * (1 + (lam > 0))
    print(f"{'PASS' if ok else 'FAIL'}  5 exact vs grid(1e-3) trial {trial}: exact {ex['obj']:.6f} grid {gr['obj']:.6f}"); fails += (not ok)
    if lam == 0.0:
        try:
            from scipy.optimize import linprog
            # variables: u_1..u_n, t_1..t_n ; min sum L_k t_k ; t_k >= +-(u_k - u_{k-1})
            c = np.concatenate((np.zeros(n), Lv)); Aub = []; bub = []
            for k in range(n):
                for sg in (1, -1):
                    row = np.zeros(2 * n); row[k] = sg; row[n + k] = -1
                    if k > 0: row[k - 1] = -sg
                    Aub.append(row); bub.append(0.0)
            r = linprog(c, A_ub=np.array(Aub), b_ub=bub, bounds=[(Ak[k], Bk[k]) for k in range(n)] + [(0, None)] * n, method="highs")
            check(f"5 exact vs LP trial {trial}", ex["obj"], r.fun, 1e-7)
        except ImportError:
            pass
# 6. burn-budget DP: best[c] equals the penalty-DP solution for c = its burn count; best[inf] equals the unconstrained optimum
for trial in range(20):
    n = int(rng.integers(4, 9)); Ak = rng.uniform(-5, 1, n); Bk = Ak + rng.uniform(0.5, 4, n); Lv = rng.uniform(0.5, 2.0, n)
    bud = lib.dp_budget(Ak, Bk, Lv, n)
    ex = lib.dp_schedule(Ak, Bk, Lv, None)
    if ex is None: continue
    check(f"6 budget DP unconstrained trial {trial}", bud[-1], ex["dv"], 1e-9)
    for lam in (0.3, 1.5):
        e2 = lib.dp_schedule(Ak, Bk, Lv, None, lam=lam)
        ok = bud[e2["burns"]] <= e2["dv"] + 1e-9
        print(f"{'PASS' if ok else 'FAIL'}  6 budget DP <= penalty DP at its burn count (lam={lam}) trial {trial}"); fails += (not ok)
print("ALL PASS" if fails == 0 else f"{fails} FAILURES")
sys.exit(1 if fails else 0)
