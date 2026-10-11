#!/usr/bin/env python3
"""
exact_fibers.py - Minimum number of extra rows, fiber by fiber (search tool)
=============================================================================
All missing combinations involve the last column (the new one). An added row has
a single value x in that column, so the problem splits into v independent
subproblems (one per value x). For each, CP-SAT finds the smallest number of rows
covering the missing combinations of that fiber.
Usage: python exact_fibers.py seed.txt seconds_per_try output.txt
Requires: ortools
"""
import sys
import itertools
import numpy as np
from ortools.sat.python import cp_model

seed_f, tmax, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
q = 5
A = np.loadtxt(seed_f, dtype=int)
N0, k = A.shape
n = k - 1

missing = []
for cols in itertools.combinations(range(k), 3):
    codes = set((A[:, cols] @ [q * q, q, 1]).tolist())
    missing += [(cols, (c // 25, (c // 5) % 5, c % 5)) for c in range(q**3) if c not in codes]
assert all(n in cols for cols, _ in missing), "all missing combinations must involve the new column"
print(f"{len(missing)} missing combinations", flush=True)


def solve(fiber, R):
    """fiber: list of ((i,j),(a,b)) to cover with R rows."""
    cols = sorted({c for (ij, _) in fiber for c in ij})
    m = cp_model.CpModel()
    y = {(c, r, a): m.NewBoolVar("") for c in cols for r in range(R) for a in range(q)}
    for c in cols:
        for r in range(R):
            m.AddExactlyOne(y[c, r, a] for a in range(q))
    for (i, j), (a, b) in fiber:
        lits = []
        for r in range(R):
            z = m.NewBoolVar("")
            m.AddImplication(z, y[i, r, a]); m.AddImplication(z, y[j, r, b])
            lits.append(z)
        m.AddBoolOr(lits)
    c0 = cols[0]                                   # symmetry breaking
    val = [sum(a * y[c0, r, a] for a in range(q)) for r in range(R)]
    for r in range(R - 1):
        m.Add(val[r] <= val[r + 1])
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tmax; s.parameters.num_workers = 2
    st = s.Solve(m)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        rows = np.zeros((R, k), dtype=int)
        for c in cols:
            for r in range(R):
                rows[r, c] = next(a for a in range(q) if s.Value(y[c, r, a]))
        return rows, "found"
    return None, s.StatusName(st)


extra = []
for x in range(q):
    fiber = []
    for cols, vals in missing:
        if vals[cols.index(n)] != x:
            continue
        others = [(c, v) for c, v in zip(cols, vals) if c != n]
        fiber.append(((others[0][0], others[1][0]), (others[0][1], others[1][1])))
    bound = max(sum(1 for f in fiber if f[0] == ij) for ij in {f[0] for f in fiber})
    R = bound
    while True:
        rows, st = solve(fiber, R)
        print(f"  fiber x={x}: {len(fiber)} combinations, lower bound {bound}, try R={R} -> {st}", flush=True)
        if rows is not None:
            rows[:, n] = x
            extra.append(rows)
            break
        R += 1
E = np.vstack(extra)
np.savetxt(out, np.vstack([A, E]), fmt="%d")
print(f"wrote {out}: {N0} + {len(E)} = {N0 + len(E)} rows")
