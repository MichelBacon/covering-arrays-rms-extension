#!/usr/bin/env python3
"""
construct_layers.py - Explicit construction of CA(2q^3 + q^2*L ; 3, q^2+q+1+m, q)
==================================================================================
q a prime power (arithmetic in GF(q), see gf.py), m >= 2 new columns, L = max(q-2, ceil((q-2)*m / q)) layers of q^2 rows.

1. RMS base: CA(2q^3-1; 3, q^2+q+1, q) (Raaphorst-Moura-Stevens 2014).
2. For i = 1..m: column n_i = copy of column c_i = i-1, plus 1 in block 2.
   The only uncovered triples are {c_i, n_i, w} with n_i - c_i = e, e in {2..q-1}.
3. Layers: each layer = q^2 rows indexed by y in GF(q)^2.
   The q+1 directions of PG(1,q) give q+1 mutually orthogonal linear functions
   y -> y.P. In a layer, at most q copies are "active": active copy i gets its own
   direction P_i (c_i = y.P_i, n_i = y.P_i + e, sums in GF(q)); all other columns get direction
   P_0 (= y_1). Each copy is active in q-2 layers, with e running over the q-2 elements of GF(q) other than 0 and 1
   (labelled 2..q-1).
   -> for each copy i, each e, each other column w: (c_i, w) runs through all
      pairs (bijection), so all (a, e, b) are covered.
4. One correction row (0 on the original columns, 1 on the new ones):
   the zero vector is not in block 2, hence the combination (0,0,1) to add.
Usage: python construct_layers.py m [q]
"""
import sys
import math
import numpy as np
from rms_base import build
from gf import GF

m = int(sys.argv[1]); q = int(sys.argv[2]) if len(sys.argv) > 2 else 5
A, _ = build(q)
F = GF(q)
assert q >= 3, "the shifted-copy extension needs q >= 3"
k0 = A.shape[1]; n0 = 1 + (q**3 - 1)
assert 2 <= m <= k0
for i in range(m):
    nc = A[:, i].copy(); nc[n0:] = F.add[nc[n0:], 1]
    A = np.column_stack([A, nc])
k = A.shape[1]

act = [(i, e) for e in range(2, q) for i in range(m)]
if m <= q:                      # every layer contains all copies
    layers = [[(i, e) for i in range(m)] for e in range(2, q)]
else:                           # sliding windows of q activations (distinct copies since m > q)
    layers = [act[l * q:(l + 1) * q] for l in range(math.ceil(len(act) / q))]
L = len(layers)
dirs = [(1, 0)] + [(0, 1)] + [(1, j) for j in range(1, q)]       # P_0 ... P_q
Y = np.array([(y1, y2) for y1 in range(q) for y2 in range(q)])
extra = []
for lay in layers:
    assert len({i for i, _ in lay}) == len(lay)          # a copy at most once per layer
    B = np.tile(F.dot(Y, np.array([dirs[0]]))[:, 0], (k, 1)).T             # all columns: direction P_0
    for slot, (i, e) in enumerate(lay):
        f = F.dot(Y, np.array([dirs[slot + 1]]))[:, 0]
        B[:, i] = f
        B[:, k0 + i] = F.add[f, e]
    extra.append(B)
# correction row: the zero vector is missing from block 2 -> combination (0,0,1)
corr = np.zeros(k, dtype=int); corr[k0:] = 1
A = np.vstack([A] + extra + [corr[None, :]])
name = f"CA_{A.shape[0]}_3_{k}_{q}.txt"
np.savetxt(name, A, fmt="%d")
print(f"m={m} copies, {L} layer(s) -> CA({A.shape[0]}; 3, {k}, {q}) in {name}")
