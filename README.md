# New upper bounds for strength-3 covering arrays

*October 2026. The reference bounds are the Colbourn covering array tables (November 2024 status), as preserved by U. Grömping in the GitHub repository `ugroempi/CAs`. The original tables are no longer maintained since 2025.*

## The problem

A covering array CA(N; t, k, v) is an N × k array over v symbols in which every N × t subarray contains every t-tuple of symbols at least once. In software testing, the rows are test cases: k parameters with v values each, such that every combination of t parameter values is tested. CAN(t, k, v), the smallest possible N, is unknown for almost all parameters, so the community keeps tables of the best known upper bounds.

## Results: 60 improved entries

| v | k | New N | Previous bound (Nov 2024) | Gain |
|---|---|---|---|---|
| 5 | 32 | **324** | 349 | −25 |
| 5 | 33–36 | **325** | 349 | −24 |
| 5 | 37 | **350** | 365 | −15 |
| 7 | 58 | **930** | 979 | −49 |
| 7 | 59–64 | **931** | 1014 | −83 |
| 7 | 65 | **980** | 1015 | −35 |
| 11 | 134 | **3750** | 3871 | −121 |
| 11 | 135–144 | **3751** | 3971 | −220 |
| 11 | 145 | **3872** | 3971 | −99 |
| 13 | 184 | **6252** | 6421 | −169 |
| 13 | 185–196 | **6253** | 6565 | −312 |
| 13 | 197 | **6422** | 6565 | −143 |
| 19 | 382 | **19854** | 20215 | −361 |
| 19 | 383–400 | **19855** | 20539 to 20575 | −684 to −720 |
| 19 | 401 | **20216** | 20575 | −359 |

The full list is in `improvements.csv`. There is no improvement for q = 3, 17 or 23.

## Construction

For prime q and m ≥ 2 added columns:

> **CA( 2q³ + q²·L ; 3, q² + q + 1 + m, q )**, with L = max(q − 2, ⌈m(q − 2)/q⌉)

For m = 1: CA(2q³ − 1 + q²(q − 2) ; 3, q² + q + 2, q).

1. **Base (Raaphorst–Moura–Stevens, 2014)**: CA(2q³ − 1 ; 3, q² + q + 1, q).
   - Columns are the points p_j of PG(2, q), ordered along a Singer cycle (powers of a root α of a primitive cubic).
   - Block B1 has rows x ≠ 0 in GF(q)³ with entries ⟨x, p_j⟩.
   - Block B2 is the same with the column order reversed.
   - There is one all-zero row.
   - By the RMS theorem, three columns collinear in B1 are not collinear in B2.
2. **Shifted copies**: new column n_i is a copy of column c_i = i − 1 on the zero row and on B1, and equals **c_i + 1** on B2.
3. **Layers**: L layers of q² rows indexed by y ∈ GF(q)².
   - The q + 1 directions P_0, …, P_q of PG(1, q) give mutually orthogonal linear functions y ↦ y·P_s.
   - In each layer, at most q copies are *active*: copy i gets its own direction (c_i = y·P_s, n_i = y·P_s + e). All other columns get P_0.
   - Each copy is active in q − 2 layers, with e = 2, …, q − 1.
4. **One correction row**: 0 on the original columns, 1 on the new columns.

### Proof sketch

- **Three original columns**: by the RMS theorem.
- **Triple {a, b, n_i} with c_i ∉ {a, b}**:
  - Replacing each new column by its source column gives three distinct original columns, which are non-collinear in B1 or in B2.
  - In that block, x ↦ (values) is a bijection of GF(q)³ followed by a translation, so every 3-tuple is covered except the image of x = 0.
  - In B1 that image is (0,0,0), covered by the zero row.
  - In B2 it has 1s on new columns and 0s on original columns, which is covered by the correction row.
- **Triple {c_i, n_i, w}**, with ε = n_i − c_i:
  - On B1 and the zero row, ε = 0 and (c_i, w) runs through all pairs (two independent linear forms).
  - On B2, ε = 1, by the same argument.
  - In the layers where copy i is active with value e, ε = e and (c_i, w) = (y·P_s, f_w(y)), with f_w along a different direction. This is a bijection of GF(q)², so every pair appears.
  - Hence ε = 0, 1, 2, …, q − 1 all occur with every pair (c_i, w), and the triple is fully covered. ∎

## Verification

Every array was checked by brute force (every triple of columns, every value combination) with checkers written independently of the construction code: `verify.py` (Python) and `verify.c` (C). Sanity check of the checkers: removing one row of CA(324;3,32,5), or changing a single cell, is correctly reported as INVALID.

All arrays in `arrays/` are verified. The q = 19 arrays (about 20 MB each) were also verified but are not included; `python construct_layers.py 20 19` regenerates them.

## Files

- `rms_base.py`: rebuilds the RMS base array.
- `construct_one_copy.py`: construction with one added column (m = 1).
- `construct_layers.py`: general construction (m ≥ 2).
- `verify.py`, `verify.c`: independent checkers.
- `arrays/CA_N_t_k_v.txt`: the arrays, one test per line, values 0 … v − 1.
- `improvements.csv`: the 60 improved entries.
- `search/`: exploration tools used along the way (simulated annealing in C, CP-SAT model).
- `R-contribution/`: an implementation for the R package `CAs` (functions `RMScopies`, `RMScopiesCA`, `N_RMScopiesCA`, with documentation and tests), plus the corresponding git patch.

## Caveats

- Novelty is established **relative to the November 2024 tables** and the published constructions I could find (Shokri, Moura and Stevens 2023–2025 start at k = 2(q² + q + 1)). A more recent unindexed result is possible.
- The code handles prime q only. Prime powers (4, 8, 9, 16, 25) would need GF(q) arithmetic.

## Acknowledgement

Found and developed with the help of Claude (Anthropic), an AI assistant.
