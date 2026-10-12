# New upper bounds for strength-3 covering arrays

*October 2026. The reference bounds are the Colbourn covering array tables (November 2024 status), as preserved by U. Grömping in the GitHub repository `ugroempi/CAs`. The original tables have not been maintained since 2025.*

## The problem

A covering array CA(N; t, k, v) is an N × k array over v symbols in which every N × t subarray contains every t-tuple of symbols at least once. In software testing, the rows are test cases: k parameters with v values each, such that every combination of t parameter values is tested. CAN(t, k, v), the smallest possible N, is unknown for almost all parameters, so the community keeps tables of the best known upper bounds.

## Results: 78 improved entries

| v | k | New N | Previous bound (Nov 2024) | Gain |
|---|---|---|---|---|
| 5 | 32 | **324** | 349 | −25 |
| 5 | 33–36 | **325** | 349 | −24 |
| 5 | 37 | **350** | 365 | −15 |
| 7 | 58 | **930** | 979 | −49 |
| 7 | 59–64 | **931** | 1014 | −83 |
| 7 | 65 | **980** | 1015 | −35 |
| 8 | 74 | **1407** | 1447 | −40 |
| 8 | 75 | **1408** | 1447 | −39 |
| 8 | 76–81 | **1408** | 1455 | −47 |
| 9 | 92 | **2024** | 2105 | −81 |
| 9 | 93–100 | **2025** | 2169 | −144 |
| 9 | 101 | **2106** | 2169 | −63 |
| 11 | 134 | **3750** | 3871 | −121 |
| 11 | 135–144 | **3751** | 3971 | −220 |
| 11 | 145 | **3872** | 3971 | −99 |
| 13 | 184 | **6252** | 6421 | −169 |
| 13 | 185–196 | **6253** | 6565 | −312 |
| 13 | 197 | **6422** | 6565 | −143 |
| 19 | 382 | **19854** | 20215 | −361 |
| 19 | 383–400 | **19855** | 20539 to 20575 | −684 to −720 |
| 19 | 401 | **20216** | 20575 | −359 |

The full list is in `improvements.csv`. The construction works for every prime power q ≥ 3, but gives no improvement for q = 3, 4, 16, 17, 23 or 25 (the tables stop at v = 25).

## Construction

For a prime power q ≥ 3 and m ≥ 2 added columns (m ≤ q² + q + 1):

> **CA( 2q³ + q²·L ; 3, q² + q + 1 + m, q )**, with L = max(q − 2, ⌈m(q − 2)/q⌉)

For m = 1: CA(2q³ − 1 + q²(q − 2) ; 3, q² + q + 2, q).

All arithmetic is in GF(q). The elements other than 0 and 1 are denoted e₂, …, e_{q−1} (in the code: the integers 2, …, q − 1).

1. **Base (Raaphorst–Moura–Stevens, 2014)**: CA(2q³ − 1 ; 3, q² + q + 1, q).
   - Columns j = 0, …, k₀ − 1 (k₀ = q² + q + 1) are the points p_j = αʲ of PG(2, q), where α is a root of a primitive cubic over GF(q) (Singer cycle).
   - Block B1 has one row per x ≠ 0 in GF(q)³, with entries ⟨x, p_j⟩.
   - Block B2 is the same with the column order reversed (column j uses p_{k₀−1−j}).
   - There is one all-zero row.
   - RMS theorem: three columns whose points are collinear in B1 have non-collinear points in B2.
2. **Shifted copies**: new column n_i (i = 1, …, m) is a copy of column c_i = i − 1 on the zero row and on B1, and equals **c_i + 1** on B2.
3. **Layers** (m ≥ 2): L layers of q² rows, indexed by y ∈ GF(q)².
   - The q + 1 directions P₀ = (1,0), P₁ = (0,1), P_{1+j} = (1, a_j) (a_j ≠ 0) of PG(1, q) give linear forms y ↦ y·P_s, any two of which form a bijection GF(q)² → GF(q)².
   - The m(q − 2) "activations" (i, e), with e ∈ {e₂, …, e_{q−1}}, are distributed over the layers so that each layer contains at most q activations, with distinct copies. If m ≤ q, layer ℓ activates every copy with e = e_{ℓ+1}. If m > q, the list of activations, ordered by e then i, is cut into consecutive windows of q; the copies in a window are distinct because m > q.
   - In a layer, the s-th active copy (s = 1, …, at most q) with value e gets c_i = y·P_s and n_i = y·P_s + e. Every other column (original or new) gets y·P₀.
4. **One correction row** (m ≥ 2): 0 on the original columns, 1 on the new columns.
5. **For m = 1**, the layers and the correction row are replaced by the q²(q − 2) rows (c₁ = a, n₁ = a + e, every other column = s), for a, s ∈ GF(q) and e ∈ {e₂, …, e_{q−1}}.

## Proof

Fix three distinct columns T. We show that every 3-tuple of values appears on T.

**Case 1: T has three original columns.** By the RMS theorem, the three points are non-collinear in B1 or in B2. In that block, x ↦ (⟨x, p_a⟩, ⟨x, p_b⟩, ⟨x, p_c⟩) is a bijection of GF(q)³. The non-zero x give every tuple except (0, 0, 0), which is on the zero row.

**Case 2: T contains new columns, but no pair {c_i, n_i}.**
- Replace each new column n_i by its source c_i. This gives three distinct original columns (two new columns have distinct sources), and therefore three points that are non-collinear in B1 or in B2.
- On B1, the values on T are exactly those of the sources. On B2, they are those of the sources plus a constant vector δ (1 on the new columns, 0 on the others).
- So on the good block, x ↦ (values on T) is a bijection of GF(q)³ followed by a translation. The non-zero x give every tuple except the image of x = 0:
  - on B1, the image of 0 is (0, 0, 0), which is on the zero row (new columns copy their sources there);
  - on B2, it is δ, which is on the correction row for m ≥ 2.
  - For m = 1, δ has a single 1 on n₁, and c₁ ∉ T. Choose a ∉ {0, 1} and e = 1 − a (then e ∉ {0, 1}). The added row with c₁ = a, n₁ = 1 and s = 0 covers δ.

**Case 3: T = {c_i, n_i, w}.** Write ε = n_i − c_i. It suffices to show that, for each ε ∈ GF(q), the pair (c_i, w) takes all q² values on the rows where n_i − c_i = ε.
- **ε = 0, on B1.** The sources of c_i and w are distinct points. So x ↦ (⟨x, p⟩, ⟨x, p′⟩) is a linear surjection GF(q)³ → GF(q)², with a one-dimensional kernel. Every pair is the image of q vectors x, and at least q − 1 ≥ 1 of them are non-zero. This also holds if w is another new column n_j.
- **ε = 1, on B2.** Same argument. The constant shift of n_j does not change surjectivity.
- **ε = e ∉ {0, 1}, m ≥ 2.** Take the layer where (i, e) is active, with direction P_s. Then w carries a different direction P_{s′}, plus possibly a constant: P₀ if w is not active in that layer, or the direction of its own activation otherwise. Hence y ↦ (c_i, w) is a bijection of GF(q)² followed by a translation, so every pair appears.
- **ε = e ∉ {0, 1}, m = 1.** The rows (c₁ = a, n₁ = a + e, w = s) give all pairs (a, s).

∎

**Run count.** The base has 1 + 2(q³ − 1) rows. The layers add q²L rows, with L = max(q − 2, ⌈m(q − 2)/q⌉), and the correction row adds 1, for a total of 2q³ + q²L. For m = 1, the added rows are q²(q − 2).

## Verification

Every array in `arrays/` was checked by brute force (every triple of columns, every value combination) with checkers written independently of the construction code: `verify.py` (Python) and `verify.c` (C).
- Sanity check of the checkers: removing one row of CA(1407;3,74,8), or changing a single cell of CA(2025;3,100,9), is correctly reported as INVALID.
- The arrays for q = 19 (about 20 MB each) and a CA(11776; 3, 289, 16) were also verified but are not included. Run `python construct_layers.py 20 19` to regenerate them.
- The R implementation in `R-contribution/` produces the same arrays, up to row order, for q = 4, 5, 7, 8, 9.
- For prime q, the GF(q) version of the code produces arrays byte-identical to the earlier prime-only version.

## Files

- `gf.py`: GF(q) arithmetic for prime powers. The field axioms are checked exhaustively when it is built.
- `rms_base.py`: rebuilds the RMS base array.
- `construct_one_copy.py`: construction with one added column (m = 1).
- `construct_layers.py`: general construction (m ≥ 2).
- `verify.py`, `verify.c`: independent checkers.
- `arrays/CA_N_t_k_v.txt`: the arrays, one test per line, values 0 … v − 1.
- `improvements.csv`: the 78 improved entries.
- `search/`: exploration tools used along the way (simulated annealing in C, CP-SAT model).
- `R-contribution/`: an implementation for the R package `CAs`, with the corresponding git patch:
  - functions `RMScopies`, `RMScopiesCA` and `N_RMScopiesCA`, for prime powers via `lhs::create_galois_field`;
  - documentation and tests.

## Caveats

- Novelty is established **relative to the November 2024 tables** and to the published constructions I could find (Shokri, Moura and Stevens 2023–2025 start at k = 2(q² + q + 1)). A more recent unindexed result is possible, and the idea may be considered folklore by specialists.
- The proof above has not yet been checked by an independent human expert.

## Acknowledgement

Found and developed with the help of Claude (Anthropic), an AI assistant.
