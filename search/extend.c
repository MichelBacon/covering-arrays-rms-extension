/*
 * extend.c - Extends a strength-3 covering array by adding columns (search tool)
 * ==============================================================================
 * Input: a valid CA(N0; 3, k0, v) (text file), target k > k0 columns,
 *        R initial extra rows.
 * Free cells: new columns in the base rows + all cells of the added rows.
 * Base cells of the old columns are frozen (they already cover all triples of
 * old columns). Targeted simulated annealing: pick an uncovered combination and
 * modify a compatible row to cover it. When everything is covered, remove one
 * added row and start again -> minimizes the number of rows.
 *
 * Build: gcc -O3 -march=native -o extend extend.c -lm
 * Usage: ./extend base.txt v k R seconds seed output.txt [kfixed] [nfrozen]
 *   kfixed : columns < kfixed are frozen in the base rows (default: base width)
 *   nfrozen: number of frozen rows; further rows in the file are free (warm start)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAXN 1200
#define MAXK 64

static int N0, k0, v, k, N;
static unsigned char A[MAXN][MAXK];
static unsigned char FREE[MAXN][MAXK];
static int *tri;                  /* tri[(a*k+b)*k+c] -> index of triple a<b<c */
static int ntri, V3;
static int *cnt;                  /* cnt[t*V3 + code] */
static int *upos;                 /* position in the uncovered list, -1 otherwise */
static int *ulist, nun;
static unsigned long long rs;

static inline unsigned long long rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double urand(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }

static inline void add_un(int id) { upos[id] = nun; ulist[nun++] = id; }
static inline void rem_un(int id) { int p = upos[id], last = ulist[--nun]; ulist[p] = last; upos[last] = p; upos[id] = -1; }

static inline int code3(int r, int a, int b, int c) { return (A[r][a] * v + A[r][b]) * v + A[r][c]; }

/* set cell (r,c) to value nv; returns the change in the number of uncovered combinations */
static int setcell(int r, int c, int nv) {
    int old = A[r][c], d = 0;
    if (old == nv) return 0;
    for (int a = 0; a < k; a++) {
        if (a == c) continue;
        for (int b = a + 1; b < k; b++) {
            if (b == c) continue;
            int x = a, y = b, z = c;            /* sort (a,b,c) */
            if (z < y) { int t = y; y = z; z = t; }
            if (y < x) { int t = x; x = y; y = t; }
            int t = tri[(x * k + y) * k + z];
            int base = t * V3;
            int id = base + code3(r, x, y, z);
            if (--cnt[id] == 0) { add_un(id); d++; }
            A[r][c] = nv;
            id = base + code3(r, x, y, z);
            if (cnt[id]++ == 0) { rem_un(id); d--; }
            A[r][c] = old;
        }
    }
    A[r][c] = nv;
    return d;
}

static void rebuild(void) {
    memset(cnt, 0, sizeof(int) * (size_t)ntri * V3);
    for (int a = 0; a < k; a++) for (int b = a + 1; b < k; b++) for (int c = b + 1; c < k; c++) {
        int t = tri[(a * k + b) * k + c];
        for (int r = 0; r < N; r++) cnt[t * V3 + code3(r, a, b, c)]++;
    }
    nun = 0;
    for (int i = 0; i < ntri * V3; i++) { upos[i] = -1; if (cnt[i] == 0) add_un(i); }
}

static int tcols[MAXK * MAXK * MAXK / 6 + 10][3];

int main(int argc, char **argv) {
    if (argc < 8) { fprintf(stderr, "usage: base v k R seconds seed output [kfixed] [nfrozen]\n"); return 1; }
    FILE *f = fopen(argv[1], "r");
    v = atoi(argv[2]); k = atoi(argv[3]); int R = atoi(argv[4]);
    double tmax = atof(argv[5]); rs = 88172645463325252ULL ^ (unsigned long long)atoll(argv[6]) * 2654435761ULL;
    const char *out = argv[7];
    int kfixed = argc > 8 ? atoi(argv[8]) : -1;
    int nfrozen = argc > 9 ? atoi(argv[9]) : -1;   /* frozen rows; further rows of the file are free */
    /* read the base array */
    char line[8192]; N0 = 0; k0 = 0;
    while (fgets(line, sizeof line, f)) {
        int c = 0; char *p = strtok(line, " \t\n");
        while (p) { A[N0][c++] = (unsigned char)atoi(p); p = strtok(NULL, " \t\n"); }
        if (c) { k0 = c; N0++; }
    }
    fclose(f);
    V3 = v * v * v; int Nfile = N0; N = Nfile + R;
    if (kfixed < 0) kfixed = k0;
    if (nfrozen >= 0) N0 = nfrozen;                    /* N0 = number of frozen rows */
    for (int r = 0; r < N; r++) for (int c = 0; c < k; c++) {
        FREE[r][c] = (r >= N0 || c >= kfixed);
        if (FREE[r][c] && !(r < Nfile && c < k0)) A[r][c] = rnd() % v;   /* keep the values read from the file */
    }
    tri = malloc(sizeof(int) * k * k * k); ntri = 0;
    for (int a = 0; a < k; a++) for (int b = a + 1; b < k; b++) for (int c = b + 1; c < k; c++) {
        tri[(a * k + b) * k + c] = ntri; tcols[ntri][0] = a; tcols[ntri][1] = b; tcols[ntri][2] = c; ntri++;
    }
    cnt = malloc(sizeof(int) * (size_t)ntri * V3);
    upos = malloc(sizeof(int) * (size_t)ntri * V3);
    ulist = malloc(sizeof(int) * (size_t)ntri * V3);
    rebuild();
    fprintf(stderr, "base %dx%d, target k=%d, start N=%d, uncovered=%d\n", N0, k0, k, N, nun);

    clock_t t0 = clock();
    double T = 1.0;
    int bestN = -1;
    long long it = 0;
    int ch_r[3], ch_c[3], ch_old[3];
    while ((double)(clock() - t0) / CLOCKS_PER_SEC < tmax) {
        if (nun == 0) {
            /* save, then remove the least useful added row */
            bestN = N;
            FILE *o = fopen(out, "w");
            for (int r = 0; r < N; r++) { for (int c = 0; c < k; c++) fprintf(o, c ? " %d" : "%d", A[r][c]); fprintf(o, "\n"); }
            fclose(o);
            fprintf(stderr, "  [%.0fs] VALID with N=%d (k=%d)\n", (double)(clock() - t0) / CLOCKS_PER_SEC, N, k);
            if (N <= N0) break;
            /* remove the chosen added row: overwrite it with row N-1 */
            int bestr = -1, bestloss = 1 << 30;
            for (int r = N0; r < N; r++) {
                int loss = 0;
                for (int t = 0; t < ntri; t++) {
                    int id = t * V3 + code3(r, tcols[t][0], tcols[t][1], tcols[t][2]);
                    if (cnt[id] == 1) loss++;
                }
                if (loss < bestloss) { bestloss = loss; bestr = r; }
            }
            memcpy(A[bestr], A[N - 1], k);
            N--;
            rebuild();
            T = 1.0;
            continue;
        }
        it++;
        /* targeted move */
        int id = ulist[rnd() % nun];
        int t = id / V3, cd = id % V3;
        int vals[3] = { cd / (v * v), (cd / v) % v, cd % v };
        int *cs = tcols[t];
        int r = -1;
        for (int tries = 0; tries < 64; tries++) {
            int rr = rnd() % N, ok = 1, nfree = 0;
            for (int j = 0; j < 3; j++) {
                if (!FREE[rr][cs[j]]) { if (A[rr][cs[j]] != vals[j]) { ok = 0; break; } }
                else nfree++;
            }
            if (ok && nfree) { r = rr; break; }
        }
        if (r < 0) continue;
        int nc = 0, d = 0;
        for (int j = 0; j < 3; j++) if (FREE[r][cs[j]] && A[r][cs[j]] != vals[j]) {
            ch_r[nc] = r; ch_c[nc] = cs[j]; ch_old[nc] = A[r][cs[j]]; nc++;
            d += setcell(r, cs[j], vals[j]);
        }
        if (d <= 0 || urand() < exp(-d / T)) { /* accepted */ }
        else for (int j = nc - 1; j >= 0; j--) setcell(ch_r[j], ch_c[j], ch_old[j]);
        if ((it & 0xFFFF) == 0) {
            T *= 0.97; if (T < 0.05) T = 0.6;
        }
        if ((it & 0x3FFFF) == 0) {
            static double lastp = 0; double el = (double)(clock() - t0) / CLOCKS_PER_SEC;
            if (el - lastp > 5) { lastp = el; fprintf(stderr, "    [%.0fs] it=%lld N=%d uncovered=%d T=%.3f\n", el, it, N, nun, T); }
        }
    }
    printf("%d\n", bestN);
    return 0;
}
