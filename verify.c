/* verify.c - Independent strength-3 covering array checker (column-major storage, fast)
 * Build: gcc -O3 -march=native -o verify verify.c
 * Usage: ./verify file.txt v
 * Checks that every triple of columns covers all v^3 value combinations. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXN 25000
#define MAXK 512
static unsigned char T[MAXK][MAXN];
int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: verify file v\n"); return 1; }
    FILE *f = fopen(argv[1], "r"); int v = atoi(argv[2]); int N = 0, k = 0;
    static char line[8192];
    while (fgets(line, sizeof line, f)) {
        int j = 0; char *p = strtok(line, " \n");
        while (p) { T[j++][N] = atoi(p); p = strtok(NULL, " \n"); }
        if (j) { k = j; N++; }
    }
    int V3 = v * v * v; unsigned int *seen = calloc(V3, 4); unsigned int stamp = 0; long miss = 0;
    static unsigned short xy[MAXN];
    for (int x = 0; x < k; x++) for (int y = x + 1; y < k; y++) {
        for (int r = 0; r < N; r++) xy[r] = (T[x][r] * v + T[y][r]) * v;
        for (int z = y + 1; z < k; z++) {
            stamp++; int cnt = 0; const unsigned char *Z = T[z];
            for (int r = 0; r < N; r++) { int id = xy[r] + Z[r]; if (seen[id] != stamp) { seen[id] = stamp; cnt++; } }
            miss += V3 - cnt;
        }
    }
    printf("%s: N=%d k=%d v=%d -> %s (%ld missing)\n", argv[1], N, k, v, miss ? "INVALID" : "VALID", miss);
    return 0;
}
