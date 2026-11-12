/**
 * test_qmatinv.c, host-side accuracy test against a double-precision reference.
 *   gcc -O2 -std=c99 -Wall -Wextra qmatinv.c test_qmatinv.c -o test_qmatinv -lm
 * License: MIT
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "qmatinv.h"

static int inv_double(const double *A, double *X, int n)
{
    double M[QMI_MAX_N][2 * QMI_MAX_N];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < 2 * n; j++)
            M[i][j] = j < n ? A[i * n + j] : (j - n == i ? 1.0 : 0.0);
    for (int p = 0; p < n; p++) {
        int mr = p;
        for (int r = p + 1; r < n; r++) if (fabs(M[r][p]) > fabs(M[mr][p])) mr = r;
        if (M[mr][p] == 0.0) return -1;
        for (int k = 0; k < 2 * n; k++) { double t = M[p][k]; M[p][k] = M[mr][k]; M[mr][k] = t; }
        double pv = M[p][p];
        for (int k = 0; k < 2 * n; k++) M[p][k] /= pv;
        for (int r = 0; r < n; r++) if (r != p) {
            double f = M[r][p];
            for (int k = 0; k < 2 * n; k++) M[r][k] -= f * M[p][k];
        }
    }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) X[i * n + j] = M[i][n + j];
    return 0;
}

static double urand(double lo, double hi) { return lo + (hi - lo) * (rand() / (double)RAND_MAX); }

int main(void)
{
    /* 1) worked example */
    float A2[4] = { -6.3f, 2.0f, 3.0f, -4.0f };
    qmi_t X2[4];
    qmi_status_t st = qmi_inverse_f32in(A2, X2, 2);
    double Ad[4] = { -6.3, 2.0, 3.0, -4.0 }, Xd[4];
    inv_double(Ad, Xd, 2);
    printf("Example 2x2 [[-6.3 2][3 -4]] status=%d\n", st);
    for (int i = 0; i < 4; i++)
        printf("  inv[%d] fixed % .8f  double % .8f  err % .2e\n", i, qmi_to_float(X2[i]), Xd[i], qmi_to_float(X2[i]) - Xd[i]);

    /* 2) random sweep */
    srand(12345);
    const int TRIALS = 20000;
    printf("\nRandom sweep: entries U(-10,10), %d matrices per size, Q%d.%d\n",
           TRIALS, 31 - QMI_FRAC, QMI_FRAC);
    printf("  n | ok     | singular | saturated | max abs err | median rel err | 99th pct rel err\n");
    for (int n = 2; n <= 6; n++) {
        int ok = 0, sing = 0, sat = 0;
        static double rel[20000];
        double maxabs = 0.0;
        for (int t = 0; t < TRIALS; t++) {
            float A[QMI_MAX_N * QMI_MAX_N]; double Ad2[QMI_MAX_N * QMI_MAX_N], Xr[QMI_MAX_N * QMI_MAX_N];
            for (int k = 0; k < n * n; k++) { A[k] = (float)urand(-10, 10); Ad2[k] = A[k]; }
            qmi_t X[QMI_MAX_N * QMI_MAX_N];
            st = qmi_inverse_f32in(A, X, n);
            if (st == QMI_SINGULAR) { sing++; continue; }
            if (st == QMI_SATURATED) { sat++; continue; }
            if (inv_double(Ad2, Xr, n)) { sing++; continue; }
            double num = 0, den = 0;
            for (int k = 0; k < n * n; k++) {
                double e = fabs((double)X[k] / (double)(1L << QMI_FRAC) - Xr[k]);
                if (e > maxabs) maxabs = e;
                num += e * e; den += Xr[k] * Xr[k];
            }
            rel[ok++] = sqrt(num / den);
        }
        /* sort for percentiles */
        for (int i = 1; i < ok; i++) { double v = rel[i]; int j = i - 1; while (j >= 0 && rel[j] > v) { rel[j + 1] = rel[j]; j--; } rel[j + 1] = v; }
        printf("  %d | %6d | %8d | %9d | %11.2e | %14.2e | %16.2e\n", n, ok, sing, sat, maxabs,
               ok ? rel[ok / 2] : 0.0, ok ? rel[(int)(ok * 0.99)] : 0.0);
    }

    /* 3) integer-only path */
    qmi_t Aq[4] = { (qmi_t)(-6.3 * (1 << QMI_FRAC)), 2 << QMI_FRAC, 3 << QMI_FRAC, -(4 << QMI_FRAC) };
    st = qmi_inverse_q(Aq, X2, 2);
    printf("\nInteger-input path, same example: status=%d inv = [% .6f % .6f; % .6f % .6f]\n", st,
           qmi_to_float(X2[0]), qmi_to_float(X2[1]), qmi_to_float(X2[2]), qmi_to_float(X2[3]));
    return 0;
}
