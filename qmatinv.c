/**
 * qmatinv.c, see qmatinv.h. License: MIT
 */
#include "qmatinv.h"

#define QMI_ONE ((int64_t)1 << QMI_FRAC)

typedef struct { uint8_t sat; } qmi_ctx_t;

static inline qmi_t sat32(qmi_ctx_t *c, int64_t x)
{
    if (x > INT32_MAX) { c->sat = 1; return INT32_MAX; }
    if (x < INT32_MIN) { c->sat = 1; return INT32_MIN; }
    return (qmi_t)x;
}

static inline qmi_t qmul(qmi_ctx_t *c, qmi_t a, qmi_t b)   /* rounded */
{
    int64_t p = (int64_t)a * (int64_t)b;
    p += (int64_t)1 << (QMI_FRAC - 1);
    return sat32(c, p >> QMI_FRAC);
}

static inline qmi_t qdiv(qmi_ctx_t *c, qmi_t a, qmi_t b)   /* rounded to nearest */
{
    int64_t num = (int64_t)a * QMI_ONE;
    int64_t half = (b > 0 ? b : -(int64_t)b) / 2;
    num += ((num >= 0) == (b > 0)) ? half : -half;
    return sat32(c, num / (int64_t)b);
}

static inline qmi_t qabs(qmi_t x) { return x == INT32_MIN ? INT32_MAX : (x < 0 ? -x : x); }

/* Core: M is n x 2n workspace = [A/2^K | I/2^K]; result left in the right half. */
static qmi_status_t gauss_jordan(qmi_ctx_t *c, qmi_t M[][2 * QMI_MAX_N], int n)
{
    for (int p = 0; p < n; p++) {
        int maxr = p;
        qmi_t maxv = qabs(M[p][p]);
        for (int r = p + 1; r < n; r++) {
            qmi_t v = qabs(M[r][p]);
            if (v > maxv) { maxv = v; maxr = r; }
        }
        if (maxv == 0) return QMI_SINGULAR;
        if (maxr != p) {
            for (int k = 0; k < 2 * n; k++) { qmi_t t = M[p][k]; M[p][k] = M[maxr][k]; M[maxr][k] = t; }
        }
        qmi_t pv = M[p][p];
        /* eliminate the other rows first, normalise the pivot row last
         * (keeps the pivot row at full magnitude while it is being used) */
        for (int r = 0; r < n; r++) {
            if (r == p) continue;
            qmi_t f = qdiv(c, M[r][p], pv);
            if (f == 0) continue;
            for (int k = 0; k < 2 * n; k++)
                M[r][k] = sat32(c, (int64_t)M[r][k] - (int64_t)qmul(c, f, M[p][k]));
        }
        for (int k = 0; k < 2 * n; k++) M[p][k] = qdiv(c, M[p][k], pv);
    }
    return QMI_OK;
}

/* K such that max|a| / 2^K < 2, from the IEEE-754 exponent only. */
static int scale_power_from_float(const float *A, int nn, int *bad)
{
    int max_exp = 0;
    for (int i = 0; i < nn; i++) {
        union { float f; uint32_t u; } v = { A[i] };
        int e = (int)((v.u & 0x7F800000u) >> 23);
        if (e == 255) { *bad = 1; return 0; }          /* Inf / NaN */
        if (e > max_exp) max_exp = e;
    }
    int true_exp = max_exp - 127;                     /* floor(log2 max|a|) */
    return (max_exp == 0 || true_exp < 0) ? 0 : true_exp + 1;
}

/* float -> integer round(x * 2^shift) using only the IEEE-754 bits (integer ops,
 * no FPU / no soft-float call). Subnormals are treated as 0. */
static qmi_t f32_bits_to_q(qmi_ctx_t *c, float x, int shift)
{
    union { float f; uint32_t u; } v = { x };
    int e = (int)((v.u >> 23) & 0xFFu);
    if (e == 0) return 0;
    int64_t m = (int64_t)((v.u & 0x7FFFFFu) | 0x800000u);   /* 24-bit mantissa */
    int s = e - 150 + shift;                                 /* value = m * 2^s */
    int64_t q;
    if (s >= 0) {
        if (s > 38) { c->sat = 1; q = INT64_MAX; }
        else q = m << s;
    } else {
        int r = -s;
        q = (r > 62) ? 0 : ((m + ((int64_t)1 << (r - 1))) >> r);
    }
    if (v.u & 0x80000000u) q = -q;
    return sat32(c, q);
}

/* K for Q input: smallest K with max|a| >> K < 2.0 in Q */
static int scale_power_from_q(const qmi_t *A, int nn)
{
    uint32_t m = 0;
    for (int i = 0; i < nn; i++) { uint32_t v = (uint32_t)qabs(A[i]); if (v > m) m = v; }
    int K = 0;
    while ((m >> K) >= (uint32_t)(2 * QMI_ONE)) K++;
    return K;
}

static qmi_status_t finish(qmi_ctx_t *c, qmi_t M[][2 * QMI_MAX_N], int n, int K, qmi_t *Ainv)
{
    qmi_t diag = (K <= QMI_FRAC) ? (qmi_t)(QMI_ONE >> K) : 0;     /* 1 / 2^K in Q */
    if (diag == 0) return QMI_SATURATED;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) M[i][n + j] = (i == j) ? diag : 0;
    qmi_status_t st = gauss_jordan(c, M, n);
    if (st != QMI_OK) return st;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) Ainv[i * n + j] = M[i][n + j];
    return c->sat ? QMI_SATURATED : QMI_OK;
}

qmi_status_t qmi_inverse_f32in(const float *A, qmi_t *Ainv, int n)
{
    if (n < 1 || n > QMI_MAX_N) return QMI_BAD_INPUT;
    int bad = 0;
    int K = scale_power_from_float(A, n * n, &bad);
    if (bad) return QMI_BAD_INPUT;

    qmi_ctx_t c = {0};
    qmi_t M[QMI_MAX_N][2 * QMI_MAX_N];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            M[i][j] = f32_bits_to_q(&c, A[i * n + j], QMI_FRAC - K);
    return finish(&c, M, n, K, Ainv);
}

qmi_status_t qmi_inverse_q(const qmi_t *A, qmi_t *Ainv, int n)
{
    if (n < 1 || n > QMI_MAX_N) return QMI_BAD_INPUT;
    int K = scale_power_from_q(A, n * n);
    qmi_ctx_t c = {0};
    qmi_t M[QMI_MAX_N][2 * QMI_MAX_N];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            qmi_t a = A[i * n + j];                      /* arithmetic shift with rounding */
            M[i][j] = K ? (qmi_t)(((int64_t)a + ((int64_t)1 << (K - 1))) >> K) : a;
        }
    return finish(&c, M, n, K, Ainv);
}
