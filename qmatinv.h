/**
 * qmatinv.h, Fixed-point (int32) N x N matrix inverse for MCUs without an FPU
 *
 * Gauss-Jordan elimination with partial pivoting on [A | I], entirely in
 * signed Q(31-QMI_FRAC).QMI_FRAC integer arithmetic (default Q5.26 = sign + 5 integer
 * + 26 fraction bits, sometimes written "Q6.26": range +/-32,
 * resolution 1.5e-8). No floating-point division or log/pow is used in the core.
 *
 * Range trick (no rescaling afterwards)
 *   K is chosen from the largest |a_ij| using only the IEEE-754 exponent bits,
 *   so that every entry of A / 2^K is below 2 in magnitude. The augmented matrix
 *   is built as [ A / 2^K | I / 2^K ]. Gauss-Jordan turns it into
 *   [ I | (A/2^K)^-1 (I/2^K) ] = [ I | A^-1 ], so the right half is already the
 *   inverse, the scale cancels exactly and no multiply-back is needed.
 *
 * Requirements / limits
 *   - every entry of A^-1 must fit the Q range (|x| < 2^(31-QMI_FRAC), i.e. < 32
 *     for Q5.26). Nearly singular matrices overflow; this is reported, not hidden.
 *   - n <= QMI_MAX_N (stack workspace n x 2n int32).
 *   - needs a 32x32->64-bit multiply and a 64/32 divide (software on Cortex-M0,
 *     single instructions/library calls on M3/M4; works on any C99 target).
 *
 * License: MIT
 */
#ifndef QMATINV_H
#define QMATINV_H

#include <stdint.h>

#ifndef QMI_FRAC
#define QMI_FRAC   26          /* fractional bits; 26 -> Q5.26 (range +/-32) */
#endif
#ifndef QMI_MAX_N
#define QMI_MAX_N  8
#endif

typedef int32_t qmi_t;

typedef enum {
    QMI_OK = 0,
    QMI_SINGULAR,       /* zero pivot after partial pivoting                */
    QMI_SATURATED,      /* an intermediate or result value hit the Q range  */
    QMI_BAD_INPUT       /* n out of range, or Inf/NaN in the input          */
} qmi_status_t;

/* Inverse of a row-major n x n float matrix A. Result Ainv is row-major n x n in
 * Q(QMI_FRAC). Float is only touched at the input boundary (bit inspection +
 * one conversion per entry); the elimination itself is pure integer. */
qmi_status_t qmi_inverse_f32in(const float *A, qmi_t *Ainv, int n);

/* Same, for input already in Q(QMI_FRAC) (fully integer path, no float at all). */
qmi_status_t qmi_inverse_q(const qmi_t *A, qmi_t *Ainv, int n);

/* Helpers */
static inline float qmi_to_float(qmi_t x) { return (float)x / (float)(1L << QMI_FRAC); }

#endif /* QMATINV_H */
