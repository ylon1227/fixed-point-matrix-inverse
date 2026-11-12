/* On-target benchmark: qmatinv vs float Gauss-Jordan (soft-float and FPU).
 * Results are written to a RAM block (section .results) and read back with a
 * debugger, so nothing (no printf/UART) runs inside the timed region.
 * Core clock: reset default HSI 16 MHz, 0 flash wait states (cycle counts are
 * what matter). License: MIT */
#include <stdint.h>
#include "../qmatinv.h"

int inv_f32_soft(const float *A, float *X, int n);
int inv_f32_hw(const float *A, float *X, int n);

#define DEMCR      (*(volatile uint32_t *)0xE000EDFC)
#define DWT_CTRL   (*(volatile uint32_t *)0xE0001000)
#define DWT_CYCCNT (*(volatile uint32_t *)0xE0001004)

#define NMIN    2
#define NMAX    6
#define TRIALS  50
#define NSIZES  (NMAX - NMIN + 1)

typedef struct {
    uint32_t magic;                       /* 0xC0FFEE01 when done            */
    uint32_t trials, nmin, nmax, frac;
    uint32_t cyc_min[3][NSIZES];          /* [0]=fixed [1]=float soft [2]=float FPU */
    uint32_t cyc_max[3][NSIZES];
    uint32_t cyc_sum[3][NSIZES];
    uint32_t status_bad[NSIZES];          /* fixed path: not QMI_OK          */
    int32_t  fixed_out[NSIZES][TRIALS][NMAX * NMAX];
    uint32_t seed_used;
} results_t;

__attribute__((section(".results"))) volatile results_t R;

static uint32_t lcg = 12345u;
static float frand(void)                  /* U(-10,10), same LCG as the host checker */
{
    lcg = lcg * 1664525u + 1013904223u;
    return ((float)(lcg >> 8) / 16777216.0f) * 20.0f - 10.0f;
}

int main(void)
{
    DEMCR |= (1u << 24);                  /* TRCENA   */
    DWT_CYCCNT = 0;
    DWT_CTRL |= 1u;                       /* CYCCNTENA */

    R.magic = 0; R.trials = TRIALS; R.nmin = NMIN; R.nmax = NMAX; R.frac = QMI_FRAC; R.seed_used = lcg;
    for (int m = 0; m < 3; m++) for (int s = 0; s < NSIZES; s++) {
        R.cyc_min[m][s] = 0xFFFFFFFFu; R.cyc_max[m][s] = 0; R.cyc_sum[m][s] = 0;
    }

    static float A[NMAX * NMAX], Xf[NMAX * NMAX];
    static qmi_t Xq[NMAX * NMAX];
    for (int n = NMIN; n <= NMAX; n++) {
        int s = n - NMIN;
        R.status_bad[s] = 0;
        for (int t = 0; t < TRIALS; t++) {
            for (int k = 0; k < n * n; k++) A[k] = frand();
            uint32_t c[3];

            uint32_t t0 = DWT_CYCCNT;
            qmi_status_t st = qmi_inverse_f32in(A, Xq, n);
            c[0] = DWT_CYCCNT - t0;

            t0 = DWT_CYCCNT; (void)inv_f32_soft(A, Xf, n); c[1] = DWT_CYCCNT - t0;
            t0 = DWT_CYCCNT; (void)inv_f32_hw(A, Xf, n);   c[2] = DWT_CYCCNT - t0;

            if (st != QMI_OK) R.status_bad[s]++;
            for (int k = 0; k < n * n; k++) R.fixed_out[s][t][k] = Xq[k];
            for (int m = 0; m < 3; m++) {
                if (c[m] < R.cyc_min[m][s]) R.cyc_min[m][s] = c[m];
                if (c[m] > R.cyc_max[m][s]) R.cyc_max[m][s] = c[m];
                R.cyc_sum[m][s] += c[m];
            }
        }
    }
    R.magic = 0xC0FFEE01u;
    for (;;) { }
}
