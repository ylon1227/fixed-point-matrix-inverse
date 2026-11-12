/* Single-precision Gauss-Jordan reference. Compiled twice:
 *   -mfloat-abi=soft   -DREF_NAME=inv_f32_soft  (what a no-FPU MCU would run)
 *   -mfloat-abi=softfp -mfpu=fpv4-sp-d16 -DREF_NAME=inv_f32_hw (Cortex-M4F FPU)
 * License: MIT */
#define MAXN 8

static float fabs_f(float x) { return x < 0.0f ? -x : x; }

int REF_NAME(const float *A, float *X, int n)
{
    float M[MAXN][2 * MAXN];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < 2 * n; j++)
            M[i][j] = j < n ? A[i * n + j] : (j - n == i ? 1.0f : 0.0f);
    for (int p = 0; p < n; p++) {
        int mr = p;
        for (int r = p + 1; r < n; r++) if (fabs_f(M[r][p]) > fabs_f(M[mr][p])) mr = r;
        if (M[mr][p] == 0.0f) return -1;
        for (int k = 0; k < 2 * n; k++) { float t = M[p][k]; M[p][k] = M[mr][k]; M[mr][k] = t; }
        float pv = M[p][p];
        for (int k = 0; k < 2 * n; k++) M[p][k] /= pv;
        for (int r = 0; r < n; r++) if (r != p) {
            float f = M[r][p];
            for (int k = 0; k < 2 * n; k++) M[r][k] -= f * M[p][k];
        }
    }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) X[i * n + j] = M[i][n + j];
    return 0;
}
