# Fixed-point matrix inverse for MCUs without an FPU (int32, Gauss-Jordan)

Invert an N×N matrix on a microcontroller **with no floating-point unit** (Cortex-M0/M0+/M3, PowerPC e200z0, 8/16-bit parts) using only 32-bit integers with 64-bit intermediates.

- Gauss-Jordan elimination with partial pivoting on `[A | I]`.
- Signed **Q5.26** by default (sign + 5 integer + 26 fraction bits, range ±32, resolution 1.5e-8); change `QMI_FRAC` for another Q format.
- **No float division, no `log2`/`pow`, no soft-float calls** in the solver. Even the optional `float` input path converts using only the IEEE-754 bit pattern.
- Overflow and singularity are **reported, not hidden** (`QMI_SATURATED`, `QMI_SINGULAR`).
- About **1.7 KB** of code on Cortex-M3 (`-O2`), 2.0 KB on Cortex-M0. Needs only `__aeabi_ldivmod` (64/32 divide) and, on M0, `__aeabi_lmul`.
- **Verified on hardware**: bit-exact between PC and STM32F411, faster than soft-float for n ≥ 3 (see below).

## The range trick (why no rescaling is needed)

Fixed-point Gauss-Jordan overflows easily: input entries can be large, and the elimination multiplies them. So the input is scaled down first:

1. Find K from the **IEEE-754 exponent bits** of the largest |a_ij| (no `log2`): every entry of `A / 2^K` is then below 2 in magnitude.
2. Build the augmented matrix as **`[ A/2^K | I/2^K ]`** (scale the identity too).
3. Gauss-Jordan turns the left half into `I`, and the right half becomes `(A/2^K)^-1 · (I/2^K) = 2^K A^-1 · 2^-K = A^-1`.

The scale **cancels exactly**, so the right half is already `A^-1`: no multiply-back step, no extra rounding.

Two smaller details that matter in fixed point:

- The pivot row is **normalised last** (after it has been used to eliminate the other rows), so it keeps its full magnitude while it is the reference.
- Multiply and divide **round to nearest** instead of truncating.

## Usage

```c
#include "qmatinv.h"

float A[4] = { -6.3f, 2.0f,
                3.0f, -4.0f };
qmi_t Ainv[4];                        /* result in Q5.26 */

if (qmi_inverse_f32in(A, Ainv, 2) == QMI_OK) {
    /* use Ainv directly in fixed point, or qmi_to_float(Ainv[k]) for display */
}

/* fully integer input (already Q5.26): */
qmi_status_t st = qmi_inverse_q(Aq, Ainv, n);
```

`n ≤ QMI_MAX_N` (default 8). The workspace is `n × 2n` int32 on the stack.

## Accuracy (host test, `test_qmatinv.c`)

Compared against a double-precision Gauss-Jordan. Random matrices with entries uniform in (−10, 10), 20 000 per size, seed 12345:

| n | solved | singular | saturated (inverse entry ≥ 32) | max abs error | median relative error | 99th percentile relative error |
|---|---|---|---|---|---|---|
| 2 | 19 234 | 0 | 766 | 1.39e-4 | 7.8e-8 | 3.0e-7 |
| 3 | 18 803 | 0 | 1 197 | 4.42e-5 | 1.2e-7 | 3.8e-7 |
| 4 | 18 631 | 0 | 1 369 | 6.5e-6 | 1.6e-7 | 4.8e-7 |
| 5 | 18 481 | 0 | 1 519 | 9.0e-6 | 2.0e-7 | 5.9e-7 |
| 6 | 18 423 | 0 | 1 577 | 1.19e-5 | 2.4e-7 | 7.2e-7 |

Relative error = ‖X_fixed − X_double‖_F / ‖X_double‖_F. "Saturated" cases are nearly singular matrices whose inverse has an entry ≥ 32 and cannot be represented in Q5.26; they are flagged, never silently wrong. Use fewer fraction bits (e.g. `QMI_FRAC 20`, range ±2048) if your inverses are larger.

Worked example `[[-6.3, 2], [3, -4]]`: fixed result `[-0.2083333, -0.1041667; -0.1562500, -0.3281250]`, largest error 3.5e-8.

Run it yourself:

```bash
gcc -O2 -std=c99 -Wall -Wextra qmatinv.c test_qmatinv.c -o test_qmatinv -lm
./test_qmatinv
```

## Measured on hardware (STM32F411, Cortex-M4F)

`target_bench/` is a bare-metal benchmark (no HAL, nothing but the call inside the timed region, DWT cycle counter). Results are stored in RAM and read back with STM32CubeProgrammer, so no UART/printf disturbs the timing. Core at reset-default HSI 16 MHz, 0 flash wait states. 50 random matrices per size, entries U(−10, 10).

| n | fixed (this library) | float Gauss-Jordan, **soft-float** (what a no-FPU MCU runs) | float Gauss-Jordan, **FPU** (FPv4-SP) |
|---|---|---|---|
| 2 | 2 379 cycles | 2 228 | 800 |
| 3 | 5 868 | 6 237 | 1 872 |
| 4 | 11 493 | 13 500 | 3 607 |
| 5 | 19 605 | 24 985 | 6 147 |
| 6 | 30 711 | 41 664 | 9 637 |

(average cycles; min-max spread is within ±3 %.)

Correctness on the target:

- **Bit-exact with the host build: 250 / 250 results identical** (same `qmatinv.c` compiled for x86-64 and for Cortex-M4).
- Against a double-precision inverse: median relative error 0.7e-7 to 2.4e-7, worst case 1.0e-6. 9 of 250 matrices were **flagged as saturated** by the library (inverse entry beyond ±32) and excluded, as intended.

How to read this:

- **If your MCU has an FPU, use float.** On Cortex-M4F the hardware-float version is about 3× faster.
- **Without an FPU, this is faster than soft-float for n ≥ 3** (6 % at 3×3, 26 % at 6×6), slightly slower at 2×2, and its result is deterministic and bit-exact across platforms. The soft-float column was measured on the same M4 core with `-mfloat-abi=soft`, which is how a no-FPU Cortex-M3 would run it; the integer instructions used (UMULL, SDIV) are the same on M3, so the ratio should carry over. On Cortex-M0 (no 64-bit multiply instruction, no divide) both paths get slower; not measured.

Reproduce: see `target_bench/` (build commands at the top of each file, `check_results.py` re-derives the matrices on the host and does the bit-exact and error checks).

## Related work

- **CMSIS-DSP** provides `arm_mat_inverse_f32` / `f64` but, as far as we know, no fixed-point (q15/q31) matrix inverse; it is the usual starting point on Cortex-M.
- **libfixmatrix** (Q16.16) offers fixed-point matrix operations including decompositions for solving linear systems.
- Power-of-two pre-scaling is the classic *block floating point* idea from DSP.

This repository is a small, tested, documented implementation of direct Gauss-Jordan inversion in fixed point. The scaled-identity detail (scale cancels, no multiply-back) is an implementation choice, not a new method. For larger or ill-conditioned matrices prefer LU/Cholesky-based solvers and avoid forming the inverse at all.

## Files

| File | |
|---|---|
| `qmatinv.h`, `qmatinv.c` | the library (portable C99, no dependencies) |
| `test_qmatinv.c` | host accuracy test against a double reference |
| `target_bench/` | bare-metal STM32F4 benchmark (startup, linker script, float references, result checker) |

Keywords: fixed point matrix inverse, Q31 Q26 matrix inversion, Gauss-Jordan fixed point, no FPU, integer only, Cortex-M0 Cortex-M3 matrix inverse, embedded linear algebra without float, CMSIS-DSP q31 alternative

## License

MIT
