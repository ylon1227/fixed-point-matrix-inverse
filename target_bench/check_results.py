"""Check results.bin read back from the target.
1) Rebuilds the same LCG test matrices on the host.
2) Runs the same qmatinv.c on the host (via a shared library) and checks the
   target output is bit-exact.
3) Compares against a double-precision inverse (cases flagged by the library
   as saturated/singular are excluded, as the library reports them).
Build the host library first:  gcc -O2 -std=c99 -shared -o build/qmatinv_host.dll ../qmatinv.c
Usage: python check_results.py build/results.bin        License: MIT"""
import ctypes, os, struct, sys
import numpy as np

lib = os.path.abspath("build/qmatinv_host.dll" if os.name == "nt" else "build/qmatinv_host.so")
dll = ctypes.CDLL(lib)
dll.qmi_inverse_f32in.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_int32), ctypes.c_int]

b = open(sys.argv[1] if len(sys.argv) > 1 else "build/results.bin", "rb").read()
magic, trials, nmin, nmax, frac = struct.unpack_from("<5I", b, 0)
assert magic == 0xC0FFEE01, "benchmark not finished: magic=%#x" % magic
ns = nmax - nmin + 1
off = 20
def take(cnt, fmt="I"):
    global off
    v = np.array(struct.unpack_from("<%d%s" % (cnt, fmt), b, off)); off += 4 * cnt; return v
cmin, cmax, csum = (take(3 * ns).reshape(3, ns) for _ in range(3))
take(ns)                                                   # status_bad (re-derived on host)
out = take(ns * trials * nmax * nmax, "i").reshape(ns, trials, nmax * nmax)

lcg = 12345
def frand():                                               # same LCG and float ops as bench.c
    global lcg
    lcg = (lcg * 1664525 + 1013904223) & 0xFFFFFFFF
    return np.float32(np.float32(np.float32(lcg >> 8) / np.float32(16777216.0)) * np.float32(20.0)) - np.float32(10.0)

print("cycles: avg (min-max); fixed = qmatinv, soft = float GJ with soft-float, FPU = float GJ on FPv4")
for s in range(ns):
    n = nmin + s; rels = []; exact = 0; flagged = 0
    for t in range(trials):
        A = np.array([frand() for _ in range(n * n)], dtype=np.float32)
        Xh = (ctypes.c_int32 * (n * n))()
        st = dll.qmi_inverse_f32in(A.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), Xh, n)
        exact += np.array_equal(np.array(Xh[:]), out[s, t, :n * n])
        if st != 0: flagged += 1; continue
        Xd = np.linalg.inv(A.reshape(n, n).astype(np.float64))
        Xq = out[s, t, :n * n].reshape(n, n) / float(1 << frac)
        rels.append(np.linalg.norm(Xq - Xd) / np.linalg.norm(Xd))
    cyc = ["%s %6.0f (%d-%d)" % (nm, csum[m, s] / trials, cmin[m, s], cmax[m, s]) for m, nm in enumerate(("fixed", "soft", "FPU"))]
    print("n=%d | %s | bit-exact %d/%d | flagged %d | rel err median %.1e max %.1e"
          % (n, " | ".join(cyc), exact, trials, flagged, np.median(rels), max(rels)))
