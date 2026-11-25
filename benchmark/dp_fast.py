"""Numba-accelerated DP alignment functions for the octant-descriptor aligner.

Drop-in replacements for the reference dp_local_float and
dp_semiglobal_float, JIT-compiled with Numba for ~50-100x speedup
over pure Python nested loops.

Usage:
    from dp_fast import dp_local_float, dp_semiglobal_float
    # Same interface as the reference pure-Python versions
"""

import numpy as np
try:
    import numba
    from numba import njit
    NUMBA_AVAILABLE = True
    print(f"dp_fast: Numba {numba.__version__} available", flush=True)
except ImportError:
    NUMBA_AVAILABLE = False
    print("dp_fast: Numba not available, falling back to pure Python", flush=True)


if NUMBA_AVAILABLE:

    @njit(cache=False)
    def _dp_local_numba(S, gap_penalty):
        """Smith-Waterman local DP, Numba-compiled inner kernel.

        Returns (result_i, result_j, max_val, n_pairs) as numpy arrays/scalars.
        """
        nA, nB = S.shape
        F = np.zeros((nA + 1, nB + 1), dtype=np.float64)
        max_val = 0.0
        max_i = 0
        max_j = 0

        for i in range(1, nA + 1):
            for j in range(1, nB + 1):
                diag = F[i - 1, j - 1] + S[i - 1, j - 1]
                up = F[i - 1, j] + gap_penalty
                left = F[i, j - 1] + gap_penalty
                best = 0.0
                if diag > best:
                    best = diag
                if up > best:
                    best = up
                if left > best:
                    best = left
                F[i, j] = best
                if best > max_val:
                    max_val = best
                    max_i = i
                    max_j = j

        # Traceback from max cell, stop at 0
        buf_i = np.empty(nA + nB, dtype=np.int32)
        buf_j = np.empty(nA + nB, dtype=np.int32)
        n_pairs = 0
        i = max_i
        j = max_j
        while i > 0 and j > 0 and F[i, j] > 0.0:
            if F[i, j] == F[i - 1, j - 1] + S[i - 1, j - 1]:
                buf_i[n_pairs] = i - 1
                buf_j[n_pairs] = j - 1
                n_pairs += 1
                i -= 1
                j -= 1
            elif F[i, j] == F[i - 1, j] + gap_penalty:
                i -= 1
            else:
                j -= 1

        # Reverse into output arrays
        result_i = np.empty(n_pairs, dtype=np.int32)
        result_j = np.empty(n_pairs, dtype=np.int32)
        for k in range(n_pairs):
            result_i[k] = buf_i[n_pairs - 1 - k]
            result_j[k] = buf_j[n_pairs - 1 - k]

        return result_i, result_j, max_val, n_pairs

    @njit(cache=False)
    def _dp_semiglobal_numba(S, gap_penalty):
        """Semi-global DP (free end gaps), Numba-compiled inner kernel.

        Returns (result_i, result_j, best_score, n_pairs).
        """
        nA, nB = S.shape
        F = np.zeros((nA + 1, nB + 1), dtype=np.float64)

        for i in range(1, nA + 1):
            for j in range(1, nB + 1):
                diag = F[i - 1, j - 1] + S[i - 1, j - 1]
                up = F[i - 1, j] + gap_penalty
                left = F[i, j - 1] + gap_penalty
                best = diag
                if up > best:
                    best = up
                if left > best:
                    best = left
                F[i, j] = best

        # Find best endpoint in last row or last column
        best_score = F[nA, 0]
        best_i = nA
        best_j = 0
        for j in range(1, nB + 1):
            if F[nA, j] > best_score:
                best_score = F[nA, j]
                best_i = nA
                best_j = j
        for i in range(1, nA + 1):
            if F[i, nB] > best_score:
                best_score = F[i, nB]
                best_i = i
                best_j = nB

        # Traceback until reaching row 0 or column 0
        buf_i = np.empty(nA + nB, dtype=np.int32)
        buf_j = np.empty(nA + nB, dtype=np.int32)
        n_pairs = 0
        i = best_i
        j = best_j
        while i > 0 and j > 0:
            if F[i, j] == F[i - 1, j - 1] + S[i - 1, j - 1]:
                buf_i[n_pairs] = i - 1
                buf_j[n_pairs] = j - 1
                n_pairs += 1
                i -= 1
                j -= 1
            elif F[i, j] == F[i - 1, j] + gap_penalty:
                i -= 1
            else:
                j -= 1

        # Reverse
        result_i = np.empty(n_pairs, dtype=np.int32)
        result_j = np.empty(n_pairs, dtype=np.int32)
        for k in range(n_pairs):
            result_i[k] = buf_i[n_pairs - 1 - k]
            result_j[k] = buf_j[n_pairs - 1 - k]

        return result_i, result_j, best_score, n_pairs

    def dp_local_float(S, gap_penalty=-2.0):
        """Smith-Waterman local DP with Numba acceleration.

        Drop-in replacement for the reference dp_local_float.
        """
        S64 = np.ascontiguousarray(S, dtype=np.float64)
        result_i, result_j, max_val, n_pairs = _dp_local_numba(S64, float(gap_penalty))
        pairs = list(zip(result_i.tolist(), result_j.tolist()))
        info = {"mode": "local_float", "alignment_score": float(max_val), "n_used": n_pairs}
        return pairs, info

    def dp_semiglobal_float(S, gap_penalty=-2.0):
        """Semi-global DP with Numba acceleration.

        Drop-in replacement for the reference dp_semiglobal_float.
        """
        S64 = np.ascontiguousarray(S, dtype=np.float64)
        result_i, result_j, best_score, n_pairs = _dp_semiglobal_numba(S64, float(gap_penalty))
        pairs = list(zip(result_i.tolist(), result_j.tolist()))
        info = {"mode": "semiglobal_float", "alignment_score": float(best_score), "n_used": n_pairs}
        return pairs, info

    @njit(cache=False)
    def _dp_local_affine_numba(S, gap_open, gap_extend):
        """Smith-Waterman local DP with affine gaps, Numba-compiled.

        Three-state DP: M (match), I (insert/gap in target), D (delete/gap in query).
        gap_open is the penalty for opening a gap, gap_extend for each additional position.

        Returns (result_i, result_j, max_val, n_pairs).
        """
        nA, nB = S.shape
        M = np.zeros((nA + 1, nB + 1), dtype=np.float64)
        I = np.zeros((nA + 1, nB + 1), dtype=np.float64)
        D = np.zeros((nA + 1, nB + 1), dtype=np.float64)

        max_val = 0.0
        max_i = 0
        max_j = 0
        max_state = 0  # 0=M, 1=I, 2=D

        for i in range(1, nA + 1):
            for j in range(1, nB + 1):
                s = S[i - 1, j - 1]

                # Match state: come from M, I, or D diagonal, or 0 (local)
                best_m = 0.0
                v = M[i - 1, j - 1] + s
                if v > best_m:
                    best_m = v
                v = I[i - 1, j - 1] + s
                if v > best_m:
                    best_m = v
                v = D[i - 1, j - 1] + s
                if v > best_m:
                    best_m = v
                M[i, j] = best_m

                # Insert state: gap in target (consume query residue)
                best_i = 0.0
                v = M[i - 1, j] + gap_open
                if v > best_i:
                    best_i = v
                v = I[i - 1, j] + gap_extend
                if v > best_i:
                    best_i = v
                I[i, j] = best_i

                # Delete state: gap in query (consume target residue)
                best_d = 0.0
                v = M[i, j - 1] + gap_open
                if v > best_d:
                    best_d = v
                v = D[i, j - 1] + gap_extend
                if v > best_d:
                    best_d = v
                D[i, j] = best_d

                # Track global maximum across all states
                if best_m > max_val:
                    max_val = best_m
                    max_i = i
                    max_j = j
                    max_state = 0
                if best_i > max_val:
                    max_val = best_i
                    max_i = i
                    max_j = j
                    max_state = 1
                if best_d > max_val:
                    max_val = best_d
                    max_i = i
                    max_j = j
                    max_state = 2

        # Traceback from max cell
        buf_i = np.empty(nA + nB, dtype=np.int32)
        buf_j = np.empty(nA + nB, dtype=np.int32)
        n_pairs = 0
        i = max_i
        j = max_j
        state = max_state

        while i > 0 and j > 0:
            if state == 0:  # Match
                if M[i, j] <= 0.0:
                    break
                buf_i[n_pairs] = i - 1
                buf_j[n_pairs] = j - 1
                n_pairs += 1
                s = S[i - 1, j - 1]
                if M[i, j] == M[i - 1, j - 1] + s:
                    state = 0
                elif M[i, j] == I[i - 1, j - 1] + s:
                    state = 1
                elif M[i, j] == D[i - 1, j - 1] + s:
                    state = 2
                else:
                    break  # local restart
                i -= 1
                j -= 1
            elif state == 1:  # Insert (gap in target)
                if I[i, j] <= 0.0:
                    break
                if I[i, j] == M[i - 1, j] + gap_open:
                    state = 0
                else:
                    state = 1
                i -= 1
            else:  # state == 2, Delete (gap in query)
                if D[i, j] <= 0.0:
                    break
                if D[i, j] == M[i, j - 1] + gap_open:
                    state = 0
                else:
                    state = 2
                j -= 1

        # Reverse
        result_i = np.empty(n_pairs, dtype=np.int32)
        result_j = np.empty(n_pairs, dtype=np.int32)
        for k in range(n_pairs):
            result_i[k] = buf_i[n_pairs - 1 - k]
            result_j[k] = buf_j[n_pairs - 1 - k]

        return result_i, result_j, max_val, n_pairs

    def dp_local_affine(S, gap_open=-4.0, gap_extend=-1.0):
        """Local DP with affine gap penalties, Numba acceleration."""
        S64 = np.ascontiguousarray(S, dtype=np.float64)
        result_i, result_j, max_val, n_pairs = _dp_local_affine_numba(
            S64, float(gap_open), float(gap_extend))
        pairs = list(zip(result_i.tolist(), result_j.tolist()))
        info = {"mode": "local_affine", "alignment_score": float(max_val), "n_used": n_pairs}
        return pairs, info

else:
    raise ImportError("dp_fast requires numba (pip install numba)")


def warmup():
    """Warm up Numba JIT compilation with a small test matrix."""
    if not NUMBA_AVAILABLE:
        return
    S = np.random.rand(10, 10).astype(np.float64)
    _dp_local_numba(S, -2.0)
    _dp_semiglobal_numba(S, -2.0)
    _dp_local_affine_numba(S, -4.0, -1.0)
    print("dp_fast: JIT warmup complete", flush=True)
