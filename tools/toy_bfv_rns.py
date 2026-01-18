# -*- coding: utf-8 -*-
"""
A tiny hand-checkable BFV(RNS) toy example in pure Python with small polynomials,
DCRT (multi-prime) and NTT-based negacyclic convolution, plus EvalMult+Relinearization.

Security: none. This is purely didactic and sized for hand verification.
Ring: R_q = Z_q[x]/(x^n+1) with n=4 and q_i chosen so 2n | (q_i-1) for NTT.
"""
from typing import List, Tuple
from math import prod

# ------------------------------ helpers -------------------------------------

def egcd(a: int, b: int):
    if b == 0:
        return (1, 0, a)
    x, y, g = egcd(b, a % b)
    return (y, x - (a // b) * y, g)


def invmod(a: int, m: int) -> int:
    a %= m
    x, y, g = egcd(a, m)
    if g != 1:
        raise ValueError("no inverse")
    return x % m


def powmod(a: int, e: int, m: int) -> int:
    r = 1
    a %= m
    while e > 0:
        if e & 1:
            r = (r * a) % m
        a = (a * a) % m
        e >>= 1
    return r

# ----------------------------- CRT on Q -------------------------------------

def crt_pack(residues: List[int], moduli: List[int]) -> int:
    Q = prod(moduli)
    x = 0
    for ri, qi in zip(residues, moduli):
        Mi = Q // qi
        Mi_inv = invmod(Mi % qi, qi)
        x = (x + ri * Mi * Mi_inv) % Q
    return x


def crt_unpack(x: int, moduli: List[int]) -> List[int]:
    return [x % qi for qi in moduli]

# --------------------------- small NTT (negacyclic) -------------------------

def bit_reverse_copy(a: List[int]) -> List[int]:
    n = len(a)
    j = 0
    b = [0] * n
    for i in range(n):
        b[j] = a[i]
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j ^= bit
    return b


def ntt(a: List[int], root: int, q: int) -> List[int]:
    # iterative Cooley¨CTukey, assumes root is n-th primitive root
    n = len(a)
    A = bit_reverse_copy([x % q for x in a])
    m = 1
    while m < n:
        mh = m
        m <<= 1
        w_m = powmod(root, n // m, q)
        for k in range(0, n, m):
            w = 1
            for j in range(mh):
                t = (w * A[k + j + mh]) % q
                u = A[k + j]
                A[k + j] = (u + t) % q
                A[k + j + mh] = (u - t) % q
                w = (w * w_m) % q
    return A


def intt(A: List[int], root: int, q: int) -> List[int]:
    n = len(A)
    root_inv = invmod(root, q)
    a = bit_reverse_copy([x % q for x in A])
    m = 1
    while m < n:
        mh = m
        m <<= 1
        w_m = powmod(root_inv, n // m, q)
        for k in range(0, n, m):
            w = 1
            for j in range(mh):
                t = (w * a[k + j + mh]) % q
                u = a[k + j]
                a[k + j] = (u + t) % q
                a[k + j + mh] = (u - t) % q
                w = (w * w_m) % q
    n_inv = invmod(n, q)
    return [(x * n_inv) % q for x in a]


def negacyclic_convolution_ntt(a: List[int], b: List[int], q: int, n: int, psi: int) -> List[int]:
    # Using Harvey's twist: psi is primitive 2n-th root; omega=psi^2 primitive n-th root.
    omega = (psi * psi) % q
    # pre-twist by psi^i
    psi_pows = [1] * n
    for i in range(1, n):
        psi_pows[i] = (psi_pows[i - 1] * psi) % q
    a_t = [(a[i] * psi_pows[i]) % q for i in range(n)]
    b_t = [(b[i] * psi_pows[i]) % q for i in range(n)]
    # NTT
    A = ntt(a_t, omega, q)
    B = ntt(b_t, omega, q)
    C = [(A[i] * B[i]) % q for i in range(n)]
    c_t = intt(C, omega, q)
    # post-twist by psi^{-i}
    psi_inv = invmod(psi, q)
    psi_inv_pows = [1] * n
    for i in range(1, n):
        psi_inv_pows[i] = (psi_inv_pows[i - 1] * psi_inv) % q
    c = [(c_t[i] * psi_inv_pows[i]) % q for i in range(n)]
    return c


def schoolbook_negacyclic(a: List[int], b: List[int], q: int, n: int) -> List[int]:
    res = [0] * n
    for i in range(n):
        for j in range(n):
            k = i + j
            v = (a[i] * b[j]) % q
            if k >= n:
                k -= n
                v = (-v) % q
            res[k] = (res[k] + v) % q
    return res

# --------------------------- DCRT small context -----------------------------

class DCRTContext:
    def __init__(self, n: int, primes: List[int]):
        # primes must satisfy 2n | (p-1)
        self.n = n
        self.qs = primes
        self.Q = prod(primes)
        # pick psi per prime: psi = g^((p-1)/(2n)) where g is primitive generator
        self.psi = []
        for p in primes:
            # find a primitive root g of Z_p*
            g = self._primitive_root(p)
            psi = powmod(g, (p - 1) // (2 * n), p)
            # ensure correct order (psi^(2n)=1 and psi^n = -1)
            assert powmod(psi, 2 * n, p) == 1
            assert (powmod(psi, n, p) + 1) % p == 0
            self.psi.append(psi)

    def _primitive_root(self, p: int) -> int:
        # naive primitive root finder for small p
        # factor p-1
        m = p - 1
        fac = []
        x = 2
        t = m
        while x * x <= t:
            if t % x == 0:
                fac.append(x)
                while t % x == 0:
                    t //= x
            x += 1
        if t > 1:
            fac.append(t)
        for g in range(2, p):
            ok = True
            for f in fac:
                if powmod(g, m // f, p) == 1:
                    ok = False
                    break
            if ok:
                return g
        raise RuntimeError("no primitive root found")

# --------------------------- parameters -------------------------------------

def toy_params():
    # Plaintext modulus (small)
    t = 13
    # n must be power of two; choose small n=4 so 2n=8
    n = 4
    # choose primes where 8 | (p-1)
    qs = [17, 41]  # Q=697
    ctx = DCRTContext(n, qs)
    return t, ctx

# ---------------------------- poly utils ------------------------------------

def poly_add_mod(a: List[int], b: List[int], q: int) -> List[int]:
    return [((x + y) % q) for x, y in zip(a, b)]


def poly_sub_mod(a: List[int], b: List[int], q: int) -> List[int]:
    return [((x - y) % q) for x, y in zip(a, b)]


def poly_mul_mod(a: List[int], b: List[int], q: int, n: int, psi: int) -> List[int]:
    # NTT-based negacyclic convolution
    return negacyclic_convolution_ntt(a, b, q, n, psi)

# ------------------------------ BFV toy -------------------------------------

def keygen(t: int, ctx: DCRTContext):
    # secret ternary polynomial s (fixed small for determinism)
    n = ctx.n
    s = [1, 0, 1, 0][:n]
    # a uniform small (fixed) per prime; e small error poly
    a = [[2, 3, 1, 4][:n] for _ in ctx.qs]
    e = [[1, 0, 0, 0][:n] for _ in ctx.qs]
    # public key per prime: b = -a*s + e  (negacyclic conv)
    b = []
    for idx, q in enumerate(ctx.qs):
        as_ = poly_mul_mod(a[idx], s, q, ctx.n, ctx.psi[idx])
        b_q = poly_add_mod([(q - x) % q for x in as_], e[idx], q)
        b.append(b_q)
    pk = (a, b)
    # relinearization key (single-digit toy): encrypt s^2 under s
    # pick a_r,e_r small fixed
    a_r = [[5, 1, 2, 0][:n] for _ in ctx.qs]
    e_r = [[0, 1, 0, 0][:n] for _ in ctx.qs]
    # s^2 (negacyclic) over integers modulo q
    s2 = []
    for idx, q in enumerate(ctx.qs):
        s2.append(poly_mul_mod(s, s, q, ctx.n, ctx.psi[idx]))
    # b_r = -a_r*s + e_r + s^2
    b_r = []
    for idx, q in enumerate(ctx.qs):
        ar_s = poly_mul_mod(a_r[idx], s, q, ctx.n, ctx.psi[idx])
        minus_ar_s = [(q - x) % q for x in ar_s]
        tmp = poly_add_mod(minus_ar_s, e_r[idx], q)
        b_r.append(poly_add_mod(tmp, s2[idx], q))
    rlk = (a_r, b_r)
    return s, pk, rlk


def encode_poly(message_coeffs: List[int], t: int, ctx: DCRTContext) -> List[List[int]]:
    # message in R_t; scale by Delta and lift to R_Q, then split to residues
    Q = ctx.Q
    Delta = Q // t
    # scale each coefficient
    lift = [(m % t) * Delta % Q for m in message_coeffs]
    # CRT-unpack per coefficient onto each prime (vectorized per prime)
    # For simplicity, pack coefficient-wise using crt_unpack per value
    residues_per_prime = [[0] * ctx.n for _ in ctx.qs]
    for i in range(ctx.n):
        r = crt_unpack(lift[i], ctx.qs)
        for pi in range(len(ctx.qs)):
            residues_per_prime[pi][i] = r[pi]
    return residues_per_prime  # list per prime of poly coefficients


def encrypt(message_coeffs: List[int], t: int, ctx: DCRTContext, pk) -> Tuple[List[List[int]], List[List[int]]]:
    a, b = pk
    # sample u,e1,e2 as tiny fixed polys
    u = [1, 0, 0, 0][:ctx.n]
    e1 = [[1, 0, 0, 0][:ctx.n] for _ in ctx.qs]
    e2 = [[0, 1, 0, 0][:ctx.n] for _ in ctx.qs]
    m_scaled = encode_poly(message_coeffs, t, ctx)
    c0, c1 = [], []
    for idx, q in enumerate(ctx.qs):
        bu = poly_mul_mod(b[idx], u, q, ctx.n, ctx.psi[idx])
        au = poly_mul_mod(a[idx], u, q, ctx.n, ctx.psi[idx])
        c0_q = poly_add_mod(poly_add_mod(bu, e1[idx], q), m_scaled[idx], q)
        c1_q = poly_add_mod(au, e2[idx], q)
        c0.append(c0_q)
        c1.append(c1_q)
    return (c0, c1)


def eval_add(ctA, ctB, ctx: DCRTContext):
    c0 = []
    c1 = []
    for idx, q in enumerate(ctx.qs):
        c0.append(poly_add_mod(ctA[0][idx], ctB[0][idx], q))
        c1.append(poly_add_mod(ctA[1][idx], ctB[1][idx], q))
    return (c0, c1)


def eval_mult_and_relin(ctA, ctB, ctx: DCRTContext, rlk):
    a_r, b_r = rlk
    # multiply to get (d0, d1, d2)
    d0, d1, d2 = [], [], []
    for idx, q in enumerate(ctx.qs):
        psi = ctx.psi[idx]
        c0a, c1a = ctA[0][idx], ctA[1][idx]
        c0b, c1b = ctB[0][idx], ctB[1][idx]
        d0_q = poly_mul_mod(c0a, c0b, q, ctx.n, psi)
        t1 = poly_mul_mod(c0a, c1b, q, ctx.n, psi)
        t2 = poly_mul_mod(c1a, c0b, q, ctx.n, psi)
        d1_q = poly_add_mod(t1, t2, q)
        d2_q = poly_mul_mod(c1a, c1b, q, ctx.n, psi)
        d0.append(d0_q)
        d1.append(d1_q)
        d2.append(d2_q)
    # relinearize with single-key toy: c0' = d0 + d2*b_r; c1' = d1 + d2*a_r
    c0p, c1p = [], []
    for idx, q in enumerate(ctx.qs):
        psi = ctx.psi[idx]
        d2_br = poly_mul_mod(d2[idx], b_r[idx], q, ctx.n, psi)
        d2_ar = poly_mul_mod(d2[idx], a_r[idx], q, ctx.n, psi)
        c0p.append(poly_add_mod(d0[idx], d2_br, q))
        c1p.append(poly_add_mod(d1[idx], d2_ar, q))
    return (c0p, c1p)


def decrypt(ct, s: List[int], t: int, ctx: DCRTContext) -> List[int]:
    # m' = round(t/Q * (c0 + c1*s)) mod t (coefficient-wise)
    m_coeffs = [0] * ctx.n
    for i in range(ctx.n):
        # reconstruct coefficient i across primes, after combining c0 + c1*s (negacyclic)
        # build per-prime value of coeff ia
        vals = []
        for idx, q in enumerate(ctx.qs):
            psi = ctx.psi[idx]
            c0_q = ct[0][idx]
            c1_q = ct[1][idx]
            # compute (c1*s) via convolution then pick coeff i
            c1s = poly_mul_mod(c1_q, s, q, ctx.n, psi)
            v = (c0_q[i] + c1s[i]) % q
            vals.append(v)
        x = crt_pack(vals, ctx.qs)
        Q = ctx.Q
        m_scaled = (x * t + Q // 2) // Q
        m_coeffs[i] = int(m_scaled % t)
    return m_coeffs

# ------------------------------ demo ----------------------------------------
if __name__ == "__main__":
    t, ctx = toy_params()
    print(f"t={t}, qs={ctx.qs}, Q={ctx.Q}, n={ctx.n}")

    s, pk, rlk = keygen(t, ctx)
    print(f"secret s={s}")

    # message polys (small): m(x) = 5 + 2x, n=4
    m1 = [5, 2, 0, 0]
    m2 = [9, 1, 0, 0]

    ct1 = encrypt(m1, t, ctx, pk)
    ct2 = encrypt(m2, t, ctx, pk)

    dec1 = decrypt(ct1, s, t, ctx)
    dec2 = decrypt(ct2, s, t, ctx)
    print("Dec(m1)=", dec1, "expect:", [c % t for c in m1])
    print("Dec(m2)=", dec2, "expect:", [c % t for c in m2])

    # EvalAdd
    ct_add = eval_add(ct1, ct2, ctx)
    dec_add = decrypt(ct_add, s, t, ctx)
    expect_add = [((m1[i] + m2[i]) % t) for i in range(ctx.n)]
    print("Dec(add)=", dec_add, "expect:", expect_add)

    # EvalMult + Relinearize
    ct_mul = eval_mult_and_relin(ct1, ct2, ctx, rlk)
    dec_mul = decrypt(ct_mul, s, t, ctx)
    # expected negacyclic product in R_t
    expect_mul = schoolbook_negacyclic(m1, m2, t, ctx.n)
    print("Dec(mul)=", dec_mul, "expect:", expect_mul)
