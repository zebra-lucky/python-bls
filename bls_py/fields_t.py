# -*- coding: utf-8 -*-
#
# Copyright 2018 Chia Network Inc
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging


bls12381_q = Q = int('0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf'
                     '6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab', 16)
bls12381_n = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
bls12381_b = 0x4
bls12381_x = -0xd201000000010000
bls12381_nx = 0xd201000000010000


# twist/untwist inverse constants
tw1 = int('0xd0088f51cbff34d258dd3db21a5d66bb23ba5c279c2895'
          'fb39869507b587b120f55ffff58a9ffffdcff7fffffffd556', 16)
tw2 = int('0xd0088f51cbff34d258dd3db21a5d66bb23ba5c279c2895'
          'fb39869507b587b120f55ffff58a9ffffdcff7fffffffd555', 16)


FQ2_ROOT = -1 % Q
FQ2_ONE_TUPLE = (1, 0)
FQ2_ZERO_TUPLE = (0, 0)
FQ6_ROOT_TUPLE = (1, 1)
FQ6_ONE_TUPLE = (1, 0, 0, 0, 0, 0)
FQ6_ZERO_TUPLE = (0, 0, 0, 0, 0, 0)
FQ12_ROOT_TUPLE = (0, 0, 1, 0, 0, 0)
FQ12_ONE_TUPLE = (1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
FQ12_ZERO_TUPLE = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
FINAL_EXP_E = (pow(Q, 4) - pow(Q, 2) + 1) // bls12381_n


def fq_invert(P, a):
    '''Ivnert int value using extended euclidian algorithm for inversion'''
    p = P
    x0, x1, y0, y1 = 1, 0, 0, 1
    while p != 0:
        q, a, p = a // p, p, a % p
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return x0 % P


def fq_pow(P, a, X):
    '''Pow a to X mod P '''
    if X == 0:
        return 1 % P
    res = 1
    while X > 0:
        if X & 1:
            res = res * a % P
        X >>= 1
        a = a * a % P
    return res


def fq_floordiv(P, a, X):
    return a * fq_invert(P, X) % P


def fq2_neg(t_a):
    '''Neg tuple t_a returning tuple'''
    a, b = t_a
    return (-a % Q, -b % Q)


def fq2_invert(t_a):
    '''Invert tuple t_a returning tuple'''
    a, b = t_a
    factor = fq_invert(Q, a * a + b * b)
    return ((a*factor) % Q, (-b*factor) % Q)


def fq2_floordiv(t_a, t_x):
    return fq2_mul(t_a, fq2_invert(t_x))


def fq2_pow(t_a, e):
    '''Pow tuple t_a returning tuple'''
    m, n = t_a
    a, b = 1, 0
    while e:
        if e & 1:
            a, b = (a*m - b*n) % Q, (a*n + b*m) % Q
        m, n = (m*m - n*n) % Q, (m*n + n*m) % Q
        e >>= 1
    return (a, b)


def fq2_qi_pow(t_x, i):
    '''Calc qi_power on t_x tuple returning tuple'''
    global frob_coeffs
    i %= 2
    if i == 0:
        return t_x
    return (t_x[0], t_x[1]*frob_coeffs[2, i, 1] % Q)


def fq2_mul_by_nonresidue(t_a):
    '''Mul by nonresidue on tuple t_a returning tuple'''
    a, b = t_a
    return ((a-b) % Q, (a+b) % Q)


def fq2_add_fq(t_a, m):
    '''Add tuple t_a and int m returning tuple'''
    a, b = t_a
    return ((a+m) % Q, b)


def fq2_add(t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b = t_a
    m, n = t_m
    return ((a+m) % Q, (b+n) % Q)


def fq_sub_fq2(a, t_m):
    '''Sub tuple t_m from int a returning tuple'''
    m, n = t_m
    return ((a-m) % Q, -n % Q)


def fq2_sub_fq(t_a, m):
    '''Sub int m from tuple t_a returning tuple'''
    a, b = t_a
    return ((a-m) % Q, b)


def fq2_sub(t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b = t_a
    m, n = t_m
    return ((a-m) % Q, (b-n) % Q)


def fq2_mul_fq(t_a, m):
    '''Multiple tuple t_a on int m returning tuple'''
    a, b = t_a
    return (a*m % Q, b*m % Q)


def fq2_mul(t_a, t_m):
    '''Multiple tuple t_a on tuple t_m returning tuple'''
    a, b = t_a
    m, n = t_m
    return ((a*m - b*n) % Q, (a*n + b*m) % Q)


def fq6_neg(t_a):
    '''Neg tuple t_a returning tuple'''
    a, b, c, d, e, f = t_a
    return (-a % Q, -b % Q, -c % Q, -d % Q, -e % Q, -f % Q)


def fq6_invert(t_x):
    '''Invert tuple t_a returning tuple'''
    a, b, c = t_x[:2], t_x[2:4], t_x[4:]
    g0 = fq2_mul(a, a)
    g0 = fq2_sub(g0, fq2_mul(b, fq2_mul_by_nonresidue(c)))
    g1 = fq2_sub(fq2_mul_by_nonresidue(fq2_mul(c, c)),
                     fq2_mul(a, b))
    g2 = fq2_sub(fq2_mul(b, b), fq2_mul(a, c))
    g0a = fq2_mul(g0, a)
    g1cpg2b = fq2_add(fq2_mul(g1, c), fq2_mul(g2, b))
    factor = fq2_invert(fq2_add(g0a, fq2_mul_by_nonresidue(g1cpg2b)))
    ar, br = fq2_mul(g0, factor)
    cr, dr = fq2_mul(g1, factor)
    er, fr = fq2_mul(g2, factor)
    return (ar % Q, br % Q, cr % Q, dr % Q, er % Q, fr % Q)


def fq6_floordiv(t_a, t_x):
    return fq6_mul(t_a, fq6_invert(t_x))


def fq6_pow(t_a, e):
    '''Pow tuple t_a returning tuple'''
    t_ans = FQ6_ONE_TUPLE
    while e:
        if e & 1:
            t_ans = fq6_mul(t_ans, t_a)
        t_a = fq6_mul(t_a, t_a)
        e >>= 1
    a, b, c, d, e, f = t_ans
    return (a % Q, b % Q, c % Q, d % Q, e % Q, f % Q)


def fq6_qi_pow(t_x, i):
    '''Calc qi_power on t_x tuple returning tuple'''
    global frob_coeffs
    i %= 6
    if i == 0:
        return t_x
    a, b = fq2_qi_pow(t_x[:2], i)
    c, d = fq2_mul(fq2_qi_pow(t_x[2:4], i), frob_coeffs[6, i, 1])
    e, f = fq2_mul(fq2_qi_pow(t_x[4:6], i), frob_coeffs[6, i, 2])
    return (a, b, c, d, e, f)


def fq6_mul_by_nonresidue(t_x):
    '''Mul by nonresidue on tuple t_a returning tuple'''
    ar, br = fq2_mul(t_x[4:], FQ6_ROOT_TUPLE)
    cr, dr = t_x[:2]
    er, fr = t_x[2:4]
    return (ar, br, cr, dr, er, fr)


def fq6_add_fq(t_a, m):
    '''Add tuple t_a and int m returning tuple'''
    a, b, c, d, e, f = t_a
    return ((a+m) % Q, b, c, d, e, f)


def fq6_add_fq2(t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b, c, d, e, f = t_a
    m, n = t_m
    return ((a+m) % Q, (b+n) % Q, c, d, e, f)


def fq6_add(t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b, c, d, e, f = t_a
    m, n, o, p, q, r = t_m
    return ((a+m) % Q, (b+n) % Q, (c+o) % Q,
            (d+p) % Q, (e+q) % Q, (f+r) % Q)


def fq_sub_fq6(a, t_m):
    '''Sub tuple t_m from int a returning tuple'''
    m, n, o, p, q, r = t_m
    return ((a-m) % Q, -n % Q, -o % Q, -p % Q, -q % Q, -r % Q)


def fq6_sub_fq(t_a, m):
    '''Sub int m from tuple t_a returning tuple'''
    a, b, c, d, e, f = t_a
    return ((a-m) % Q, b, c, d, e, f)


def fq2_sub_fq6(t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b = t_a
    m, n, o, p, q, r = t_m
    return ((a-m) % Q, (b-n) % Q, -o % Q, -p % Q, -q % Q, -r % Q)


def fq6_sub_fq2(t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f = t_a
    m, n = t_m
    return ((a-m) % Q, (b-n) % Q, c, d, e, f)


def fq6_sub(t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f = t_a
    m, n, o, p, q, r = t_m
    return ((a-m) % Q, (b-n) % Q, (c-o) % Q,
            (d-p) % Q, (e-q) % Q, (f-r) % Q)


def fq6_mul_fq(t_a, m):
    '''Multiple tuple t_a on int m returning tuple'''
    a, b, c, d, e, f = t_a
    return (a*m % Q, b*m % Q, c*m % Q, d*m % Q, e*m % Q, f*m % Q)


def fq6_mul_fq2(t_a, t_m):
    '''Multiple tuple t_a on tuple t_m returning tuple'''
    a, b, c, d, e, f = t_a
    m, n = t_m
    return ((a*m - b*n) % Q, (a*n + b*m) % Q,
            (c*m - d*n) % Q, (c*n + d*m) % Q,
            (e*m - f*n) % Q, (e*n + f*m) % Q)


def fq6_mul(t_a, t_m):
    '''Multiple tuple t_a on tuple t_m returning tuple'''
    a, b, c, d, e, f = t_a
    m, n, o, p, q, r = t_m
    #0 am - bn + cq - cr - dr - dq + eo - ep - fp - fo
    #1 an + bm + cq + cr - dr + dq + eo + ep - fp + fo
    #2 ao - bp + cm - dn + eq - er - fr - fq
    #3 ap + bo + cn + dm + eq + er - fr + fq
    #4 aq - br + co - dp + em - fn
    #5 ar + bq + cp + do + en + fm

    am = a*m; an = a*n; ao = a*o; ap = a*p; aq = a*q; ar = a*r;
    bn = b*n; bm = b*m; bp = b*p; bo = b*o; br = b*r; bq = b*q;
    cq = c*q; cm = c*m; cn = c*n; co = c*o; cp = c*p; cr = c*r;
    dn = d*n; dm = d*m; dp = d*p; do = d*o; dr = d*r; dq = d*q;
    eq = e*q; em = e*m; en = e*n; eo = e*o; ep = e*p; er = e*r;
    fn = f*n; fm = f*m; fq = f*q; fp = f*p; fo = f*o; fr = f*r

    r0 = am - bn + cq - cr - dr - dq + eo - ep - fp - fo
    r1 = an + bm + cq + cr - dr + dq + eo + ep - fp + fo
    r2 = ao - bp + cm - dn + eq - er - fr - fq
    r3 = ap + bo + cn + dm + eq + er - fr + fq
    r4 = aq - br + co - dp + em - fn
    r5 = ar + bq + cp + do + en + fm

    return (r0 % Q, r1 % Q, r2 % Q, r3 % Q, r4 % Q, r5 % Q)


def fq12_neg(t_a):
    '''Neg tuple t_a returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    return (-a % Q, -b % Q, -c % Q, -d % Q, -e % Q, -f % Q,
            -g % Q, -h % Q, -i % Q, -j % Q, -k % Q, -l % Q)


def fq12_invert(t_x):
    '''Invert tuple t_a returning tuple'''
    a, b = t_x[:6], t_x[6:12]
    aa = fq6_mul(a, a)
    bb = fq6_mul(b, b)
    factor = fq6_invert(fq6_sub(aa, fq6_mul_by_nonresidue(bb)))
    ar, br, cr, dr, er, fr = fq6_mul(a, factor)
    gr, hr, ir, jr, kr, lr = fq6_mul(fq6_neg(b), factor)
    return (ar % Q, br % Q, cr % Q, dr % Q, er % Q, fr % Q,
            gr % Q, hr % Q, ir % Q, jr % Q, kr % Q, lr % Q)


def fq12_floordiv(t_a, t_x):
    return fq12_mul(t_a, fq12_invert(t_x))


def fq12_pow(t_a, e):
    '''Pow tuple t_a returning tuple'''
    t_ans = FQ12_ONE_TUPLE
    while e:
        if e & 1:
            t_ans = fq12_mul(t_ans, t_a)
        t_a = fq12_mul(t_a, t_a)
        e >>= 1
    return t_ans


def fq12_qi_pow(t_x, i):
    '''Calc qi_power on t_x tuple returning tuple'''
    global frob_coeffs
    i %= 12
    if i == 0:
        return t_x
    a, b, c, d, e, f = fq6_qi_pow(t_x[:6], i)
    g, h, i, j, k, l = fq6_mul(fq6_qi_pow(t_x[6:12], i),
                                   frob_coeffs[12, i, 1])
    return (a, b, c, d, e, f, g, h, i, j, k, l)


def fq12_add_fq(t_a, m):
    '''Add tuple t_a and int m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    return ((a+m) % Q, b, c, d, e, f, g, h, i, j, k, l)


def fq12_add_fq2(t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n = t_m
    return ((a+m) % Q, (b+n) % Q, c, d, e, f, g, h, i, j, k, l)


def fq12_add_fq6(t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n, o, p, q, r = t_m
    return ((a+m) % Q, (b+n) % Q, (c+o) % Q, (d+p) % Q, (e+q) % Q, (f+r) % Q,
            g, h, i, j, k, l)


def fq12_add(t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n, o, p, q, r, s, t, u, v, w, x = t_m
    return ((a+m) % Q, (b+n) % Q, (c+o) % Q, (d+p) % Q, (e+q) % Q, (f+r) % Q,
            (g+s) % Q, (h+t) % Q, (i+u) % Q, (j+v) % Q, (k+w) % Q, (l+x) % Q)


def fq12_sub_fq(t_a, m):
    '''Sub int m from tuple t_a returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    return ((a-m) % Q, b, c, d, e, f, g, h, i, j, k, l)


def fq_sub_fq12(a, t_m):
    '''Sub tuple t_m from int a returning tuple'''
    m, n, o, p, q, r, s, t, u, v, w, x = t_m
    return ((a-m) % Q, -n % Q, -o % Q, -p % Q, -q % Q, -r % Q,
            -s % Q, -t % Q, -u % Q, -v % Q, -w % Q, -x % Q)


def fq12_sub_fq2(t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n = t_m
    return ((a-m) % Q, (b-n) % Q, c, d, e, f, g, h, i, j, k, l)


def fq2_sub_fq12(t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b = t_a
    m, n, o, p, q, r, s, t, u, v, w, x = t_m
    return ((a-m) % Q, (b-n) % Q, -o % Q, -p % Q, -q % Q, -r % Q,
            -s % Q, -t % Q, -u % Q, -v % Q, -w % Q, -x % Q)


def fq12_sub_fq6(t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n, o, p, q, r = t_m
    return ((a-m) % Q, (b-n) % Q, (c-o) % Q, (d-p) % Q, (e-q) % Q, (f-r) % Q,
            g, h, i, j, k, l)


def fq6_sub_fq12(t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f = t_a
    m, n, o, p, q, r, s, t, u, v, w, x = t_m
    return ((a-m) % Q, (b-n) % Q, (c-o) % Q, (d-p) % Q, (e-q) % Q, (f-r) % Q,
            -s % Q, -t % Q, -u % Q, -v % Q, -w % Q, -x % Q)


def fq12_sub(t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n, o, p, q, r, s, t, u, v, w, x = t_m
    return ((a-m) % Q, (b-n) % Q, (c-o) % Q, (d-p) % Q, (e-q) % Q, (f-r) % Q,
            (g-s) % Q, (h-t) % Q, (i-u) % Q, (j-v) % Q, (k-w) % Q, (l-x) % Q)


def fq12_mul_fq(t_a, m):
    '''Multiple tuple t_a on int m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    return (a*m % Q, b*m % Q, c*m % Q, d*m % Q, e*m % Q, f*m % Q,
            g*m % Q, h*m % Q, i*m % Q, j*m % Q, k*m % Q, l*m % Q)


def fq12_mul_fq2(t_a, t_m):
    '''Multiple tuple t_a on tuple t_m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n = t_m
    mul_a = a*m - b*n; mul_b = a*n + b*m
    mul_c = c*m - d*n; mul_d = c*n + d*m
    mul_e = e*m - f*n; mul_f = e*n + f*m
    mul_g = g*m - h*n; mul_h = g*n + h*m
    mul_i = i*m - j*n; mul_j = i*n + j*m
    mul_k = k*m - l*n; mul_l = k*n + l*m
    return (mul_a % Q, mul_b % Q, mul_c % Q, mul_d % Q, mul_e % Q, mul_f % Q,
            mul_g % Q, mul_h % Q, mul_i % Q, mul_j % Q, mul_k % Q, mul_l % Q)


def fq12_mul_fq6(t_a, t_m):
    '''Multiple tuple t_a on tuple t_m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n, o, p, q, r = t_m

    am = a*m; an = a*n; ao = a*o; ap = a*p; aq = a*q; ar = a*r;
    bm = b*m; bn = b*n; bo = b*o; bp = b*p; bq = b*q; br = b*r;
    cm = c*m; cn = c*n; co = c*o; cp = c*p; cq = c*q; cr = c*r;
    dm = d*m; dn = d*n; do = d*o; dp = d*p; dq = d*q; dr = d*r;
    em = e*m; en = e*n; eo = e*o; ep = e*p; eq = e*q; er = e*r;
    fm = f*m; fn = f*n; fo = f*o; fp = f*p; fq = f*q; fr = f*r;
    gm = g*m; gn = g*n; go = g*o; gp = g*p; gq = g*q; gr = g*r;
    hm = h*m; hn = h*n; ho = h*o; hp = h*p; hq = h*q; hr = h*r;
    im = i*m; in_ = i*n; io = i*o; ip = i*p; iq = i*q; ir = i*r;
    jm = j*m; jn = j*n; jo = j*o; jp = j*p; jq = j*q; jr = j*r;
    km = k*m; kn = k*n; ko = k*o; kp = k*p; kq = k*q; kr = k*r;
    lm = l*m; ln = l*n; lo = l*o; lp = l*p; lq = l*q; lr = l*r;

    r0 = am - bn + cq - cr - dr - dq + eo - ep - fp - fo
    r1 = an + bm + cq + cr - dr + dq + eo + ep - fp + fo
    r2 = ao - bp + cm - dn + eq - fr - er - fq
    r3 = ap + bo + cn + dm + eq - fr + er + fq
    r4 = aq - br + co - dp + em - fn
    r5 = ar + bq + cp + do + en + fm
    r6 = gm - hn + iq - ir - jr - jq + ko - kp - lp - lo
    r7 = gn + hm + iq + ir - jr + jq + ko + kp - lp + lo
    r8 = go - hp + im - jn + kq - kr - lr - lq
    r9 = gp + ho + in_ + jm + kq + kr - lr + lq
    r10 = gq - hr + io - jp + km - ln
    r11 = gr + hq + ip + jo + kn + lm
    return (r0 % Q, r1 % Q, r2 % Q, r3 % Q, r4 % Q, r5 % Q,
            r6 % Q, r7 % Q, r8 % Q, r9 % Q, r10 % Q, r11 % Q)


def fq12_mul(t_a, t_m):
    '''Multiple tuple t_a on tuple t_m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n, o, p, q, r, s, t, u, v, w, x = t_m

    am = a*m; an = a*n; ao = a*o; ap = a*p; aq = a*q; ar = a*r;
    as_ = a*s; at = a*t; au = a*u; av = a*v; aw = a*w; ax = a*x;
    bm = b*m; bn = b*n; bo = b*o; bp = b*p; bq = b*q; br = b*r;
    bs = b*s; bt = b*t; bu = b*u; bv = b*v; bw = b*w; bx = b*x;
    cm = c*m; cn = c*n; co = c*o; cp = c*p; cq = c*q; cr = c*r;
    cs = c*s; ct = c*t; cu = c*u; cv = c*v; cw = c*w; cx = c*x;
    dm = d*m; dn = d*n; do = d*o; dp = d*p; dq = d*q; dr = d*r;
    ds = d*s; dt = d*t; du = d*u; dv = d*v; dw = d*w; dx = d*x;
    em = e*m; en = e*n; eo = e*o; ep = e*p; eq = e*q; er = e*r;
    es = e*s; et = e*t; eu = e*u; ev = e*v; ew = e*w; ex = e*x;
    fm = f*m; fn = f*n; fo = f*o; fp = f*p; fq = f*q; fr = f*r;
    fs = f*s; ft = f*t; fu = f*u; fv = f*v; fw = f*w; fx = f*x;
    gm = g*m; gn = g*n; go = g*o; gp = g*p; gq = g*q; gr = g*r;
    gs = g*s; gt = g*t; gu = g*u; gv = g*v; gw = g*w; gx = g*x;
    hm = h*m; hn = h*n; ho = h*o; hp = h*p; hq = h*q; hr = h*r;
    hs = h*s; ht = h*t; hu = h*u; hv = h*v; hw = h*w; hx = h*x;
    im = i*m; in_ = i*n; io = i*o; ip = i*p; iq = i*q; ir = i*r;
    is_ = i*s; it = i*t; iu = i*u; iv = i*v; iw = i*w; ix = i*x;
    jm = j*m; jn = j*n; jo = j*o; jp = j*p; jq = j*q; jr = j*r;
    js = j*s; jt = j*t; ju = j*u; jv = j*v; jw = j*w; jx = j*x;
    km = k*m; kn = k*n; ko = k*o; kp = k*p; kq = k*q; kr = k*r;
    ks = k*s; kt = k*t; ku = k*u; kv = k*v; kw = k*w; kx = k*x;
    lm = l*m; ln = l*n; lo = l*o; lp = l*p; lq = l*q; lr = l*r;
    ls = l*s; lt = l*t; lu = l*u; lv = l*v; lw = l*w; lx = l*x;

    r0 = (am - bn + cq - cr - dr - dq + eo - ep - fp - fo +
          gw - gx - hx - hw + iu - iv - jv - ju + ks - kt - lt - ls)
    r1 = (an + bm + cq + cr - dr + dq + eo + ep - fp + fo +
          gw + gx - hx + hw + iu + iv - jv + ju + ks + kt - lt + ls)
    r2 = (ao - bp + cm - dn + eq - er - fr - fq + gs - ht +
          iw - ix - jx - jw + ku - kv - lv - lu)
    r3 = (ap + bo + cn + dm + eq + er - fr + fq + gt + hs +
          iw + ix - jx + jw + ku + kv - lv + lu)
    r4 = aq - br + co - dp + em - fn + gu - hv + is_ - jt + kw - kx - lx - lw
    r5 = ar + bq + cp + do + en + fm + gv + hu + it + js + kw + kx - lx + lw
    r6 = (as_ - bt + cw - cx - dx - dw + eu - ev - fv - fu +
          gm - hn + iq - ir - jr - jq + ko - kp - lp - lo)
    r7 = (at + bs + cw + cx - dx + dw + eu + ev - fv + fu +
          gn + hm + iq + ir - jr + jq + ko + kp - lp + lo)
    r8 = (au - bv + cs - dt + ew - ex - fx - fw + go - hp +
          im - jn + kq - kr - lr - lq)
    r9 = (av + bu + ct + ds + ew + ex - fx + fw + gp + ho +
          in_ + jm + kq + kr - lr + lq)
    r10 = aw - bx + cu - dv + es - ft + gq - hr + io - jp + km - ln
    r11 = ax + bw + cv + du + et + fs + gr + hq + ip + jo + kn + lm
    return (r0 % Q, r1 % Q, r2 % Q, r3 % Q, r4 % Q, r5 % Q,
            r6 % Q, r7 % Q, r8 % Q, r9 % Q, r10 % Q, r11 % Q)


# ec.py methods
def fq_is_on_curve_affine(px, py, pinf):
    if pinf:
        return True
    left = fq_pow(Q, py, 2)
    right = (fq_pow(Q, px, 3) + bls12381_b) % Q
    return left == right


def fq2_is_on_curve_affine(px, py, pinf):
    if pinf:
        return True
    left = fq2_pow(py, 2)
    right = fq2_add(fq2_pow(px, 3), (bls12381_b, )*2)
    return left == right


def fq12_is_on_curve_affine(px, py, pinf):
    if pinf:
        return True
    left = fq12_pow(py, 2)
    right = fq12_add(fq12_pow(px, 3), (bls12381_b, )*12)
    return left == right


def fq_is_on_curve_jacobian(px, py, pz, pinf):
    xr, yr, infr = fq_to_affine(px, py, pz, pinf)
    return fq_is_on_curve_affine(xr, yr, infr)


def fq2_is_on_curve_jacobian(px, py, pz, pinf):
    xr, yr, infr = fq2_to_affine(px, py, pz, pinf)
    return fq2_is_on_curve_affine(xr, yr, infr)


def fq12_is_on_curve_jacobian(px, py, pz, pinf):
    xr, yr, infr = fq12_to_affine(px, py, pz, pinf)
    return fq12_is_on_curve_affine(xr, yr, infr)


def fq_to_jacobian(px, py, pinf):
    return px, py, 1, pinf


def fq2_to_jacobian(px, py, pinf):
    return px, py, FQ2_ONE_TUPLE, pinf


def fq12_to_jacobian(px, py, pinf):
    return px, py, FQ12_ONE_TUPLE, pinf


def fq_to_affine(px, py, pz, pinf):
    if pinf:
        return 0, 0, pinf
    rx = px * fq_invert(Q, fq_pow(Q, pz, 2)) % Q
    ry = py * fq_invert(Q, fq_pow(Q, pz, 3)) % Q
    return rx, ry, pinf


def fq2_to_affine(px, py, pz, pinf):
    if pinf:
        return FQ2_ZERO_TUPLE, FQ2_ZERO_TUPLE, pinf
    rx = fq2_mul(px, fq2_invert(fq2_pow(pz, 2)))
    ry = fq2_mul(py, fq2_invert(fq2_pow(pz, 3)))
    return rx, ry, pinf


def fq12_to_affine(px, py, pz, pinf):
    if pinf:
        return FQ12_ZERO_TUPLE, FQ12_ZERO_TUPLE, pinf
    rx = fq12_mul(px, fq12_invert(fq12_pow(pz, 2)))
    ry = fq12_mul(py, fq12_invert(fq12_pow(pz, 3)))
    return rx, ry, pinf


def fq_double_point(px, py, pinf):
    left = 3 * fq_pow(Q, px, 2)
    s = left * fq_invert(Q, 2*py)
    xr = (fq_pow(Q, s, 2) - 2*px) % Q
    yr = (s * (px - xr) - py) % Q
    return xr, yr, False


def fq2_double_point(px, py, pinf):
    left = fq2_mul_fq(fq2_pow(px, 2), 3)
    s = fq2_mul(left, fq2_invert(fq2_mul_fq(py, 2)))
    xr = fq2_sub(fq2_pow(s, 2), fq2_mul_fq(px, 2))
    yr = fq2_sub(fq2_mul(s, fq2_sub(px, xr)), py)
    return xr, yr, False


def fq12_double_point(px, py, pinf):
    left = fq12_mul_fq(fq12_pow(px, 2), 3)
    s = fq12_mul(left, fq12_invert(fq12_mul_fq(py, 2)))
    xr = fq12_sub(fq12_pow(s, 2), fq12_mul_fq(px, 2))
    yr = fq12_sub(fq12_mul(s, fq12_sub(px, xr)), py)
    return xr, yr, False


def fq_add_points(x1, y1, inf1, x2, y2, inf2):
    if inf1:
        return x2, y2, inf2
    if inf2:
        return x1, y1, inf1
    if x1 == x2 and y1 == y2:
        return fq_double_point(x1, y1, inf1)
    if x1 == x2:
        return 0, 0, True

    s = (y2 - y1) * fq_invert(Q, x2 - x1)
    xr = (fq_pow(Q, s, 2) - x1 - x2) % Q
    yr = (s * (x1 - xr) - y1) % Q
    return xr, yr, False


def fq2_add_points(x1, y1, inf1, x2, y2, inf2):
    if inf1:
        return x2, y2, inf2
    if inf2:
        return x1, y1, inf1
    if x1 == x2 and y1 == y2:
        return fq2_double_point(x1, y1, inf1)
    if x1 == x2:
        return FQ2_ZERO_TUPLE, FQ2_ZERO_TUPLE, True

    s = fq2_mul(fq2_sub(y2, y1), fq2_invert(fq2_sub(x2, x1)))
    xr = fq2_sub(fq2_sub(fq2_pow(s, 2), x1), x2)
    yr = fq2_sub(fq2_mul(s, fq2_sub(x1, xr)), y1)
    return xr, yr, False


def fq12_add_points(x1, y1, inf1, x2, y2, inf2):
    if inf1:
        return x2, y2, inf2
    if inf2:
        return x1, y1, inf1
    if x1 == x2 and y1 == y2:
        return fq12_double_point(x1, y1, inf1)
    if x1 == x2:
        return FQ12_ZERO_TUPLE, FQ12_ZERO_TUPLE, True

    s = fq12_mul(fq12_sub(y2, y1), fq2_invert(fq12_sub(x2, x1)))
    xr = fq12_sub(fq12_sub(fq12_pow(s, 2), x1), x2)
    yr = fq12_sub(fq12_mul(s, fq12_sub(x1, xr)), y1)
    return xr, yr, False


def fq_scalar_mult_jacobian(c, x1, y1, z1, inf1):
    xr = 1
    yr = 1
    zr = 0
    infr = True
    if inf1 or c % Q == 0:
        return xr, yr, zr, infr
    # addend = p1
    while c > 0:
        if c & 1:
            # result += addend
            xr, yr, zr, infr = fq_add_points_jacobian(xr, yr, zr, infr,
                                                      x1, y1, z1, inf1)
        # double point
        x1, y1, z1 = fq_double_point_jacobian(x1, y1, z1)
        c = c >> 1
    return xr, yr, zr, infr


def fq2_scalar_mult_jacobian(c, x1, y1, z1, inf1):
    xr = FQ2_ONE_TUPLE
    yr = FQ2_ONE_TUPLE
    zr = FQ2_ZERO_TUPLE
    infr = True
    if inf1 or c % Q == 0:
        return xr, yr, zr, infr
    # addend = p1
    while c > 0:
        if c & 1:
            # result += addend
            xr, yr, zr, infr = fq2_add_points_jacobian(xr, yr, zr, infr,
                                                       x1, y1, z1, inf1)
        # double point
        x1, y1, z1 = fq2_double_point_jacobian(x1, y1, z1)
        c = c >> 1
    return xr, yr, zr, infr


def fq12_scalar_mult_jacobian(c, x1, y1, z1, inf1):
    xr = FQ12_ONE_TUPLE
    yr = FQ12_ONE_TUPLE
    zr = FQ12_ZERO_TUPLE
    infr = True
    if inf1 or c % Q == 0:
        return xr, yr, zr, infr
    # addend = p1
    while c > 0:
        if c & 1:
            # result += addend
            xr, yr, zr, infr = fq12_add_points_jacobian(xr, yr, zr, infr,
                                                        x1, y1, z1, inf1)
        # double point
        x1, y1, z1 = fq12_double_point_jacobian(x1, y1, z1)
        c = c >> 1
    return xr, yr, zr, infr


def fq_add_points_jacobian(x1, y1, z1, inf1, x2, y2, z2, inf2):
    if inf1:
        return x2, y2, z2, inf2
    if inf2:
        return x1, y1, z1, inf1
    # u1 = x1*z2^2
    u1 = x1*z2*z2 % Q
    # u2 = x2*z1^2
    u2 = x2*z1*z1 % Q
    # s1 = y1*z2^3
    s1 = y1*z2*z2*z2 % Q
    # s2 = y2*z1^3
    s2 = y2*z1*z1*z1 % Q
    if u1 == u2:
        if s1 != s2:
            xr = yr = 1
            zr = 0
            infr = True
        else:
            xr, yr, zr = fq_double_point_jacobian(x1, y1, z1, Q)
            infr = False
    else:
        # h = u2 - u1
        h = (u2-u1) % Q
        # r = s2 - s1
        r = (s2-s1) % Q
        h_sq = h*h % Q
        h_cu = h*h_sq % Q
        # xr = r^2 - h^3 - 2*u1*h^2
        xr = (r*r % Q - h_cu - 2*u1*h_sq % Q) % Q
        # yr = r*(u1*h^2 - xr) - s1*h^3
        yr = (r*(u1*h_sq % Q - xr) % Q - s1*h_cu % Q) % Q
        # zr = h*z1*z2
        zr = h*z1*z2 % Q
        infr = False
    return xr, yr, zr, infr


def fq2_add_points_jacobian(x1, y1, z1, inf1, x2, y2, z2, inf2):
    if inf1:
        return x2, y2, z2, inf2
    if inf2:
        return x1, y1, z1, inf1
    u1, u2, s1, s2 = fqx_calc_u1_u2_s1_s2(fq2_mul, x1, y1, z1,
                                          x2, y2, z2)
    if u1 == u2:
        if s1 != s2:
            xr = yr = FQ2_ONE_TUPLE
            zr = FQ2_ZERO_TUPLE
            infr = True
        else:
            xr, yr, zr = fq2_double_point_jacobian(x1, y1, z1)
            infr = False
    else:
        func_t = (fq2_mul, fq2_sub, fq2_mul_fq)
        xr, yr, zr = fqx_calc_jp_on_us(func_t, u1, u2, s1, s2, z1, z2)
        infr = False
    return xr, yr, zr, infr


def fq12_add_points_jacobian(x1, y1, z1, inf1, x2, y2, z2, inf2):
    if inf1:
        return x2, y2, z2, inf2
    if inf2:
        return x1, y1, z1, inf1
    u1, u2, s1, s2 = fqx_calc_u1_u2_s1_s2(fq12_mul, x1, y1, z1,
                                          x2, y2, z2)
    if u1 == u2:
        if s1 != s2:
            xr = yr = FQ12_ONE_TUPLE
            zr = FQ12_ZERO_TUPLE
            infr = True
        else:
            xr, yr, zr = fq12_double_point_jacobian(x1, y1, z1)
            infr = False
    else:
        func_t = (fq12_mul, fq12_sub, fq12_mul_fq)
        xr, yr, zr = fqx_calc_jp_on_us(func_t, u1, u2, s1, s2, z1, z2)
        infr = False
    return xr, yr, zr, infr


def fqx_calc_u1_u2_s1_s2(mul_f, x1_t, y1_t, z1_t, x2_t, y2_t, z2_t):
    '''x, y, z inputs of type fq2, returning tuple of fq2 tuples'''
    # u1 = x1*z2^2
    u1 = mul_f(mul_f(x1_t, z2_t), z2_t)
    # u2 = x2*z1^2
    u2 = mul_f(mul_f(x2_t, z1_t), z1_t)
    # s1 = y1*z2^3
    s1 = mul_f(mul_f(mul_f(y1_t, z2_t), z2_t), z2_t)
    # s2 = y2*z1^3
    s2 = mul_f(mul_f(mul_f(y2_t, z1_t), z1_t), z1_t)
    return u1, u2, s1, s2


def fqx_calc_jp_on_us(func_t, u1, u2, s1, s2, z1, z2):
    '''calc jacobian point with tuples u1, u2, s1, s2, z1, z2,
    using func_t functions tuple for operations'''
    mul_f, sub_f, muli_f = func_t
    # h = u2 - u1
    h = sub_f(u2, u1)
    # r = s2 - s1
    r = sub_f(s2, s1)
    h_sq = mul_f(h, h)
    h_cu = mul_f(h, h_sq)
    # x3 = r^2 - h^3 - 2*u1*h^2
    x3 = sub_f(sub_f(mul_f(r, r), h_cu),
               mul_f(h_sq, muli_f(u1, 2)))
    # y3 = r*(u1*h^2 - x3) - s1*h^3
    y3 = sub_f(mul_f(r, sub_f(mul_f(u1, h_sq), x3)),
               mul_f(s1, h_cu))
    # z3 = h*z1*z2
    z3 = mul_f(mul_f(z1, z2), h)
    return x3, y3, z3


def fq_double_point_jacobian(X, Y, Z):
    '''dobule point with fq int X, Y, Z, returning tuple'''
    # S = 4*X*Y^2
    S = 4*X*Y*Y % Q

    Z_sq = Z*Z % Q
    Y_sq = Y*Y % Q
    Y_4th = Y_sq*Y_sq % Q

    # M = 3*X^2 + a*Z^4
    # A is 0 for bls12-381
    M = 3*X*X % Q

    # X' = M^2 - 2*S
    X_p = (M*M % Q - 2*S % Q) % Q
    # Y' = M*(S - X') - 8*Y^4
    Y_p = (M*((S - X_p) % Q) % Q - 8*Y_4th % Q) % Q
    # Z' = 2*Y*Z
    Z_p = 2*Y*Z % Q
    return X_p, Y_p, Z_p


def fq2_double_point_jacobian(X, Y, Z):
    '''dobule point with fq2 tuples X, Y, Z, returning tuple of tuples'''
    func_t = (fq2_mul, fq2_mul_fq, fq2_add, fq2_sub)
    return fqx_double_point_jacobian(func_t, X, Y, Z)


def fq12_double_point_jacobian(X, Y, Z):
    '''dobule point with fq12 tuples X, Y, Z, returning tuple of tuples'''
    func_t = (fq12_mul, fq12_mul_fq, fq12_add, fq12_sub)
    return fqx_double_point_jacobian(func_t, X, Y, Z)


def fqx_double_point_jacobian(func_t, X, Y, Z):
    '''dobule point with tuples X, Y, Z, returning tuple of tuples,
    using func_t functions tuple for operations'''
    mul_f, mul_i_f, add_f, sub_f = func_t
    # S = 4*X*Y^2
    S = mul_f(mul_f(mul_i_f(X, 4), Y), Y)

    Z_sq = mul_f(Z, Z)
    Y_sq = mul_f(Y, Y)
    Y_4th = mul_f(Y_sq, Y_sq)

    # M = 3*X^2 + a*Z^4
    # A is 0 for bls12-381
    M = mul_f(mul_i_f(X, 3), X)

    # X' = M^2 - 2*S
    X_p = sub_f(mul_f(M, M), mul_i_f(S, 2))
    # Y' = M*(S - X') - 8*Y^4
    Y_p = sub_f(mul_f(M, sub_f(S, X_p)), mul_i_f(Y_4th, 8))
    # Z' = 2*Y*Z
    Z_p = mul_f(mul_i_f(Y, 2), Z)
    return X_p, Y_p, Z_p


def fq2_untwist(x_t, y_t):
    m, n = x_t
    new_x = (0, 0, 0, 0, (tw1*m - tw2*n) % Q, (tw1*n + tw2*m) % Q,
             0, 0, 0, 0, 0, 0)
    m, n = y_t
    new_y = (0, 0, 0, 0, 0, 0,
             0, 0, (tw1*m - tw2*n) % Q, (tw1*n + tw2*m) % Q, 0, 0)
    return new_x, new_y


def fq12_untwist(x_t, y_t):
    m, n, o, p, q, r, s, t, u, v, w, x = x_t
    e = tw1
    f = tw2
    eo = e*o; ep = e*p; fp = f*p; fo = f*o;
    eq = e*q; er = e*r; fr = f*r; fq = f*q;
    em = e*m; fn = f*n; en = e*n; fm = f*m;
    eu = e*u; ev = e*v; fv = f*v; fu = f*u;
    ew = e*w; ex = e*x; fx = f*x; fw = f*w;
    es = e*s; ft = f*t; et = e*t; fs = f*s;
    #0  eo - ep - fp - fo
    #1  eo + ep - fp + fo
    #2  eq - er - fr - fq
    #3  eq + er - fr + fq
    #4  em - fn
    #5  en + fm
    #6  eu - ev - fv - fu
    #7  eu + ev - fv + fu
    #8  ew - ex - fx - fw
    #9  ew + ex - fx + fw
    #10 es - ft
    #11 et + fs
    new_x = ((eo - ep - fp - fo) % Q,
             (eo + ep - fp + fo) % Q,
             (eq - er - fr - fq) % Q,
             (eq + er - fr + fq) % Q,
             (em - fn) % Q,
             (en + fm) % Q,
             (eu - ev - fv - fu) % Q,
             (eu + ev - fv + fu) % Q,
             (ew - ex - fx - fw) % Q,
             (ew + ex - fx + fw) % Q,
             (es - ft) % Q,
             (et + fs) % Q)
    m, n, o, p, q, r, s, t, u, v, w, x = y_t
    i = tw1
    j = tw2
    iu = i*u; iv = i*v; jv = j*v; ju = j*u;
    iw = i*w; ix = i*x; jx = j*x; jw = j*w;
    is_ = i*s; jt = j*t; it = i*t; js = j*s;
    iq = i*q; ir = i*r; jr = j*r; jq = j*q;
    im = i*m; jn = j*n; in_ = i*n; jm = j*m;
    io = i*o; jp = j*p; ip = i*p; jo = j*o;
    #0  iu - iv - jv - ju
    #1  iu + iv - jv + ju
    #2  iw - ix - jx - jw
    #3  iw + ix - jx + jw
    #4  is - jt
    #5  it + js
    #6  iq - ir - jr - jq
    #7  iq + ir - jr + jq
    #8  im - jn
    #9  in + jm
    #10 io - jp
    #11 ip + jo
    new_y = ((iu - iv - jv - ju) % Q,
             (iu + iv - jv + ju) % Q,
             (iw - ix - jx - jw) % Q,
             (iw + ix - jx + jw) % Q,
             (is_ - jt) % Q,
             (it + js) % Q,
             (iq - ir - jr - jq) % Q,
             (iq + ir - jr + jq) % Q,
             (im - jn) % Q,
             (in_ + jm) % Q,
             (io - jp) % Q,
             (ip + jo) % Q)
    return new_x, new_y


def fq2_twist(x_t, y_t):
    m, n = x_t
    new_x = (0, 0, m, n, 0, 0, 0, 0, 0, 0, 0, 0)
    m, n = y_t
    new_y = (0, 0, 0, 0, 0, 0, 0, 0, m, n, 0, 0)
    return new_x, new_y


def fq12_twist(x_t, y_t):
    m, n, o, p, q, r, s, t, u, v, w, x = x_t
    new_x = ((q - r) % Q, (q + r) % Q, m, n, o, p,
             (w - x) % Q, (w + x) % Q, s, t, u, v)
    m, n, o, p, q, r, s, t, u, v, w, x = y_t
    new_y = ((u - v) % Q, (u + v) % Q, (w - x) % Q, (w + x) % Q, s, t,
             (q - r) % Q, (q + r) % Q, m, n, o, p)
    return new_x, new_y


# pairing.py methods
def fq2_double_line_eval(rx_t, ry_t, px, py):
    #R12 = untwist(R)
    r12_x, r12_y = fq2_untwist(rx_t, ry_t)

    #slope = (3 * pow(R12.x, 2)) / (2 * R12.y)
    slope = fq12_mul_fq(fq12_mul(r12_x, r12_x), 3)
    slope = fq12_mul(slope, fq12_invert(fq12_mul_fq(r12_y, 2)))

    #v = R12.y - slope * R12.x
    v_t = fq12_sub(r12_y, fq12_mul(slope, r12_x))

    #res = P.y - P.x * slope - v
    res_t = fq_sub_fq12(py, fq12_mul_fq(slope, px))
    res_t = fq12_sub(res_t, v_t)
    return res_t


def fq2_add_line_eval(rx_t, ry_t, qx_t, qy_t, px, py):
    #R12 = untwist(R)
    #Q12 = untwist(Q)
    r12x_t, r12y_t = fq2_untwist(rx_t, ry_t)
    q12x_t, q12y_t = fq2_untwist(qx_t, qy_t)

    # This is the case of a vertical line, where the denominator
    # will be 0.
    #if R12 == Q12.negate():
    #    return P.x - R12.x
    nq12x_t = fq12_neg(q12x_t)
    nq12y_t = fq12_neg(q12y_t)
    if r12x_t == nq12x_t and r12y_t == nq12y_t:
        return fq_sub_fq12(px, r12x_t)

    #slope = (Q12.y - R12.y) / (Q12.x - R12.x)
    slope = fq12_mul(fq12_sub(q12y_t, r12y_t),
                          fq12_invert(fq12_sub(q12x_t, r12x_t)))
    #v = (Q12.y * R12.x - R12.y * Q12.x) / (R12.x - Q12.x)
    v_t = fq12_mul(fq12_sub(fq12_mul(q12y_t, r12x_t),
                                     fq12_mul(r12y_t, q12x_t)),
                        fq12_invert(fq12_sub(r12x_t, q12x_t)))

    #res = P.y - P.x * slope - v
    res_t = fq_sub_fq12(py, fq12_mul_fq(slope, px))
    res_t = fq12_sub(res_t, v_t)
    return res_t


def int_to_bits(i):
    if i < 1:
        return [0]
    bits = []
    while i != 0:
        bits.append(i % 2)
        i = i // 2
    return list(reversed(bits))


def fq_miller_loop(px, py, pinf, qx_t, qy_t, qinf):
    T = bls12381_nx
    T_bits = int_to_bits(T)
    rx_t = qx_t
    ry_t = qy_t
    rinf = qinf
    f = FQ12_ONE_TUPLE
    for i in range(1, len(T_bits)):
        # Compute sloped line lrr
        llr = fq2_double_line_eval(rx_t, ry_t, px, py)
        f = fq12_mul(fq12_pow(f, 2), llr)
        # R = 2 * R
        rx_t, ry_t, rinf = fq2_double_point(rx_t, ry_t, rinf)
        if T_bits[i] == 1:
            # Compute sloped line lrq
            lrq = fq2_add_line_eval(rx_t, ry_t, qx_t, qy_t, px, py)
            f = fq12_mul(f, lrq)
            # R = R + Q
            rx_t, ry_t, rinf = fq2_add_points(rx_t, ry_t, rinf,
                                              qx_t, qy_t, qinf)
    return f


def fq_ate_pairing_multi(Ps, Qs):
    prod = FQ12_ONE_TUPLE
    for i in range(len(Qs)):
        px, py, pinf = Ps[i]
        qx, qy, qinf = Qs[i]
        ml_res = fq_miller_loop(px, py, pinf, qx, qy, qinf)
        prod = fq12_mul(prod, ml_res)
    return fq12_final_exp(prod)


def fq12_final_exp(t_x):
    ans = fq12_pow(t_x, FINAL_EXP_E)
    ans = fq12_mul(fq12_qi_pow(ans, 2), ans)
    ans = fq12_floordiv(fq12_qi_pow(ans, 6), ans)
    return ans


# Frobenius coefficients for raising elements to q**i -th powers
# These are specific to this given bls12381_q
frob_coeffs = {
    (2, 1, 1): -1 % bls12381_q,
    (6, 1, 1): (0,
                int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                    '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                    'fd8bfd00000000aaac', 16)),
    (6, 1, 2): (int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                    '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                    'fd8bfd00000000aaad', 16), 0),
    (6, 2, 1): (int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                    '3be6f89688de17d813620a00022e01fffffffeff'
                    'fe', 16), 0),
    (6, 2, 2): (int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                    '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                    'fd8bfd00000000aaac', 16), 0),
    (6, 3, 1): (0, 1),
    (6, 3, 2): (int('0x1a0111ea397fe69a4b1ba7b6434bacd764774b'
                    '84f38512bf6730d2a0f6b0f6241eabfffeb153ff'
                    'ffb9feffffffffaaaa', 16), 0),
    (6, 4, 1): (int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                    '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                    'fd8bfd00000000aaac', 16), 0),
    (6, 4, 2): (int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                    '3be6f89688de17d813620a00022e01fffffffeff'
                    'fe', 16), 0),
    (6, 5, 1): (0,
                int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                    '3be6f89688de17d813620a00022e01fffffffeff'
                    'fe', 16)),
    (6, 5, 2): (int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                    '3be6f89688de17d813620a00022e01fffffffeff'
                    'ff', 16), 0),
    (12, 1, 1): (int('0x1904d3bf02bb0667c231beb4202c0d1f0fd603'
                     'fd3cbd5f4f7b2443d784bab9c4f67ea53d63e781'
                     '3d8d0775ed92235fb8', 16),
                 int('0xfc3e2b36c4e03288e9e902231f9fb854a14787'
                     'b6c7b36fec0c8ec971f63c5f282d5ac14d6c7ec2'
                     '2cf78a126ddc4af3', 16), 0, 0, 0, 0),
    (12, 2, 1): (int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                     '3be6f89688de17d813620a00022e01fffffffeff'
                     'ff', 16), 0, 0, 0, 0, 0),
    (12, 3, 1): (int('0x135203e60180a68ee2e9c448d77a2cd91c3ded'
                     'd930b1cf60ef396489f61eb45e304466cf3e67fa'
                     '0af1ee7b04121bdea2', 16),
                 int('0x6af0e0437ff400b6831e36d6bd17ffe48395da'
                     'bc2d3435e77f76e17009241c5ee67992f72ec05f'
                     '4c81084fbede3cc09', 16), 0, 0, 0, 0),
    (12, 4, 1): (int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                     '3be6f89688de17d813620a00022e01fffffffeff'
                     'fe', 16), 0, 0, 0, 0, 0),
    (12, 5, 1): (int('0x144e4211384586c16bd3ad4afa99cc9170df35'
                     '60e77982d0db45f3536814f0bd5871c1908bd478'
                     'cd1ee605167ff82995', 16),
                 int('0x5b2cfd9013a5fd8df47fa6b48b1e045f398162'
                     '40c0b8fee8beadf4d8e9c0566c63a3e6e257f873'
                     '29b18fae980078116', 16), 0, 0, 0, 0),
    (12, 6, 1): (int('0x1a0111ea397fe69a4b1ba7b6434bacd764774b'
                     '84f38512bf6730d2a0f6b0f6241eabfffeb153ff'
                     'ffb9feffffffffaaaa', 16), 0, 0, 0, 0, 0),
    (12, 7, 1): (int('0xfc3e2b36c4e03288e9e902231f9fb854a14787'
                     'b6c7b36fec0c8ec971f63c5f282d5ac14d6c7ec2'
                     '2cf78a126ddc4af3', 16),
                 int('0x1904d3bf02bb0667c231beb4202c0d1f0fd603'
                     'fd3cbd5f4f7b2443d784bab9c4f67ea53d63e781'
                     '3d8d0775ed92235fb8', 16), 0, 0, 0, 0),
    (12, 8, 1): (int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                     '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                     'fd8bfd00000000aaac', 16), 0, 0, 0, 0, 0),
    (12, 9, 1): (int('0x6af0e0437ff400b6831e36d6bd17ffe48395da'
                     'bc2d3435e77f76e17009241c5ee67992f72ec05f'
                     '4c81084fbede3cc09', 16),
                 int('0x135203e60180a68ee2e9c448d77a2cd91c3ded'
                     'd930b1cf60ef396489f61eb45e304466cf3e67fa'
                     '0af1ee7b04121bdea2', 16), 0, 0, 0, 0),
    (12, 10, 1): (int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                      '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                      'fd8bfd00000000aaad', 16), 0, 0, 0, 0, 0),
    (12, 11, 1): (int('0x5b2cfd9013a5fd8df47fa6b48b1e045f398162'
                      '40c0b8fee8beadf4d8e9c0566c63a3e6e257f873'
                      '29b18fae980078116', 16),
                  int('0x144e4211384586c16bd3ad4afa99cc9170df35'
                      '60e77982d0db45f3536814f0bd5871c1908bd478'
                      'cd1ee605167ff82995', 16), 0, 0, 0, 0),
}

try:
    from .fields_t_c import (
        fq_invert,
        fq_floordiv,
        fq_pow,
    )
    from .fields_t_c import (
        fq2_invert,
        fq2_floordiv,
        fq2_qi_pow,
        fq2_pow,
    )
    from .fields_t_c import (
        fq6_invert,
        fq6_floordiv,
        fq6_qi_pow,
        fq6_add,
        fq6_mul,
    )
    from .fields_t_c import (
        fq12_invert,
        fq12_floordiv,
        fq12_qi_pow,
        fq12_pow,
        fq12_mul_fq,
        fq12_add,
        fq12_mul,
    )
    from .fields_t_c import (
        fq2_to_affine,
        fq2_double_point,
        fq2_add_points,
        fq_double_point_jacobian,
        fq2_double_point_jacobian,
        fq12_double_point_jacobian,
        fq_add_points_jacobian,
        fq2_add_points_jacobian,
        fq12_add_points_jacobian,
        fq2_scalar_mult_jacobian,
        fq2_double_line_eval,
        fq2_add_line_eval,
        fq2_untwist,
        fq_miller_loop,
        fq12_final_exp,
        fq_ate_pairing_multi,
    )
except ImportError as e:
    logging.error(f'Can not import from fields_t_c: {e}')
