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

bls12381_q = int('0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf'
                 '6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab', 16)


FQ2_ROOT = -1 % bls12381_q
FQ2_ONE_TUPLE = (1, 0)
FQ2_ZERO_TUPLE = (0, 0)
FQ6_ROOT_TUPLE = (1, 1)
FQ6_ONE_TUPLE = (1, 0, 0, 0, 0, 0)
FQ6_ZERO_TUPLE = (0, 0, 0, 0, 0, 0)
FQ12_ROOT_TUPLE = (0, 0, 1, 0, 0, 0)
FQ12_ONE_TUPLE = (1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
FQ12_ZERO_TUPLE = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


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


def fq2_neg(P, t_a):
    '''Neg tuple t_a returning tuple'''
    a, b = t_a
    return (-a % P, -b % P)


def fq2_invert(P, t_a):
    '''Invert tuple t_a returning tuple'''
    a, b = t_a
    factor = fq_invert(P, a * a + b * b)
    return ((a*factor) % P, (-b*factor) % P)


def fq2_pow(P, t_a, e):
    '''Pow tuple t_a returning tuple'''
    m, n = t_a
    a, b = 1, 0
    fq2r = FQ2_ROOT
    while e:
        if e & 1:
            a, b = (a*m + b*n*fq2r) % P, (a*n + b*m) % P
        m, n = (m*m + n*n*fq2r) % P, (m*n + n*m) % P
        e >>= 1
    return (a, b)


def fq2_qi_pow(P, t_x, i):
    '''Calc qi_power on t_x tuple returning tuple'''
    global bls12381_q, frob_coeffs
    if P != bls12381_q:
        raise NotImplementedError
    i %= 2
    if i == 0:
        return t_x
    return (t_x[0], t_x[1]*frob_coeffs[2, i, 1] % P)


def fq2_mul_by_nonresidue(P, t_a):
    '''Mul by nonresidue on tuple t_a returning tuple'''
    a, b = t_a
    return ((a-b) % P, (a+b) % P)


def fq2_add_fq(P, t_a, m):
    '''Add tuple t_a and int m returning tuple'''
    a, b = t_a
    return ((a+m) % P, b)


def fq2_add_fq2(P, t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b = t_a
    m, n = t_m
    return ((a+m) % P, (b+n) % P)


def fq_sub_fq2(P, a, t_m):
    '''Sub tuple t_m from int a returning tuple'''
    m, n = t_m
    return ((a-m) % P, -n % P)


def fq2_sub_fq(P, t_a, m):
    '''Sub int m from tuple t_a returning tuple'''
    a, b = t_a
    return ((a-m) % P, b)


def fq2_sub_fq2(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b = t_a
    m, n = t_m
    return ((a-m) % P, (b-n) % P)


def fq2_mul_fq(P, t_a, m):
    '''Multiple tuple t_a on int m returning tuple'''
    a, b = t_a
    return (a*m % P, b*m % P)


def fq2_mul_fq2(P, t_a, t_m):
    '''Multiple tuple t_a on tuple t_m returning tuple'''
    a, b = t_a
    m, n = t_m
    return ((a*m + b*n*FQ2_ROOT) % P, (a*n + b*m) % P)


def fq6_neg(P, t_a):
    '''Neg tuple t_a returning tuple'''
    a, b, c, d, e, f = t_a
    return (-a % P, -b % P, -c % P, -d % P, -e % P, -f % P)


def fq6_invert(P, t_x):
    '''Invert tuple t_a returning tuple'''
    a, b, c = t_x[:2], t_x[2:4], t_x[4:]
    g0 = fq2_mul_fq2(P, a, a)
    g0 = fq2_sub_fq2(P, g0, fq2_mul_fq2(P, b, fq2_mul_by_nonresidue(P, c)))
    g1 = fq2_sub_fq2(P, fq2_mul_by_nonresidue(P, fq2_mul_fq2(P, c, c)),
                     fq2_mul_fq2(P, a, b))
    g2 = fq2_sub_fq2(P, fq2_mul_fq2(P, b, b), fq2_mul_fq2(P, a, c))
    g0a = fq2_mul_fq2(P, g0, a)
    g1cpg2b = fq2_add_fq2(P, fq2_mul_fq2(P, g1, c), fq2_mul_fq2(P, g2, b))
    factor = fq2_invert(P, fq2_add_fq2(P, g0a,
                        fq2_mul_by_nonresidue(P, g1cpg2b)))
    ar, br = fq2_mul_fq2(P, g0, factor)
    cr, dr = fq2_mul_fq2(P, g1, factor)
    er, fr = fq2_mul_fq2(P, g2, factor)
    return (ar % P, br % P, cr % P, dr % P, er % P, fr % P)


def fq6_pow(P, t_a, e):
    '''Pow tuple t_a returning tuple'''
    t_ans = FQ6_ONE_TUPLE
    while e:
        if e & 1:
            t_ans = fq6_mul_fq6(P, t_ans, t_a)
        t_a = fq6_mul_fq6(P, t_a, t_a)
        e >>= 1
    a, b, c, d, e, f = t_ans
    return (a % P, b % P, c % P, d % P, e % P, f % P)


def fq6_qi_pow(P, t_x, i):
    '''Calc qi_power on t_x tuple returning tuple'''
    global bls12381_q, frob_coeffs
    if P != bls12381_q:
        raise NotImplementedError
    i %= 6
    if i == 0:
        return t_x
    a, b = fq2_qi_pow(P, t_x[:2], i)
    c, d = fq2_mul_fq2(P, fq2_qi_pow(P, t_x[2:4], i), frob_coeffs[6, i, 1])
    e, f = fq2_mul_fq2(P, fq2_qi_pow(P, t_x[4:6], i), frob_coeffs[6, i, 2])
    return (a, b, c, d, e, f)


def fq6_mul_by_nonresidue(P, t_x):
    '''Mul by nonresidue on tuple t_a returning tuple'''
    ar, br = fq2_mul_fq2(P, t_x[4:], FQ6_ROOT_TUPLE)
    cr, dr = t_x[:2]
    er, fr = t_x[2:4]
    return (ar, br, cr, dr, er, fr)


def fq6_add_fq(P, t_a, m):
    '''Add tuple t_a and int m returning tuple'''
    a, b, c, d, e, f = t_a
    return ((a+m) % P, b, c, d, e, f)


def fq6_add_fq2(P, t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b, c, d, e, f = t_a
    m, n = t_m
    return ((a+m) % P, (b+n) % P, c, d, e, f)


def fq6_add_fq6(P, t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b, c, d, e, f = t_a
    m, n, o, p, q, r = t_m
    return ((a+m) % P, (b+n) % P, (c+o) % P,
            (d+p) % P, (e+q) % P, (f+r) % P)


def fq_sub_fq6(P, a, t_m):
    '''Sub tuple t_m from int a returning tuple'''
    m, n, o, p, q, r = t_m
    return ((a-m) % P, -n % P, -o % P, -p % P, -q % P, -r % P)


def fq6_sub_fq(P, t_a, m):
    '''Sub int m from tuple t_a returning tuple'''
    a, b, c, d, e, f = t_a
    return ((a-m) % P, b, c, d, e, f)


def fq2_sub_fq6(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b = t_a
    m, n, o, p, q, r = t_m
    return ((a-m) % P, (b-n) % P, -o % P, -p % P, -q % P, -r % P)


def fq6_sub_fq2(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f = t_a
    m, n = t_m
    return ((a-m) % P, (b-n) % P, c, d, e, f)


def fq6_sub_fq6(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f = t_a
    m, n, o, p, q, r = t_m
    return ((a-m) % P, (b-n) % P, (c-o) % P,
            (d-p) % P, (e-q) % P, (f-r) % P)


def fq6_mul_fq(P, t_a, m):
    '''Multiple tuple t_a on int m returning tuple'''
    a, b, c, d, e, f = t_a
    return (a*m % P, b*m % P, c*m % P, d*m % P, e*m % P, f*m % P)


def fq6_mul_fq2(P, t_a, t_m):
    '''Multiple tuple t_a on tuple t_m returning tuple'''
    a, b, c, d, e, f = t_a
    m, n = t_m
    fq2r = FQ2_ROOT
    return ((a*m + b*n*fq2r) % P, (a*n + b*m) % P,
            (c*m + d*n*fq2r) % P, (c*n + d*m) % P,
            (e*m + f*n*fq2r) % P, (e*n + f*m) % P)


def fq6_mul_fq6(P, t_a, t_m):
    '''Multiple tuple t_a on tuple t_m returning tuple'''
    a, b, c, d, e, f = t_a
    m, n, o, p, q, r = t_m
    fq2r = FQ2_ROOT
    dnfq2r = d*n*fq2r; drfq2r = d*r*fq2r
    fpfq2r = f*p*fq2r; frfq2r = f*r*fq2r
    cq = c*q; cr = c*r; dq = d*q; eo = e*o; ep = e*p; fo = f*o
    cm = c*m; cn = c*n; dm = d*m; eq = e*q; er = e*r; fq = f*q
    mul_a = (a*m + b*n*fq2r + cq + drfq2r + (cr + dq)*fq2r +
             eo + fpfq2r + (ep + fo)*fq2r)
    mul_b = (a*n + b*m + cq + drfq2r + cr + dq + eo + fpfq2r + ep + fo)
    mul_c = (a*o + b*p*fq2r + cm + dnfq2r + eq + frfq2r + (er + fq)*fq2r)
    mul_d = (a*p + b*o + cn + dm + eq + frfq2r + er + fq)
    mul_e = (a*q + b*r*fq2r + c*o + d*p*fq2r + e*m + f*n*fq2r)
    mul_f = (a*r + b*q + c*p + d*o + e*n + f*m)
    return (mul_a % P, mul_b % P, mul_c % P, mul_d % P, mul_e % P, mul_f % P)


def fq12_neg(P, t_a):
    '''Neg tuple t_a returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    return (-a % P, -b % P, -c % P, -d % P, -e % P, -f % P,
            -g % P, -h % P, -i % P, -j % P, -k % P, -l % P)


def fq12_invert(P, t_x):
    '''Invert tuple t_a returning tuple'''
    a, b = t_x[:6], t_x[6:12]
    aa = fq6_mul_fq6(P, a, a)
    bb = fq6_mul_fq6(P, b, b)
    factor = fq6_invert(P, fq6_sub_fq6(P, aa, fq6_mul_by_nonresidue(P, bb)))
    ar, br, cr, dr, er, fr = fq6_mul_fq6(P, a, factor)
    gr, hr, ir, jr, kr, lr = fq6_mul_fq6(P, fq6_neg(P, b), factor)
    return (ar % P, br % P, cr % P, dr % P, er % P, fr % P,
            gr % P, hr % P, ir % P, jr % P, kr % P, lr % P)


def fq12_pow(P, t_a, e):
    '''Pow tuple t_a returning tuple'''
    t_ans = FQ12_ONE_TUPLE
    while e:
        if e & 1:
            t_ans = fq12_mul_fq12(P, t_ans, t_a)
        t_a = fq12_mul_fq12(P, t_a, t_a)
        e >>= 1
    return t_ans


def fq12_qi_pow(P, t_x, i):
    '''Calc qi_power on t_x tuple returning tuple'''
    global bls12381_q, frob_coeffs
    if P != bls12381_q:
        raise NotImplementedError
    i %= 12
    if i == 0:
        return t_x
    a, b, c, d, e, f = fq6_qi_pow(P, t_x[:6], i)
    g, h, i, j, k, l = fq6_mul_fq6(P, fq6_qi_pow(P, t_x[6:12], i),
                                   frob_coeffs[12, i, 1])
    return (a, b, c, d, e, f, g, h, i, j, k, l)


def fq12_add_fq(P, t_a, m):
    '''Add tuple t_a and int m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    return ((a+m) % P, b, c, d, e, f, g, h, i, j, k, l)


def fq12_add_fq2(P, t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n = t_m
    return ((a+m) % P, (b+n) % P, c, d, e, f, g, h, i, j, k, l)


def fq12_add_fq6(P, t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n, o, p, q, r = t_m
    return ((a+m) % P, (b+n) % P, (c+o) % P, (d+p) % P, (e+q) % P, (f+r) % P,
            g, h, i, j, k, l)


def fq12_add_fq12(P, t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n, o, p, q, r, s, t, u, v, w, x = t_m
    return ((a+m) % P, (b+n) % P, (c+o) % P, (d+p) % P, (e+q) % P, (f+r) % P,
            (g+s) % P, (h+t) % P, (i+u) % P, (j+v) % P, (k+w) % P, (l+x) % P)


def fq12_sub_fq(P, t_a, m):
    '''Sub int m from tuple t_a returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    return ((a-m) % P, b, c, d, e, f, g, h, i, j, k, l)


def fq_sub_fq12(P, a, t_m):
    '''Sub tuple t_m from int a returning tuple'''
    m, n, o, p, q, r, s, t, u, v, w, x = t_m
    return ((a-m) % P, -n % P, -o % P, -p % P, -q % P, -r % P,
            -s % P, -t % P, -u % P, -v % P, -w % P, -x % P)


def fq12_sub_fq2(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n = t_m
    return ((a-m) % P, (b-n) % P, c, d, e, f, g, h, i, j, k, l)


def fq2_sub_fq12(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b = t_a
    m, n, o, p, q, r, s, t, u, v, w, x = t_m
    return ((a-m) % P, (b-n) % P, -o % P, -p % P, -q % P, -r % P,
            -s % P, -t % P, -u % P, -v % P, -w % P, -x % P)


def fq12_sub_fq6(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n, o, p, q, r = t_m
    return ((a-m) % P, (b-n) % P, (c-o) % P, (d-p) % P, (e-q) % P, (f-r) % P,
            g, h, i, j, k, l)


def fq6_sub_fq12(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f = t_a
    m, n, o, p, q, r, s, t, u, v, w, x = t_m
    return ((a-m) % P, (b-n) % P, (c-o) % P, (d-p) % P, (e-q) % P, (f-r) % P,
            -s % P, -t % P, -u % P, -v % P, -w % P, -x % P)


def fq12_sub_fq12(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n, o, p, q, r, s, t, u, v, w, x = t_m
    return ((a-m) % P, (b-n) % P, (c-o) % P, (d-p) % P, (e-q) % P, (f-r) % P,
            (g-s) % P, (h-t) % P, (i-u) % P, (j-v) % P, (k-w) % P, (l-x) % P)


def fq12_mul_fq(P, t_a, m):
    '''Multiple tuple t_a on int m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    return (a*m % P, b*m % P, c*m % P, d*m % P, e*m % P, f*m % P,
            g*m % P, h*m % P, i*m % P, j*m % P, k*m % P, l*m % P)


def fq12_mul_fq2(P, t_a, t_m):
    '''Multiple tuple t_a on tuple t_m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n = t_m
    fq2r = FQ2_ROOT
    mul_a = a*m + b*n*fq2r; mul_b = a*n + b*m
    mul_c = c*m + d*n*fq2r; mul_d = c*n + d*m
    mul_e = e*m + f*n*fq2r; mul_f = e*n + f*m
    mul_g = g*m + h*n*fq2r; mul_h = g*n + h*m
    mul_i = i*m + j*n*fq2r; mul_j = i*n + j*m
    mul_k = k*m + l*n*fq2r; mul_l = k*n + l*m
    return (mul_a % P, mul_b % P, mul_c % P, mul_d % P, mul_e % P, mul_f % P,
            mul_g % P, mul_h % P, mul_i % P, mul_j % P, mul_k % P, mul_l % P)


def fq12_mul_fq6(P, t_a, t_m):
    '''Multiple tuple t_a on tuple t_m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n, o, p, q, r = t_m
    fq2r = FQ2_ROOT
    drfq2r = d*r*fq2r; fpfq2r = f*p*fq2r; frfq2r = f*r*fq2r
    jrfq2r = j*r*fq2r; lpfq2r = l*p*fq2r; lrfq2r = l*r*fq2r
    cq = c*q; cr = c*r; dq = d*q; eo = e*o; ep = e*p
    fo = f*o; eq = e*q; er = e*r; fq = f*q
    iq = i*q; ir = i*r; jq = j*q; ko = k*o; kp = k*p
    lo = l*o; kq = k*q; kr = k*r; lq = l*q
    mul_a = (a*m + b*n*fq2r + cq + drfq2r + (cr + dq)*fq2r + eo +
             fpfq2r + (ep + fo)*fq2r)
    mul_b = a*n + b*m + cq + drfq2r + cr + dq + eo + fpfq2r + ep + fo
    mul_c = (a*o + b*p*fq2r + c*m + d*n*fq2r + eq + frfq2r +
             (er + fq)*fq2r)
    mul_d = a*p + b*o + c*n + d*m + eq + frfq2r + er + fq
    mul_e = a*q + b*r*fq2r + c*o + d*p*fq2r + e*m + f*n*fq2r
    mul_f = a*r + b*q + c*p + d*o + e*n + f*m
    mul_g = (g*m + h*n*fq2r + iq + jrfq2r + (ir + jq)*fq2r + ko +
             lpfq2r + (kp + lo)*fq2r)
    mul_h = g*n + h*m + iq + jrfq2r + ir + jq + ko + lpfq2r + kp + lo
    mul_i = (g*o + h*p*fq2r + i*m + j*n*fq2r + kq + lrfq2r +
             (kr + lq)*fq2r)
    mul_j = g*p + h*o + i*n + j*m + kq + lrfq2r + kr + lq
    mul_k = g*q + h*r*fq2r + i*o + j*p*fq2r + k*m + l*n*fq2r
    mul_l = g*r + h*q + i*p + j*o + k*n + l*m
    return (mul_a % P, mul_b % P, mul_c % P, mul_d % P, mul_e % P, mul_f % P,
            mul_g % P, mul_h % P, mul_i % P, mul_j % P, mul_k % P, mul_l % P)


def fq12_mul_fq12(P, t_a, t_m):
    '''Multiple tuple t_a on tuple t_m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n, o, p, q, r, s, t, u, v, w, x = t_m
    fq2r = FQ2_ROOT
    cq = c*q; cr = c*r; dq = d*q; drfq2r = d*r*fq2r
    eo = e*o; ep = e*p; eq = e*q; er = e*r; fo = f*o
    fpfq2r = f*p*fq2r; fq = f*q; frfq2r = f*r*fq2r
    iw = i*w; ix = i*x; jw = j*w; jxfq2r = j*x*fq2r
    ku = k*u; kv = k*v; kw = k*w; kx = k*x; lu = l*u
    lvfq2r = l*v*fq2r; lw = l*w; lxfq2r = l*x*fq2r
    cw = c*w; cx = c*x; dw = d*w; dxfq2r = d*x*fq2r
    eu = e*u; ev = e*v; ew = e*w; ex = e*x; fu = f*u
    fvfq2r = f*v*fq2r; fw = f*w; fxfq2r = f*x*fq2r
    iq = i*q; ir = i*r; jq = j*q; jrfq2r = j*r*fq2r
    ko = k*o; kp = k*p; kq = k*q; kr = k*r; lo = l*o
    lpfq2r = l*p*fq2r; lq = l*q; lrfq2r = l*r*fq2r
    mul_a1 = (a*m + b*n*fq2r + cq + drfq2r + (cr + dq)*fq2r + eo +
              fpfq2r + (ep + fo)*fq2r)
    mul_a2 = (g*w + h*x*fq2r + i*u + j*v*fq2r + k*s + l*t*fq2r +
              (g*x + h*w + i*v + j*u + k*t + l*s)*fq2r)
    mul_b1 = a*n + b*m + cq + drfq2r + cr + dq + eo + fpfq2r + ep + fo
    mul_b2 = (g*w + h*x*fq2r + i*u + j*v*fq2r + k*s + l*t*fq2r + g*x +
              h*w + i*v + j*u + k*t + l*s)
    mul_c1 = (a*o + b*p*fq2r + c*m + d*n*fq2r + eq + frfq2r +
              (er + fq)*fq2r)
    mul_c2 = (g*s + h*t*fq2r + iw + jxfq2r + (ix + jw)*fq2r + ku +
              lvfq2r + (kv + lu)*fq2r)
    mul_d1 = a*p + b*o + c*n + d*m + eq + frfq2r + er + fq
    mul_d2 = g*t + h*s + iw + jxfq2r + ix + jw + ku + lvfq2r + kv + lu
    mul_e1 = a*q + b*r*fq2r + c*o + d*p*fq2r + e*m + f*n*fq2r
    mul_e2 = (g*u + h*v*fq2r + i*s + j*t*fq2r + kw + lxfq2r +
              (kx + lw)*fq2r)
    mul_f1 = a*r + b*q + c*p + d*o + e*n + f*m
    mul_f2 = g*v + h*u + i*t + j*s + kw + lxfq2r + kx + lw
    mul_g1 = (a*s + b*t*fq2r + cw + dxfq2r + (cx + dw)*fq2r + eu +
              fvfq2r + (ev + fu)*fq2r)
    mul_g2 = (g*m + h*n*fq2r + iq + jrfq2r + (ir + jq)*fq2r + ko +
              lpfq2r + (kp + lo)*fq2r)
    mul_h1 = a*t + b*s + cw + dxfq2r + cx + dw + eu + fvfq2r + ev + fu
    mul_h2 = g*n + h*m + iq + jrfq2r + ir + jq + ko + lpfq2r + kp + lo
    mul_i1 = (a*u + b*v*fq2r + c*s + d*t*fq2r + ew + fxfq2r +
              (ex + fw)*fq2r)
    mul_i2 = (g*o + h*p*fq2r + i*m + j*n*fq2r + kq + lrfq2r +
              (kr + lq)*fq2r)
    mul_j1 = a*v + b*u + c*t + d*s + ew + fxfq2r + ex + fw
    mul_j2 = g*p + h*o + i*n + j*m + kq + lrfq2r + kr + lq
    mul_k1 = a*w + b*x*fq2r + c*u + d*v*fq2r + e*s + f*t*fq2r
    mul_k2 = g*q + h*r*fq2r + i*o + j*p*fq2r + k*m + l*n*fq2r
    mul_l1 = a*x + b*w + c*v + d*u + e*t + f*s
    mul_l2 = g*r + h*q + i*p + j*o + k*n + l*m
    return ((mul_a1 + mul_a2) % P, (mul_b1 + mul_b2) % P,
            (mul_c1 + mul_c2) % P, (mul_d1 + mul_d2) % P,
            (mul_e1 + mul_e2) % P, (mul_f1 + mul_f2) % P,
            (mul_g1 + mul_g2) % P, (mul_h1 + mul_h2) % P,
            (mul_i1 + mul_i2) % P, (mul_j1 + mul_j2) % P,
            (mul_k1 + mul_k2) % P, (mul_l1 + mul_l2) % P)


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
