from copy import deepcopy


bls12381_q = int('0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf'
                 '6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab', 16)


FQ2_ROOT = -1%bls12381_q
FQ2_ONE_TUPLE = (1, 0)
FQ2_ZERO_TUPLE = (0, 0)
FQ6_ROOT_TUPLE = (1, 1)
FQ6_ONE_TUPLE = (1, 0, 0, 0, 0, 0)
FQ6_ZERO_TUPLE = (0, 0, 0, 0, 0, 0)
FQ12_ROOT_TUPLE = (0, 0, 1, 0, 0, 0)
FQ12_ONE_TUPLE = (1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
FQ12_ZERO_TUPLE = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


def fq_int_invert(P, a):
    '''Ivnert int value using extended euclidian algorithm for inversion'''
    p = P
    x0, x1, y0, y1 = 1, 0, 0, 1
    while p != 0:
        q, a, p = a // p, p, a % p
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return x0%P


def fq_int_pow(P, a, X):
    '''Pow a to X mod P '''
    if X == 0:
        return 1%P
    res = 1
    while X > 0:
        if X & 1:
            res = res * a % P
        X >>=  1
        a = a * a % P
    return res


def fq_int_floordiv(P, a, X):
    return a * fq_int_invert(P, X) % P


def fq2_t_neg(P, t_a):
    '''Neg tuple t_a returning tuple'''
    a, b = t_a
    return (-a%P, -b%P)


def fq2_t_invert(P, t_a):
    '''Invert tuple t_a returning tuple'''
    a, b = t_a
    factor = fq_int_invert(P, a * a + b * b)
    return ((a*factor)%P, (-b*factor)%P)


def fq2_t_pow(P, t_a, e):
    '''Pow tuple t_a returning tuple'''
    m, n = t_a
    a, b = 1, 0
    fq2r = FQ2_ROOT
    while e:
        if e & 1:
            a, b = (a*m + b*n*fq2r)%P, (a*n + b*m)%P
        m, n = (m*m + n*n*fq2r)%P, (m*n + n*m)%P
        e >>= 1
    return (a, b)


def fq2_t_qi_pow(P, t_x, i):
    '''Calc qi_power on t_x tuple returning tuple'''
    global bls12381_q, frob_coeffs
    if P != bls12381_q:
        raise NotImplementedError
    i %= 2
    if i == 0:
        return t_x
    return (t_x[0], t_x[1]*frob_coeffs[2, i, 1].Z%P)


def fq2_t_mul_by_nonresidue(P, t_a):
    '''Mul by nonresidue on tuple t_a returning tuple'''
    a, b = t_a
    return ((a-b)%P, (a+b)%P)


def fq2_t_add_fq_int(P, t_a, m):
    '''Add tuple t_a and int m returning tuple'''
    a, b = t_a
    return ((a+m)%P, b)


def fq2_t_add_fq2_t(P, t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b = t_a
    m, n = t_m
    return ((a+m)%P, (b+n)%P)


def fq_int_sub_fq2_t(P, a, t_m):
    '''Sub tuple t_m from int a returning tuple'''
    m, n = t_m
    return ((a-m)%P, -n%P)


def fq2_t_sub_fq_int(P, t_a, m):
    '''Sub int m from tuple t_a returning tuple'''
    a, b = t_a
    return ((a-m)%P, b)


def fq2_t_sub_fq2_t(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b = t_a
    m, n = t_m
    return ((a-m)%P, (b-n)%P)


def fq2_t_mul_fq_int(P, t_a, m):
    '''Multiple tuple t_a on int m returning tuple'''
    a, b = t_a
    return (a*m%P, b*m%P)


def fq2_t_mul_fq2_t(P, t_a, t_m):
    '''Multiple tuple t_a on tuple t_m returning tuple'''
    a, b = t_a
    m, n = t_m
    return ((a*m + b*n*FQ2_ROOT)%P, (a*n + b*m)%P)


def fq6_t_neg(P, t_a):
    '''Neg tuple t_a returning tuple'''
    a, b, c, d, e, f = t_a
    return (-a%P, -b%P, -c%P, -d%P, -e%P, -f%P)


def fq6_t_invert(P, t_x):
    '''Invert tuple t_a returning tuple'''
    a, b, c = t_x[:2], t_x[2:4], t_x[4:]
    g0 = fq2_t_mul_fq2_t(P, a, a)
    g0 = fq2_t_sub_fq2_t(P, g0,
                         fq2_t_mul_fq2_t(P, b,
                                         fq2_t_mul_by_nonresidue(P, c)))
    g1 = fq2_t_sub_fq2_t(P,
                         fq2_t_mul_by_nonresidue(P,
                                                 fq2_t_mul_fq2_t(P, c, c)),
                         fq2_t_mul_fq2_t(P, a, b))
    g2 = fq2_t_sub_fq2_t(P,
                         fq2_t_mul_fq2_t(P, b, b),
                         fq2_t_mul_fq2_t(P, a, c))
    g0a = fq2_t_mul_fq2_t(P, g0, a)
    g1cpg2b = fq2_t_add_fq2_t(P,
                              fq2_t_mul_fq2_t(P, g1, c),
                              fq2_t_mul_fq2_t(P, g2, b))
    factor = fq2_t_invert(P,
                          fq2_t_add_fq2_t(P,
                                          g0a,
                                          fq2_t_mul_by_nonresidue(P, g1cpg2b)))
    ar, br = fq2_t_mul_fq2_t(P, g0, factor)
    cr, dr = fq2_t_mul_fq2_t(P, g1, factor)
    er, fr = fq2_t_mul_fq2_t(P, g2, factor)
    return (ar%P, br%P, cr%P, dr%P, er%P, fr%P)


def fq6_t_pow(P, t_a, e):
    '''Pow tuple t_a returning tuple'''
    t_ans = FQ6_ONE_TUPLE
    while e:
        if e & 1:
            t_ans = fq6_t_mul_fq6_t(P, t_ans, t_a)
        t_a = fq6_t_mul_fq6_t(P, t_a, t_a)
        e >>= 1
    a, b, c, d, e, f = t_ans
    return (a%P, b%P, c%P, d%P, e%P, f%P)


def fq6_t_qi_pow(P, t_x, i):
    '''Calc qi_power on t_x tuple returning tuple'''
    global bls12381_q, frob_coeffs
    if P != bls12381_q:
        raise NotImplementedError
    i %= 6
    if i == 0:
        return t_x
    a, b = fq2_t_qi_pow(P, t_x[:2], i)
    c, d = fq2_t_mul_fq2_t(P, fq2_t_qi_pow(P, t_x[2:4], i),
                           frob_coeffs[6, i, 1].ZT)
    e, f = fq2_t_mul_fq2_t(P, fq2_t_qi_pow(P, t_x[4:6], i),
                           frob_coeffs[6, i, 2].ZT)
    return (a, b, c, d, e, f)


def fq6_t_mul_by_nonresidue(P, t_x):
    '''Mul by nonresidue on tuple t_a returning tuple'''
    ar, br = fq2_t_mul_fq2_t(P, t_x[4:], FQ6_ROOT_TUPLE)
    cr, dr = t_x[:2]
    er, fr = t_x[2:4]
    return (ar, br, cr, dr, er, fr)


def fq6_t_add_fq_int(P, t_a, m):
    '''Add tuple t_a and int m returning tuple'''
    a, b, c, d, e, f = t_a
    return ((a+m)%P, b, c, d, e, f)


def fq6_t_add_fq2_t(P, t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b, c, d, e, f = t_a
    m, n = t_m
    return ((a+m)%P, (b+n)%P, c, d, e, f)


def fq6_t_add_fq6_t(P, t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b, c, d, e, f = t_a
    m, n, o, p, q, r = t_m
    return ((a+m)%P, (b+n)%P, (c+o)%P,
            (d+p)%P, (e+q)%P, (f+r)%P)


def fq_int_sub_fq6_t(P, a, t_m):
    '''Sub tuple t_m from int a returning tuple'''
    m, n, o, p, q, r = t_m
    return ((a-m)%P, -n%P, -o%P, -p%P, -q%P, -r%P)


def fq6_t_sub_fq_int(P, t_a, m):
    '''Sub int m from tuple t_a returning tuple'''
    a, b, c, d, e, f = t_a
    return ((a-m)%P, b, c, d, e, f)


def fq2_t_sub_fq6_t(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b = t_a
    m, n, o, p, q, r = t_m
    return ((a-m)%P, (b-n)%P, -o%P, -p%P, -q%P, -r%P)


def fq6_t_sub_fq2_t(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f = t_a
    m, n = t_m
    return ((a-m)%P, (b-n)%P, c, d, e, f)


def fq6_t_sub_fq6_t(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f = t_a
    m, n, o, p, q, r = t_m
    return ((a-m)%P, (b-n)%P, (c-o)%P,
            (d-p)%P, (e-q)%P, (f-r)%P)


def fq6_t_mul_fq_int(P, t_a, m):
    '''Multiple tuple t_a on int m returning tuple'''
    a, b, c, d, e, f = t_a
    return (a*m%P, b*m%P, c*m%P, d*m%P, e*m%P, f*m%P)


def fq6_t_mul_fq2_t(P, t_a, t_m):
    '''Multiple tuple t_a on tuple t_m returning tuple'''
    a, b, c, d, e, f = t_a
    m, n = t_m
    fq2r = FQ2_ROOT
    return ((a*m + b*n*fq2r)%P, (a*n + b*m)%P,
            (c*m + d*n*fq2r)%P, (c*n + d*m)%P,
            (e*m + f*n*fq2r)%P, (e*n + f*m)%P)


def fq6_t_mul_fq6_t(P, t_a, t_m):
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
    return (mul_a%P, mul_b%P, mul_c%P, mul_d%P, mul_e%P, mul_f%P)


def fq12_t_neg(P, t_a):
    '''Neg tuple t_a returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    return (-a%P, -b%P, -c%P, -d%P, -e%P, -f%P,
            -g%P, -h%P, -i%P, -j%P, -k%P, -l%P)


def fq12_t_invert(P, t_x):
    '''Invert tuple t_a returning tuple'''
    a, b = t_x[:6], t_x[6:12]
    aa = fq6_t_mul_fq6_t(P, a, a)
    bb = fq6_t_mul_fq6_t(P, b, b)
    factor = fq6_t_invert(P,
                          fq6_t_sub_fq6_t(P,
                                          aa,
                                          fq6_t_mul_by_nonresidue(P, bb)))
    ar, br, cr, dr, er, fr = fq6_t_mul_fq6_t(P, a, factor)
    gr, hr, ir, jr, kr, lr = fq6_t_mul_fq6_t(P,
                                             fq6_t_neg(P, b),
                                             factor)
    return (ar%P, br%P, cr%P, dr%P, er%P, fr%P,
            gr%P, hr%P, ir%P, jr%P, kr%P, lr%P)


def fq12_t_pow(P, t_a, e):
    '''Pow tuple t_a returning tuple'''
    t_ans = FQ12_ONE_TUPLE
    while e:
        if e & 1:
            t_ans = fq12_t_mul_fq12_t(P, t_ans, t_a)
        t_a = fq12_t_mul_fq12_t(P, t_a, t_a)
        e >>= 1
    return t_ans


def fq12_t_qi_pow(P, t_x, i):
    '''Calc qi_power on t_x tuple returning tuple'''
    global bls12381_q, frob_coeffs
    if P != bls12381_q:
        raise NotImplementedError
    i %= 12
    if i == 0:
        return t_x
    a, b, c, d, e, f = fq6_t_qi_pow(P, t_x[:6], i)
    g, h, i, j, k, l = fq6_t_mul_fq6_t(P,
                                       fq6_t_qi_pow(P, t_x[6:12], i),
                                       frob_coeffs[12, i, 1].ZT)
    return (a, b, c, d, e, f, g, h, i, j, k, l)


def fq12_t_add_fq_int(P, t_a, m):
    '''Add tuple t_a and int m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    return ((a+m)%P, b, c, d, e, f, g, h, i, j, k, l)


def fq12_t_add_fq2_t(P, t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n = t_m
    return ((a+m)%P, (b+n)%P, c, d, e, f, g, h, i, j, k, l)


def fq12_t_add_fq6_t(P, t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n, o, p, q, r = t_m
    return ((a+m)%P, (b+n)%P, (c+o)%P, (d+p)%P, (e+q)%P, (f+r)%P,
            g, h, i, j, k, l)


def fq12_t_add_fq12_t(P, t_a, t_m):
    '''Add tuple t_a and tuple t_m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n, o, p, q, r, s, t, u, v, w, x = t_m
    return ((a+m)%P, (b+n)%P, (c+o)%P, (d+p)%P, (e+q)%P, (f+r)%P,
            (g+s)%P, (h+t)%P, (i+u)%P, (j+v)%P, (k+w)%P, (l+x)%P)


def fq12_t_sub_fq_int(P, t_a, m):
    '''Sub int m from tuple t_a returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    return ((a-m)%P, b, c, d, e, f, g, h, i, j, k, l)


def fq_int_sub_fq12_t(P, a, t_m):
    '''Sub tuple t_m from int a returning tuple'''
    m, n, o, p, q, r, s, t, u, v, w, x = t_m
    return ((a-m)%P, -n%P, -o%P, -p%P, -q%P, -r%P,
            -s%P, -t%P, -u%P, -v%P, -w%P, -x%P)


def fq12_t_sub_fq2_t(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n = t_m
    return ((a-m)%P, (b-n)%P, c, d, e, f, g, h, i, j, k, l)


def fq2_t_sub_fq12_t(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b = t_a
    m, n, o, p, q, r, s, t, u, v, w, x = t_m
    return ((a-m)%P, (b-n)%P, -o%P, -p%P, -q%P, -r%P,
            -s%P, -t%P, -u%P, -v%P, -w%P, -x%P)


def fq12_t_sub_fq6_t(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n, o, p, q, r = t_m
    return ((a-m)%P, (b-n)%P, (c-o)%P, (d-p)%P, (e-q)%P, (f-r)%P,
            g, h, i, j, k, l)


def fq6_t_sub_fq12_t(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f = t_a
    m, n, o, p, q, r, s, t, u, v, w, x = t_m
    return ((a-m)%P, (b-n)%P, (c-o)%P, (d-p)%P, (e-q)%P, (f-r)%P,
            -s%P, -t%P, -u%P, -v%P, -w%P, -x%P)


def fq12_t_sub_fq12_t(P, t_a, t_m):
    '''Sub tuple t_m from tuple t_a returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    m, n, o, p, q, r, s, t, u, v, w, x = t_m
    return ((a-m)%P, (b-n)%P, (c-o)%P, (d-p)%P, (e-q)%P, (f-r)%P,
            (g-s)%P, (h-t)%P, (i-u)%P, (j-v)%P, (k-w)%P, (l-x)%P)


def fq12_t_mul_fq_int(P, t_a, m):
    '''Multiple tuple t_a on int m returning tuple'''
    a, b, c, d, e, f, g, h, i, j, k, l = t_a
    return (a*m%P, b*m%P, c*m%P, d*m%P, e*m%P, f*m%P,
            g*m%P, h*m%P, i*m%P, j*m%P, k*m%P, l*m%P)


def fq12_t_mul_fq2_t(P, t_a, t_m):
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
    return (mul_a%P, mul_b%P, mul_c%P, mul_d%P, mul_e%P, mul_f%P,
            mul_g%P, mul_h%P, mul_i%P, mul_j%P, mul_k%P, mul_l%P)


def fq12_t_mul_fq6_t(P, t_a, t_m):
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
    return (mul_a%P, mul_b%P, mul_c%P, mul_d%P, mul_e%P, mul_f%P,
            mul_g%P, mul_h%P, mul_i%P, mul_j%P, mul_k%P, mul_l%P)


def fq12_t_mul_fq12_t(P, t_a, t_m):
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
    return ((mul_a1 + mul_a2)%P, (mul_b1 + mul_b2)%P,
            (mul_c1 + mul_c2)%P, (mul_d1 + mul_d2)%P,
            (mul_e1 + mul_e2)%P, (mul_f1 + mul_f2)%P,
            (mul_g1 + mul_g2)%P, (mul_h1 + mul_h2)%P,
            (mul_i1 + mul_i2)%P, (mul_j1 + mul_j2)%P,
            (mul_k1 + mul_k2)%P, (mul_l1 + mul_l2)%P)


class Fq:
    ''' Represents an element of a finite field mod a prime q.'''
    extension = 1

    @staticmethod
    def zero(Q):
        if Q > 0:
            global FQ_ZERO
            return FQ_ZERO
        else:
            return Fq(Q, 0)

    @staticmethod
    def one(Q):
        if Q > 1:
            global FQ_ONE
            return FQ_ONE
        else:
            return Fq(Q, 1)

    @classmethod
    def from_fq(cls, Q, fq):
        return fq

    def __init__(self, Q, X):
        tx = type(X)
        if tx == int:
            Z = X % Q
        elif tx == Fq:
            Z = X.Z
        else:
            raise TypeError('Fq must be consturcted from Fq or int')
        self.Z = Z
        self.Q = Q

    def __str__(self):
        s = hex(self.Z)
        s = '%s..%s' % (s[0:7], s[-5:]) if len(s) > 10 else s
        return 'Fq(Q, %s)' % s

    def __repr__(self):
        return 'Fq(Q, %s)' % hex(self.Z)

    def __iter__(self):
        return self

    def __int__(self):
        return self.Z

    def __deepcopy__(self, memo):
        return Fq(self.Q, self.Z)

    def serialize(self):
        return self.Z.to_bytes(48, 'big')

    def qi_power(self, i):
        return self

    def __mod__(self, X):
        tx = type(X)
        if tx == int:
            return self.Z % X
        else:
            return NotImplemented

    def __neg__(self):
        return Fq(self.Q, -self.Z)

    def __invert__(self):
        return Fq(self.Q, fq_int_invert(self.Q, self.Z))

    def __pow__(self, X):
        return Fq(self.Q, fq_int_pow(self.Q, self.Z, X))

    def __add__(self, X):
        tx = type(X)
        if tx == Fq:
            return Fq(self.Q, self.Z + X.Z)
        elif tx == int:
            return Fq(self.Q, self.Z + X)
        else:
            return NotImplemented

    def __radd__(self, X):
        global T_FQ_INT
        tx = type(X)
        if tx in T_FQ_INT:
            return self + X
        else:
            return NotImplemented

    def __sub__(self, X):
        tx = type(X)
        if tx == Fq:
            return Fq(self.Q, self.Z - X.Z)
        elif tx == int:
            return Fq(self.Q, self.Z - X)
        else:
            return NotImplemented

    def __rsub__(self, X):
        tx = type(X)
        if tx == Fq:
            return Fq(self.Q, X.Z - self.Z)
        elif tx == int:
            return Fq(self.Q, X - self.Z)
        else:
            return NotImplemented

    def __mul__(self, X):
        tx = type(X)
        if tx == Fq:
            return Fq(self.Q, self.Z * X.Z)
        elif tx == int:
            return Fq(self.Q, self.Z * X)
        else:
            return NotImplemented

    def __rmul__(self, X):
        global T_FQ_INT
        tx = type(X)
        if tx in T_FQ_INT:
            return self * X
        else:
            return NotImplemented

    def __floordiv__(self, X):
        tx = type(X)
        if tx == Fq:
            return Fq(self.Q, fq_int_floordiv(self.Q, self.Z, X.Z))
        elif tx == int:
            return Fq(self.Q, fq_int_floordiv(self.Q, self.Z, X))
        else:
            return NotImplemented

    __truediv__ = __floordiv__

    def __eq__(self, X):
        tx = type(X)
        if tx == Fq:
            return self.Z == X.Z and self.Q == X.Q
        elif tx == int:
            return self.Z == X
        else:
            return NotImplemented

    def __lt__(self, X):
        tx = type(X)
        if tx == Fq:
            return self.Z < X.Z
        elif tx == int:
            return self.Z < X
        else:
            return NotImplemented

    def __gt__(self, X):
        tx = type(X)
        if tx == Fq:
            return self.Z > X.Z
        elif tx == int:
            return self.Z > X
        else:
            return NotImplemented

    def modsqrt(self):
        if self.Z == 0:
            return self
        if pow(self.Z, (self.Q - 1) // 2, self.Q) != 1:
            raise ValueError('No sqrt exists')
        if self.Q % 4 == 3:
            return Fq(self.Q, pow(self.Z, (self.Q + 1) // 4, self.Q))
        if self.Q % 8 == 5:
            return Fq(self.Q, pow(self.Z, (self.Q + 3) // 8, self.Q))

        # p % 8 == 1. Tonelli Shanks algorithm for finding square root
        S = 0
        q = self.Q - 1

        while q % 2 == 0:
            q = q // 2
            S += 1

        z = 0
        for i in range(self.Q):
            euler = pow(i, (self.Q - 1)//2, self.Q)
            if euler == -1 % self.Q:
                z = i
                break

        M = S
        c = pow(z, q, self.Q)
        t = pow(self.Z, q, self.Q)
        R = pow(self.Z, (q + 1) // 2, self.Q)

        while True:
            if t == 0:
                return Fq(self.Q, 0)
            if t == 1:
                return Fq(self.Q, R)
            i = 0
            f = t
            while f != 1:
                f = pow(f, 2, self.Q)
                i += 1
            b = pow(c, pow(2, M - i - 1, self.Q), self.Q)
            M = i
            c = pow(b, 2, self.Q)
            t = (t * c) % self.Q
            R = (R * b) % self.Q
    

class FieldExtBase:
    ''' Represents an extension of a field (or extension of an extension).
    The elements of the tuple can be other FieldExtBase or they can be
    Fq elements. For example, Fq2 = (Fq, Fq). Fq12 = (Fq6, Fq6), etc.'''
    @classmethod
    def from_fq(cls, Q, fq):
        ret = super().__new__(cls)
        ret.Q = Q
        ret.ZT = tuple(0 if i else fq.Z for i in range(cls.extension))
        return ret

    def __init__(self, Q):
        global bls12381_q
        if Q != bls12381_q:
            raise TypeError('Q must be bls12381 q')
        self.Q = Q

    def __bool__(self):
        return any(z for z in self.ZT)

    def __deepcopy__(self, memo):
        cls = type(self)
        ret = super().__new__(cls)
        ret.Q = self.Q
        ret.ZT = self.ZT[:]
        return ret

    def serialize(self):
        # Returns the concatenated coordinates in big endian bytes
        sum_bytes = bytes([])
        for z in self.ZT:
            sum_bytes += z.to_bytes(48, 'big')
        return sum_bytes

    def __radd__(self, other):
        return self.__add__(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __eq__(self, other):
        global T_FE_BASED, T_FQ_INT
        tx = type(other)
        if tx == type(self):
            return self.ZT == other.ZT
        elif tx in T_FE_BASED:
            if self.extension < other.extension:
                return NotImplemented
            if sum(self.ZT[other.extension:]) != 0:
                return False
            for i in range(other.extension):
                if self.ZT[i] != other.ZT[i]:
                    return False
            return True
        elif tx in T_FQ_INT:
            if sum(self.ZT[1:]) != 0:
                return False
            if tx == Fq:
                return self.ZT[0] == other.Z
            else:
                return self.ZT[0] == other
        else:
            return NotImplemented

    def __neq__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        # Reverse the order for comparison (3i + 1 > 2i + 7)
        return self.ZT[::-1].__lt__(other.ZT[::-1])

    def __gt__(self, other):
        return self.ZT[::-1].__gt__(other.ZT[::-1])


class Fq2(FieldExtBase):
    # Fq2 is constructed as Fq(u) / (u2 - β) where β = -1
    extension = 2
    embedding = 2
    basefield = Fq
    root = Fq(bls12381_q, -1)

    @staticmethod
    def zero(Q):
        global FQ2_ZERO
        return FQ2_ZERO

    @staticmethod
    def one(Q):
        global FQ2_ONE
        return FQ2_ONE

    def __init__(self, Q, *args):
        super(Fq2, self).__init__(Q)
        if len(args) == 1 and type(args[0]) == tuple:
            self.ZT = args[0]
            return
        if len(args) != 2:
            raise Exception('Invalid number of arguments')
        a, b = args
        tx = type(a)
        if tx == int:
            a %= Q
        elif tx == Fq:
            a = a.Z
        else:
            raise Exception('Invalid arguments type, must be Fq or int')
        tx = type(b)
        if tx == int:
            b %= Q
        elif tx == Fq:
            b = b.Z
        else:
            raise Exception('Invalid arguments type, must be Fq or int')
        self.ZT = (a, b)

    def __iter__(self):
        for z in self.ZT:
            yield Fq(self.Q, z)

    def __getitem__(self, idx):
        return Fq(self.Q, self.ZT[idx])

    def __str__(self):
        return ('Fq2(Q, %s)' % ', '.join(str(fq) for fq in self))

    def __repr__(self):
        return ('Fq2(Q, %s)' % ', '.join(repr(fq) for fq in self))

    def __neg__(self):
       return Fq2(self.Q, fq2_t_neg(self.Q, self.ZT))

    def __invert__(self):
        return Fq2(self.Q, fq2_t_invert(self.Q, self.ZT))

    def __pow__(self, e):
        return Fq2(self.Q, fq2_t_pow(self.Q, self.ZT, e))

    def qi_power(self, i):
        global bls12381_q, frob_coeffs
        if self.Q != bls12381_q:
            raise NotImplementedError
        i %= 2
        if i == 0:
            return self
        return Fq2(self.Q, fq2_t_qi_pow(self.Q, self.ZT, i))

    def mul_by_nonresidue(self):
        # multiply by u + 1
        return Fq2(self.Q, fq2_t_mul_by_nonresidue(self.Q, self.ZT))

    def __add__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq2:
            return Fq2(Q, fq2_t_add_fq2_t(Q, self.ZT, other.ZT))
        elif tx == int:
            return Fq2(Q, fq2_t_add_fq_int(Q, self.ZT, other))
        elif tx == Fq:
            return Fq2(Q, fq2_t_add_fq_int(Q, self.ZT, other.Z))
        else:
            return NotImplemented

    def __sub__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq2:
            return Fq2(Q, fq2_t_sub_fq2_t(Q, self.ZT, other.ZT))
        elif tx == int:
            return Fq2(Q, fq2_t_sub_fq_int(Q, self.ZT, other))
        elif tx == Fq:
            return Fq2(Q, fq2_t_sub_fq_int(Q, self.ZT, other.Z))
        else:
            return NotImplemented

    def __rsub__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq2:
            return Fq2(Q, fq2_t_sub_fq2_t(Q, other.ZT, self.ZT))
        elif tx == int:
            return Fq2(Q, fq_int_sub_fq2_t(Q, other, self.ZT))
        elif tx == Fq:
            return Fq2(Q, fq_int_sub_fq2_t(Q, other.Z, self.ZT))
        else:
            return NotImplemented

    def __mul__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq2:
            return Fq2(Q, fq2_t_mul_fq2_t(Q, self.ZT, other.ZT))
        elif tx == int:
            return Fq2(Q, fq2_t_mul_fq_int(Q, self.ZT, other))
        elif tx == Fq:
            return Fq2(Q, fq2_t_mul_fq_int(Q, self.ZT, other.Z))
        else:
            return NotImplemented

    def __floordiv__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq12:
            return Fq12(Q, fq12_t_mul_fq2_t(Q, fq12_t_invert(Q, other.ZT),
                                            self.ZT))
        elif tx == Fq6:
            return Fq6(Q, fq6_t_mul_fq2_t(Q, fq2_t_invert(Q, other.ZT),
                                          self.ZT))
        elif tx == Fq2:
            return Fq2(Q, fq2_t_mul_fq2_t(Q, self.ZT,
                                          fq2_t_invert(Q, other.ZT)))
        elif tx == int:
            return Fq2(Q, fq2_t_mul_fq_int(Q, self.ZT,
                                           fq_int_invert(Q, other)))
        elif tx == Fq:
            return Fq2(Q, fq2_t_mul_fq_int(Q, self.ZT,
                                           fq_int_invert(Q, other.Z)))
        else:
            return NotImplemented

    __truediv__ = __floordiv__

    def modsqrt(self):
        '''Using algorithm 8 (complex method) for square roots in
        https://eprint.iacr.org/2012/685.pdf
        This is necessary for computing y value given an x value.'''
        a0, a1 = self
        if a1 == Fq.zero(self.Q):
            return a0.modsqrt()
        alpha = pow(a0, 2) + pow(a1, 2)
        gamma = pow(alpha, (self.Q - 1)//2)
        if (gamma == Fq(self.Q, -1)):
            raise ValueError('No sqrt exists')
        alpha = alpha.modsqrt()
        delta = (a0 + alpha) * ~Fq(self.Q, 2)
        gamma = pow(delta, (self.Q - 1)//2)
        if (gamma == Fq(self.Q, -1)):
            delta = (a0 - alpha) * ~Fq(self.Q, 2)

        x0 = delta.modsqrt()
        x1 = a1 * ~(2*x0)
        return Fq2(self.Q, x0, x1)


class Fq6(FieldExtBase):
    # Fq6 is constructed as Fq2(v) / (v3 - ξ) where ξ = u + 1
    extension = 6
    embedding = 3
    basefield = Fq2
    root = Fq2(bls12381_q, FQ6_ROOT_TUPLE)

    @staticmethod
    def zero(Q):
        global FQ6_ZERO
        return FQ6_ZERO

    @staticmethod
    def one(Q):
        global FQ6_ONE
        return FQ6_ONE

    def __init__(self, Q, *args):
        super(Fq6, self).__init__(Q)
        if len(args) == 1 and type(args[0]) == tuple:
            self.ZT = args[0][:]
            return
        if len(args) != 3:
            raise Exception('Invalid number of arguments')
        a, b, c = args
        if type(a) != Fq2 or type(b) != Fq2 or type(c) !=  Fq2:
            raise Exception('Invalid arguments type, must be Fq2')
        self.ZT = a.ZT + b.ZT + c.ZT

    def __iter__(self):
        for i in range(0, 6, 2):
            yield Fq2(self.Q, self.ZT[i:i+2])

    def __getitem__(self, idx):
        return Fq2(self.Q, self.ZT[idx*2:idx*2+2])

    def __str__(self):
        return ('Fq6(Q, %s)' % ', '.join(str(fq2) for fq2 in self))

    def __repr__(self):
        return ('Fq6(Q, %s)' % ', '.join(repr(fq2) for fq2 in self))

    def __neg__(self):
       return Fq6(self.Q, fq6_t_neg(self.Q, self.ZT))

    def __invert__(self):
        return Fq6(self.Q, fq6_t_invert(self.Q, self.ZT))

    def __pow__(self, e):
        return Fq6(self.Q, fq6_t_pow(self.Q, self.ZT, e))

    def qi_power(self, i):
        global bls12381_q, frob_coeffs
        if self.Q != bls12381_q:
            raise NotImplementedError
        i %= 6
        if i == 0:
            return self
        return Fq6(self.Q, fq6_t_qi_pow(self.Q, self.ZT, i))

    def mul_by_nonresidue(self):
        # multiply by v
        return Fq6(self.Q, fq6_t_mul_by_nonresidue(self.Q, self.ZT))

    def __add__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq6:
            return Fq6(Q, fq6_t_add_fq6_t(Q, self.ZT, other.ZT))
        elif tx == Fq2:
            return Fq6(Q, fq6_t_add_fq2_t(Q, self.ZT, other.ZT))
        elif tx == Fq:
            return Fq6(Q, fq6_t_add_fq_int(Q, self.ZT, other.Z))
        elif tx == int:
            return Fq6(Q, fq6_t_add_fq_int(Q, self.ZT, other))
        else:
            return NotImplemented

    def __sub__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq6:
            return Fq6(Q, fq6_t_sub_fq6_t(Q, self.ZT, other.ZT))
        elif tx == Fq2:
            return Fq6(Q, fq6_t_sub_fq2_t(Q, self.ZT, other.ZT))
        elif tx == Fq:
            return Fq6(Q, fq6_t_sub_fq_int(Q, self.ZT, other.Z))
        elif tx == int:
            return Fq6(Q, fq6_t_sub_fq_int(Q, self.ZT, other))
        else:
            return NotImplemented

    def __rsub__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq6:
            return Fq6(Q, fq6_t_sub_fq6_t(Q, other.ZT, self.ZT))
        elif tx == Fq2:
            return Fq6(Q, fq2_t_sub_fq6_t(Q, other.ZT, self.ZT))
        elif tx == Fq:
            return Fq6(Q, fq_int_sub_fq6_t(Q, other.Z, self.ZT))
        elif tx == int:
            return Fq6(Q, fq_int_sub_fq6_t(Q, other, self.ZT))
        else:
            return NotImplemented

    def __mul__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq6:
            return Fq6(Q, fq6_t_mul_fq6_t(Q, self.ZT, other.ZT))
        elif tx == Fq2:
            return Fq6(Q, fq6_t_mul_fq2_t(Q, self.ZT, other.ZT))
        elif tx == Fq:
            return Fq6(Q, fq6_t_mul_fq_int(Q, self.ZT, other.Z))
        elif tx == int:
            return Fq6(Q, fq6_t_mul_fq_int(Q, self.ZT, other))
        else:
            return NotImplemented

    def __floordiv__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq12:
            return Fq12(Q, fq12_t_mul_fq6_t(Q, fq12_t_invert(Q, other.ZT),
                                            self.ZT))
        elif tx == Fq6:
            return Fq6(Q, fq6_t_mul_fq6_t(Q, self.ZT,
                                          fq6_t_invert(Q, other.ZT)))
        elif tx == Fq2:
            return Fq6(Q, fq6_t_mul_fq2_t(Q, self.ZT,
                                          fq2_t_invert(Q,other.ZT)))
        elif tx == Fq:
            return Fq6(Q, fq6_t_mul_fq_int(Q, self.ZT,
                                           fq_int_invert(Q, other.Z)))
        elif tx == int:
            return Fq6(Q, fq6_t_mul_fq_int(Q, self.ZT,
                                           fq_int_invert(Q, other)))
        else:
            return NotImplemented

    __truediv__ = __floordiv__


class Fq12(FieldExtBase):
    # Fq12 is constructed as Fq6(w) / (w2 - γ) where γ = v
    extension = 12
    embedding = 2
    basefield = Fq6
    root = Fq6(bls12381_q, FQ12_ROOT_TUPLE)

    @staticmethod
    def zero(Q):
        global FQ12_ZERO
        return FQ12_ZERO

    @staticmethod
    def one(Q):
        global FQ12_ONE
        return FQ12_ONE

    def __init__(self, Q, *args):
        super(Fq12, self).__init__(Q)
        if len(args) == 1 and type(args[0]) == tuple:
            self.ZT = args[0][:]
            return
        if len(args) != 2:
            raise Exception('Invalid number of arguments')
        a, b = args
        if type(a) != Fq6 or type(b) != Fq6:
            raise Exception('Invalid arguments type, must be Fq6')
        self.ZT = a.ZT + b.ZT

    def __iter__(self):
        for i in range(0, 12, 6):
            yield Fq6(self.Q, self.ZT[i:i+6])

    def __getitem__(self, idx):
        return Fq6(self.Q, self.ZT[idx*6:idx*6+6])

    def __str__(self):
        return ('Fq12(Q, %s)' % ', '.join(str(fq6) for fq6 in self))

    def __repr__(self):
        return ('Fq12(Q, %s)' % ', '.join(repr(fq6) for fq6 in self))

    def __neg__(self):
       return Fq12(self.Q, fq12_t_neg(self.Q, self.ZT))

    def __invert__(self):
        return Fq12(self.Q, fq12_t_invert(self.Q, self.ZT))

    def __pow__(self, e):
        return Fq12(self.Q, fq12_t_pow(self.Q, self.ZT, e))

    def qi_power(self, i):
        global bls12381_q, frob_coeffs
        if self.Q != bls12381_q:
            raise NotImplementedError
        i %= 12
        if i == 0:
            return self
        return Fq12(self.Q, fq12_t_qi_pow(self.Q, self.ZT, i))

    def __add__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq12:
            return Fq12(Q, fq12_t_add_fq12_t(Q, self.ZT, other.ZT))
        elif tx == Fq:
            return Fq12(Q, fq12_t_add_fq_int(Q, self.ZT, other.Z))
        elif tx == Fq6:
            return Fq12(Q, fq12_t_add_fq6_t(Q, self.ZT, other.ZT))
        elif tx == Fq2:
            return Fq12(Q, fq12_t_add_fq2_t(Q, self.ZT, other.ZT))
        elif tx == int:
            return Fq12(Q, fq12_t_add_fq_int(Q, self.ZT, other))
        else:
            return NotImplemented

    def __sub__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq12:
            return Fq12(Q, fq12_t_sub_fq12_t(Q, self.ZT, other.ZT))
        elif tx == Fq:
            return Fq12(Q, fq12_t_sub_fq_int(Q, self.ZT, other.Z))
        elif tx == Fq6:
            return Fq12(Q, fq12_t_sub_fq6_t(Q, self.ZT, other.ZT))
        elif tx == Fq2:
            return Fq12(Q, fq12_t_sub_fq2_t(Q, self.ZT, other.ZT))
        elif tx == int:
            return Fq12(Q, fq12_t_sub_fq_int(Q, self.ZT, other))
        else:
            return NotImplemented

    def __rsub__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq12:
            return Fq12(Q, fq12_t_sub_fq12_t(Q, other.ZT, self.ZT))
        elif tx == Fq:
            return Fq12(Q, fq_int_sub_fq12_t(Q, other.Z, self.ZT))
        elif tx == Fq6:
            return Fq12(Q, fq6_t_sub_fq12_t(Q, other.ZT, self.ZT))
        elif tx == Fq2:
            return Fq12(Q, fq2_t_sub_fq12_t(Q, other.ZT, self.ZT))
        elif tx == int:
            return Fq12(Q, fq_int_sub_fq12_t(Q, other, self.ZT))
        else:
            return NotImplemented

    def __mul__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq12:
            return Fq12(Q, fq12_t_mul_fq12_t(Q, self.ZT, other.ZT))
        elif tx == int:
            return Fq12(Q, fq12_t_mul_fq_int(Q, self.ZT, other))
        elif tx == Fq:
            return Fq12(Q, fq12_t_mul_fq_int(Q, self.ZT, other.Z))
        elif tx == Fq2:
            return Fq12(Q, fq12_t_mul_fq2_t(Q, self.ZT, other.ZT))
        elif tx == Fq6:
            return Fq12(Q, fq12_t_mul_fq6_t(Q, self.ZT, other.ZT))
        else:
            return NotImplemented

    def __floordiv__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq12:
            return Fq12(Q, fq12_t_mul_fq12_t(Q, self.ZT,
                                             fq12_t_invert(Q, other.ZT)))
        elif tx == int:
            return Fq12(Q, fq12_t_mul_fq12_t(Q, self.ZT,
                                             fq_int_invert(Q, other)))
        elif tx == Fq:
            return Fq12(Q, fq12_t_mul_fq12_t(Q, self.ZT,
                                             fq_int_invert(Q, other.Z)))
        elif tx == Fq2:
            return Fq12(Q, fq12_t_mul_fq2_t(Q, self.ZT,
                                            fq2_t_invert(Q, other.ZT)))
        elif tx == Fq6:
            return Fq12(Q, fq12_t_mul_fq6_t(Q, self.ZT,
                                            fq6_t_invert(Q, other.ZT)))
        else:
            return NotImplemented

    __truediv__ = __floordiv__


T_FQ_INT = (Fq, int)
T_FE_BASED = (Fq2, Fq6, Fq12)


FQ_ZERO = Fq(bls12381_q, 0)
FQ_ONE = Fq(bls12381_q, 1)
FQ2_ZERO = Fq2(bls12381_q, FQ2_ZERO_TUPLE)
FQ2_ONE = Fq2(bls12381_q, FQ2_ONE_TUPLE)
FQ6_ZERO = Fq6(bls12381_q, FQ6_ZERO_TUPLE)
FQ6_ONE = Fq6(bls12381_q, FQ6_ONE_TUPLE)
FQ12_ZERO = Fq12(bls12381_q, FQ12_ZERO_TUPLE)
FQ12_ONE = Fq12(bls12381_q, FQ12_ONE_TUPLE)


# Frobenius coefficients for raising elements to q**i -th powers
# These are specific to this given q
q = bls12381_q
frob_coeffs = {
    (2, 1, 1) : Fq(q, -1),
    (6, 1, 1) : Fq2(q,
                    Fq(q, 0x0),
                    Fq(q, int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                              '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                              'fd8bfd00000000aaac', 16))),
    (6, 1, 2) : Fq2(q,
                    Fq(q, int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                              '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                              'fd8bfd00000000aaad', 16)), Fq(q, 0x0)),
    (6, 2, 1) : Fq2(q,
                    Fq(q, int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                              '3be6f89688de17d813620a00022e01fffffffeff'
                              'fe', 16)),
                    Fq(q, 0x0)),
    (6, 2, 2) : Fq2(q,
                    Fq(q, int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                              '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                              'fd8bfd00000000aaac', 16)),
                    Fq(q, 0x0)),
    (6, 3, 1) : Fq2(q,
                    Fq(q, 0x0),
                    Fq(q, 0x1)),
    (6, 3, 2) : Fq2(q,
                    Fq(q, int('0x1a0111ea397fe69a4b1ba7b6434bacd764774b'
                              '84f38512bf6730d2a0f6b0f6241eabfffeb153ff'
                              'ffb9feffffffffaaaa', 16)),
                    Fq(q, 0x0)),
    (6, 4, 1) : Fq2(q,
                    Fq(q, int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                              '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                              'fd8bfd00000000aaac', 16)),
                    Fq(q, 0x0)),
    (6, 4, 2) : Fq2(q,
                    Fq(q, int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                              '3be6f89688de17d813620a00022e01fffffffeff'
                              'fe', 16)),
                    Fq(q, 0x0)),
    (6, 5, 1) : Fq2(q,
                    Fq(q, 0x0),
                    Fq(q, int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                              '3be6f89688de17d813620a00022e01fffffffeff'
                              'fe', 16))),
    (6, 5, 2) : Fq2(q,
                    Fq(q, int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                              '3be6f89688de17d813620a00022e01fffffffeff'
                              'ff', 16)),
                    Fq(q, 0x0)),
    (12, 1, 1) : Fq6(q,
                     Fq2(q,
                         Fq(q, int('0x1904d3bf02bb0667c231beb4202c0d1f0fd603'
                                   'fd3cbd5f4f7b2443d784bab9c4f67ea53d63e781'
                                   '3d8d0775ed92235fb8', 16)),
                         Fq(q, int('0xfc3e2b36c4e03288e9e902231f9fb854a14787'
                                   'b6c7b36fec0c8ec971f63c5f282d5ac14d6c7ec2'
                                   '2cf78a126ddc4af3', 16))),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0)),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0))),
    (12, 2, 1) : Fq6(q,
                     Fq2(q,
                         Fq(q, int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                                   '3be6f89688de17d813620a00022e01fffffffeff'
                                   'ff', 16)),
                         Fq(q, 0x0)),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0)),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0))),
    (12, 3, 1) : Fq6(q,
                     Fq2(q,
                         Fq(q, int('0x135203e60180a68ee2e9c448d77a2cd91c3ded'
                                   'd930b1cf60ef396489f61eb45e304466cf3e67fa'
                                   '0af1ee7b04121bdea2', 16)),
                         Fq(q, int('0x6af0e0437ff400b6831e36d6bd17ffe48395da'
                                   'bc2d3435e77f76e17009241c5ee67992f72ec05f'
                                   '4c81084fbede3cc09', 16))),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0)),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0))),
    (12, 4, 1) : Fq6(q,
                     Fq2(q,
                         Fq(q, int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                                   '3be6f89688de17d813620a00022e01fffffffeff'
                                   'fe', 16)),
                         Fq(q, 0x0)),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0)),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0))),
    (12, 5, 1) : Fq6(q,
                     Fq2(q,
                         Fq(q, int('0x144e4211384586c16bd3ad4afa99cc9170df35'
                                   '60e77982d0db45f3536814f0bd5871c1908bd478'
                                   'cd1ee605167ff82995', 16)),
                         Fq(q, int('0x5b2cfd9013a5fd8df47fa6b48b1e045f398162'
                                   '40c0b8fee8beadf4d8e9c0566c63a3e6e257f873'
                                   '29b18fae980078116', 16))),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0)),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0))),
    (12, 6, 1) : Fq6(q,
                     Fq2(q,
                         Fq(q, int('0x1a0111ea397fe69a4b1ba7b6434bacd764774b'
                                   '84f38512bf6730d2a0f6b0f6241eabfffeb153ff'
                                   'ffb9feffffffffaaaa', 16)),
                         Fq(q, 0x0)),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0)),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0))),
    (12, 7, 1) : Fq6(q,
                     Fq2(q,
                         Fq(q, int('0xfc3e2b36c4e03288e9e902231f9fb854a14787'
                                   'b6c7b36fec0c8ec971f63c5f282d5ac14d6c7ec2'
                                   '2cf78a126ddc4af3', 16)),
                         Fq(q, int('0x1904d3bf02bb0667c231beb4202c0d1f0fd603'
                                   'fd3cbd5f4f7b2443d784bab9c4f67ea53d63e781'
                                   '3d8d0775ed92235fb8', 16))),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0)),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0))),
    (12, 8, 1) : Fq6(q,
                     Fq2(q,
                         Fq(q, int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                                   '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                                   'fd8bfd00000000aaac', 16)),
                         Fq(q, 0x0)),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0)),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0))),
    (12, 9, 1) : Fq6(q,
                     Fq2(q,
                         Fq(q, int('0x6af0e0437ff400b6831e36d6bd17ffe48395da'
                                   'bc2d3435e77f76e17009241c5ee67992f72ec05f'
                                   '4c81084fbede3cc09', 16)),
                         Fq(q, int('0x135203e60180a68ee2e9c448d77a2cd91c3ded'
                                   'd930b1cf60ef396489f61eb45e304466cf3e67fa'
                                   '0af1ee7b04121bdea2', 16))),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0)),
                     Fq2(q,
                         Fq(q, 0x0),
                         Fq(q, 0x0))),
    (12, 10, 1) : Fq6(q,
                      Fq2(q,
                          Fq(q, int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                                    '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                                    'fd8bfd00000000aaad', 16)),
                          Fq(q, 0x0)),
                      Fq2(q,
                          Fq(q, 0x0),
                          Fq(q, 0x0)),
                      Fq2(q,
                          Fq(q, 0x0),
                      Fq(q, 0x0))),
    (12, 11, 1) : Fq6(q,
                      Fq2(q,
                          Fq(q, int('0x5b2cfd9013a5fd8df47fa6b48b1e045f398162'
                                    '40c0b8fee8beadf4d8e9c0566c63a3e6e257f873'
                                    '29b18fae980078116', 16)),
                          Fq(q, int('0x144e4211384586c16bd3ad4afa99cc9170df35'
                                    '60e77982d0db45f3536814f0bd5871c1908bd478'
                                    'cd1ee605167ff82995', 16))),
                      Fq2(q,
                          Fq(q, 0x0),
                          Fq(q, 0x0)),
                      Fq2(q,
                          Fq(q, 0x0),
                          Fq(q, 0x0))),
}

'''
Copyright 2018 Chia Network Inc

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''
