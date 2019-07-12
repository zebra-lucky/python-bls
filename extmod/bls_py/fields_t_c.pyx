# cython: language_level=3
# cython: profile=False
# -*- coding: utf-8 -*-

import atexit
from queue import Queue

from cpython.object cimport Py_SIZE
from cpython.long cimport PyLong_FromLong
from cpython.longintrepr cimport _PyLong_New, py_long, digit, PyLong_SHIFT
from libc.stdlib cimport malloc, free

from .mpz cimport *
from .types cimport *


# Various functions to deal with conversion mpz <-> Python int/long
# AUTHORS:
# - Gonzalo Tornaria (2006): initial version
# - David Harvey (2007-08-18): added ``mpz_get_pyintlong`` function
# - Jeroen Demeyer (2015-02-24): moved from c_lib, rewritten using
#  ``mpz_export`` and ``mpz_import``
cdef extern from *:
    Py_ssize_t* Py_SIZE_PTR "&Py_SIZE"(object)
    int limb_bits "(8 * sizeof(mp_limb_t))"


# Unused bits in every PyLong digit
cdef size_t PyLong_nails = 8*sizeof(digit) - PyLong_SHIFT


cdef mpz_get_pylong_large(mpz_srcptr z):
    """
    Convert a non-zero ``mpz`` to a Python ``long``.
    """
    cdef size_t nbits = mpz_sizeinbase(z, 2)
    cdef size_t pylong_size = (nbits + PyLong_SHIFT - 1) // PyLong_SHIFT
    L = _PyLong_New(pylong_size)
    mpz_export(L.ob_digit, NULL,
            -1, sizeof(digit), 0, PyLong_nails, z)
    if mpz_sgn(z) < 0:
        # Set correct size (use a pointer to hack around Cython's
        # non-support for lvalues).
        sizeptr = Py_SIZE_PTR(L)
        sizeptr[0] = -pylong_size
    return L


cdef mpz_get_pylong(mpz_srcptr z):
    """
    Convert an ``mpz`` to a Python ``long``.
    """
    if mpz_fits_slong_p(z):
        return PyLong_FromLong(mpz_get_si(z))
    return mpz_get_pylong_large(z)


cdef int mpz_set_pylong(mpz_ptr z, L) except -1:
    """
    Convert a Python ``long`` `L` to an ``mpz``.
    """
    cdef Py_ssize_t pylong_size = Py_SIZE(L)
    if pylong_size < 0:
        pylong_size = -pylong_size
    mpz_import(z, pylong_size, -1, sizeof(digit), 0, PyLong_nails,
            (<py_long>L).ob_digit)
    if Py_SIZE(L) < 0:
        mpz_neg(z, z)


# utils
cdef str mpz_get_pystr(mpz_t x, base):
    cdef char * res_char
    cdef bytes res_b
    res_char = mpz_get_str(NULL, base, x)
    try:
        res_b = res_char
    finally:
        free(res_char)
    return res_b.decode('utf-8')


# fq_t operations
cdef void fq_t_init(fq_t rop):
    mpz_init2(rop, INIT_BITS)


cdef void fq_t_init_set(fq_t rop, fq_t a):
    mpz_init2(rop, INIT_BITS)
    mpz_set(rop, a)


cdef int fq_t_init_set_pylong(fq_t rop, L) except -1:
    mpz_init2(rop, INIT_BITS)
    mpz_set_pylong(rop, L)


cdef int fq_t_init_set_pylong_mod(fq_t rop, mpz_t Q, L) except -1:
    mpz_init2(rop, INIT_BITS)
    mpz_set_pylong(rop, L)
    mpz_fdiv_r(rop, rop, Q)

cdef int fq_t_set_pylong(fq_t rop, L) except -1:
    mpz_set_pylong(rop, L)


cdef fq_t_get_pylong(fq_t a):
    return mpz_get_pylong(a)


cdef str fq_t_get_pystr(fq_t x, base):
    return 'Fq(%s)' % mpz_get_pystr(x, base)


cdef unsigned int fq_t_eq(fq_t a, fq_t b):
    if mpz_cmp(a, b) == 0:
        return 1
    else:
        return 0


cdef void fq_t_invert(fq_t rop, fq_t a, mpz_t z):
    mpz_invert(rop, a, z)


cdef void fq_t_floordiv(fq_t rop, fq_t a, fq_t x, mpz_t z):
    cdef mpz_ptr res
    cdef int _res
    res = fq_t_get(&_res)
    mpz_invert(res, x, z)
    mpz_mul(res, res, a)
    mpz_fdiv_r(rop, res, z)
    fq_t_release(_res)


cdef void fq_t_pow(fq_t rop, fq_t a, fq_t e, mpz_t z):
    mpz_powm(rop, a, e, z)


cdef void fqt_neg(fq_t rop, fq_t a):
    mpz_neg(rop, a)
    mpz_fdiv_r(rop, rop, Q)


cdef void fqt_add(fq_t rop, fq_t a, fq_t b):
    mpz_add(rop, a, b)
    mpz_fdiv_r(rop, rop, Q)


cdef void fqt_sub(fq_t rop, fq_t a, fq_t b):
    mpz_sub(rop, a, b)
    mpz_fdiv_r(rop, rop, Q)


cdef void fqt_mul(fq_t rop, fq_t a, fq_t b):
    mpz_mul(rop, a, b)
    mpz_fdiv_r(rop, rop, Q)


# fq2_t operations
cdef void fq2_t_init(fq2_t rop):
    for i in range(2):
        fq_t_init(&rop[i])


cdef void fq2_t_set(fq2_t rop, fq2_t a):
    for i in range(2):
        mpz_set(&rop[i], &a[i])


cdef void fq2_t_init_set_fq_t(fq2_t rop, fq_t a):
    fq_t_init_set(&rop[0], a)
    fq_t_init(&rop[1])


cdef void fq2_t_init_set_fq2(fq2_t rop, a):
    for i in range(2):
        fq_t_init_set_pylong(&rop[i], a[i])


cdef void fq2_t_set_fq2(fq2_t rop, a):
    for i in range(2):
        mpz_set_pylong(&rop[i], a[i])


cdef tuple fq2_t_get_fq2(fq2_t a):
    res = []
    for i in range(2):
        res.append(mpz_get_pylong(&a[i]))
    return tuple(res)


cdef str fq2_t_get_pystr(fq2_t x, base):
    fq_str_list = [fq_t_get_pystr(&x[i], base) for i in range(2)]
    return f'Fq2({", ".join(fq_str_list)})'


cdef unsigned int fq2_t_eq(fq2_t a, fq2_t b):
    eq_cnt = 0
    for i in range(2):
        if mpz_cmp(&a[i], &b[i]) == 0:
            eq_cnt += 1
    if eq_cnt == 2:
        return 1
    else:
        return 0

cdef void fq2_t_invert(fq2_t rop, fq2_t x_op):
    cdef mpz_ptr res
    cdef int _res
    res = fq2_t_get(&_res)
    fq2_t_set(res, x_op)
    # factor = inv(a*a + b*b)
    mpz_mul(&res[0], &res[0], &res[0])
    mpz_mul(&res[1], &res[1], &res[1])
    mpz_add(&res[0], &res[0], &res[1])
    mpz_invert(&res[0], &res[0], Q)
    # res0 = a*factor
    mpz_mul(&rop[0], &x_op[0], &res[0])
    mpz_fdiv_r(&rop[0], &rop[0], Q)
    # res1 = -b*factor
    mpz_neg(&res[0], &res[0])
    mpz_mul(&rop[1], &x_op[1], &res[0])
    mpz_fdiv_r(&rop[1], &rop[1], Q)
    fq2_t_release(_res)


cdef void fq2_t_floordiv(fq2_t rop, fq2_t a_op, fq2_t x_op):
    cdef mpz_ptr res
    cdef int _res
    res = fq2_t_get(&_res)
    fq2_t_invert(res, x_op)
    fq2_t_mul(rop, a_op, res)
    fq2_t_release(_res)


cdef void fq2_t_qi_pow(fq2_t rop, fq2_t x_op, unsigned int i):
    i %= 2
    if i == 0:
        fq2_t_set(rop, x_op)
        return
    mpz_mul(&rop[1], &x_op[1], fq2_t_root)  # frob coeff equals to fq2_t_root
    mpz_fdiv_r(&rop[1], &rop[1], Q)


cdef void fq2_t_pow(fq2_t rop, fq2_t a_op, mpz_t e_op):
    cdef mpz_ptr res, tmul
    cdef mp_bitcnt_t bits_left, bit_n
    cdef int _tmul, _res
    tmul = fq2_t_get(&_tmul)
    res = fq2_t_get(&_res)
    fq2_t_set(res, fq2_t_one)
    fq2_t_set(tmul, a_op)

    bit_n = 0
    bits_left = mpz_popcount(e_op)
    while bits_left > 0:
        if mpz_tstbit(e_op, bit_n):
            fq2_t_mul(res, res, tmul)
            bits_left -= 1
        fq2_t_mul(tmul, tmul, tmul)
        bit_n += 1

    fq2_t_set(rop, res)
    fq2_t_release(_tmul)
    fq2_t_release(_res)


cdef void fq2_t_mul_by_nonresidue(fq2_t rop, fq2_t a_op):
    '''(a - b) % Q, (a + b) % Q'''
    cdef mpz_ptr res
    cdef int _res
    res = fq2_t_get(&_res)

    mpz_sub(&res[0], &a_op[0], &a_op[1])
    mpz_fdiv_r(&res[0], &res[0], Q)
    mpz_add(&res[1], &a_op[0], &a_op[1])
    mpz_fdiv_r(&res[1], &res[1], Q)

    fq2_t_set(rop, res)
    fq2_t_release(_res)


cdef void fq2_t_add(fq2_t rop, fq2_t a_op, fq2_t m_op):
    for i in range(2):
        fqt_add(&rop[i], &a_op[i], &m_op[i])


cdef void fq2_t_sub(fq2_t rop, fq2_t a_op, fq2_t m_op):
    for i in range(2):
        fqt_sub(&rop[i], &a_op[i], &m_op[i])


cdef void fq2_t_mul_fq_t(fq2_t rop, fq2_t a_op, fq_t m_op):
    for i in range(2):
        fqt_mul(&rop[i], &a_op[i], m_op)


cdef void fq2_t_mul(fq2_t rop, fq2_t a_op, fq2_t m_op):
    '''(a, b)*(m, n) = (a*m - b*n) % Q, (a*n + b*m) % Q'''
    cdef mpz_ptr res, tmul
    cdef int _tmul, _res
    tmul = fq_t_get(&_tmul)
    res = fq2_t_get(&_res)
    mpz_mul(&res[0], &a_op[0], &m_op[0])
    mpz_mul(tmul, &a_op[1], &m_op[1])
    mpz_sub(&res[0], &res[0], tmul)
    mpz_fdiv_r(&res[0], &res[0], Q)
    mpz_mul(&res[1], &a_op[0], &m_op[1])
    mpz_mul(tmul, &a_op[1], &m_op[0])
    mpz_add(&res[1], &res[1], tmul)
    mpz_fdiv_r(&res[1], &res[1], Q)
    fq2_t_set(rop, res)
    fq_t_release(_tmul)
    fq2_t_release(_res)


# fq6_t operations
cdef void fq6_t_init(fq6_t rop):
    for i in range(6):
        fq_t_init(&rop[i])


cdef void fq6_t_set(fq6_t rop, fq6_t a):
    for i in range(6):
        mpz_set(&rop[i], &a[i])


cdef void fq6_t_init_set_fq_t(fq6_t rop, fq_t a):
    fq_t_init_set(&rop[0], a)
    for i in range(1, 6):
        fq_t_init(&rop[i])


cdef void fq6_t_init_set_fq6(fq6_t rop, a):
    for i in range(6):
        fq_t_init_set_pylong(&rop[i], a[i])


cdef void fq6_t_set_fq6(fq6_t rop, a):
    for i in range(6):
        mpz_set_pylong(&rop[i], a[i])


cdef tuple fq6_t_get_fq6(fq6_t a):
    res = []
    for i in range(6):
        res.append(mpz_get_pylong(&a[i]))
    return tuple(res)


cdef str fq6_t_get_pystr(fq6_t x, base):
    fq2_str_list = [fq2_t_get_pystr(&x[i*2], base) for i in range(3)]
    return f'Fq6({", ".join(fq2_str_list)})'


cdef void fq6_t_neg(fq6_t rop, fq6_t x_op):
    for i in range(6):
        mpz_neg(&rop[i], &x_op[i])
        mpz_fdiv_r(&rop[i], &rop[i], Q)


cdef void fq6_t_invert(fq6_t rop, fq6_t x_op):
    cdef mpz_ptr tmul, res
    cdef int _tmul, _res
    tmul = fq2_t_get(&_tmul)
    res = fq6_t_get(&_res)
    # g0 = a*a - b*nor(c)
    fq2_t_mul(&res[0], &x_op[0], &x_op[0])
    fq2_t_mul_by_nonresidue(tmul, &x_op[4])
    fq2_t_mul(tmul, tmul, &x_op[2])
    fq2_t_sub(&res[0], &res[0], tmul)
    # g1 = nor(c*c) - a*b
    fq2_t_mul(&res[2], &x_op[4], &x_op[4])
    fq2_t_mul_by_nonresidue(&res[2], &res[2])
    fq2_t_mul(tmul, &x_op[0], &x_op[2])
    fq2_t_sub(&res[2], &res[2], tmul)
    # g2 = b*b - a*c
    fq2_t_mul(&res[4], &x_op[2], &x_op[2])
    fq2_t_mul(tmul, &x_op[0], &x_op[4])
    fq2_t_sub(&res[4], &res[4], tmul)
    # factor = inv(g0*a + nor(g1*c + g2*b))
    fq2_t_mul(&x_op[0], &x_op[0], &res[0])
    fq2_t_mul(&x_op[4], &x_op[4], &res[2])
    fq2_t_mul(&x_op[2], &x_op[2], &res[4])
    fq2_t_add(tmul, &x_op[4], &x_op[2])
    fq2_t_mul_by_nonresidue(tmul, tmul)
    fq2_t_add(tmul, tmul, &x_op[0])
    fq2_t_invert(tmul, tmul)
    # res0 = g0 * factor
    # res1 = g1 * factor
    # res2 = g2 * factor
    fq2_t_mul(&rop[0], &res[0], tmul)
    fq2_t_mul(&rop[2], &res[2], tmul)
    fq2_t_mul(&rop[4], &res[4], tmul)
    fq2_t_release(_tmul)
    fq6_t_release(_res)


cdef void fq6_t_floordiv(fq6_t rop, fq6_t a_op, fq6_t x_op):
    cdef mpz_ptr res
    cdef int _res
    res = fq6_t_get(&_res)
    fq6_t_invert(res, x_op)
    fq6_t_mul(rop, res, a_op)
    fq6_t_release(_res)


cdef void fq6_t_qi_pow(fq6_t rop, fq6_t x_op, unsigned int i):
    i %= 6
    if i == 0:
        fq6_t_set(rop, x_op)
        return
    fq2_t_qi_pow(&rop[0], &x_op[0], i)
    fq2_t_qi_pow(&rop[2], &x_op[2], i)
    fq2_t_mul(&rop[2], &rop[2], fc_6[i-1][0])
    fq2_t_qi_pow(&rop[4], &x_op[4], i)
    fq2_t_mul(&rop[4], &rop[4], fc_6[i-1][1])


cdef void fq6_t_mul_by_nonresidue(fq6_t rop, fq6_t a_op):
    '''(e - f) % Q, (e + f) %Q, a, b, c, d'''
    cdef mpz_ptr res
    cdef int _res
    res = fq6_t_get(&_res)

    mpz_sub(&res[0], &a_op[4], &a_op[5])
    mpz_fdiv_r(&res[0], &res[0], Q)
    mpz_add(&res[1], &a_op[4], &a_op[5])
    mpz_fdiv_r(&res[1], &res[1], Q)
    mpz_set(&res[2], &a_op[0])
    mpz_set(&res[3], &a_op[1])
    mpz_set(&res[4], &a_op[2])
    mpz_set(&res[5], &a_op[3])

    fq6_t_set(rop, res)
    fq6_t_release(_res)


cdef void fq6_t_add(fq6_t rop, fq6_t a_op, fq6_t m_op):
    for i in range(6):
        fqt_add(&rop[i], &a_op[i], &m_op[i])


cdef void fq6_t_sub(fq6_t rop, fq6_t a_op, fq6_t m_op):
    for i in range(6):
        fqt_sub(&rop[i], &a_op[i], &m_op[i])


cdef void fq6_t_mul(fq6_t rop, fq6_t a_op, fq6_t m_op):
    '''
    r0 = am - bn + cq - cr - dr - dq + eo - ep - fp - fo
    r1 = an + bm + cq + cr - dr + dq + eo + ep - fp + fo
    r2 = ao - bp + cm - dn + eq - er - fr - fq
    r3 = ap + bo + cn + dm + eq + er - fr + fq
    r4 = aq - br + co - dp + em - fn
    r5 = ar + bq + cp + do + en + fm
    '''
    cdef mpz_ptr a, b, c, d, e, f, m, n, o, p, q, r, r0, r1, r2, r3, r4, r5
    cdef mpz_ptr tmul, res
    cdef int _tmul, _res
    tmul = fq_t_get(&_tmul)
    res = fq6_t_get(&_res)

    a = &a_op[0]; b = &a_op[1]; c = &a_op[2];
    d = &a_op[3]; e = &a_op[4]; f = &a_op[5];
    m = &m_op[0]; n = &m_op[1]; o = &m_op[2];
    p = &m_op[3]; q = &m_op[4]; r = &m_op[5];
    r0 = &res[0]; r1 = &res[1]; r2 = &res[2];
    r3 = &res[3]; r4 = &res[4]; r5 = &res[5];

    mpz_mul(r0, a, m); mpz_mul(r1, a, n); mpz_mul(r2, a, o);
    mpz_mul(r3, a, p); mpz_mul(r4, a, q); mpz_mul(r5, a, r);

    mpz_mul(tmul, b, n); mpz_sub(r0, r0, tmul);
    mpz_mul(tmul, b, m); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, b, p); mpz_sub(r2, r2, tmul);
    mpz_mul(tmul, b, o); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, b, r); mpz_sub(r4, r4, tmul);
    mpz_mul(tmul, b, q); mpz_add(r5, r5, tmul);

    mpz_mul(tmul, c, q); mpz_add(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, c, m); mpz_add(r2, r2, tmul);
    mpz_mul(tmul, c, n); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, c, o); mpz_add(r4, r4, tmul);
    mpz_mul(tmul, c, p); mpz_add(r5, r5, tmul);

    mpz_mul(tmul, c, r); mpz_sub(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, d, n); mpz_sub(r2, r2, tmul);
    mpz_mul(tmul, d, m); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, d, p); mpz_sub(r4, r4, tmul);
    mpz_mul(tmul, d, o); mpz_add(r5, r5, tmul);

    mpz_mul(tmul, d, r); mpz_sub(r0, r0, tmul); mpz_sub(r1, r1, tmul);
    mpz_mul(tmul, e, q); mpz_add(r2, r2, tmul); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, e, m); mpz_add(r4, r4, tmul);
    mpz_mul(tmul, e, n); mpz_add(r5, r5, tmul);

    mpz_mul(tmul, d, q); mpz_sub(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, e, r); mpz_sub(r2, r2, tmul); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, f, n); mpz_sub(r4, r4, tmul);
    mpz_mul(tmul, f, m); mpz_add(r5, r5, tmul);

    mpz_mul(tmul, e, o); mpz_add(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, f, r); mpz_sub(r3, r3, tmul); mpz_sub(r2, r2, tmul);

    mpz_mul(tmul, e, p); mpz_sub(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, f, q); mpz_sub(r2, r2, tmul); mpz_add(r3, r3, tmul);

    mpz_mul(tmul, f, p); mpz_sub(r0, r0, tmul); mpz_sub(r1, r1, tmul);
    mpz_mul(tmul, f, o); mpz_sub(r0, r0, tmul); mpz_add(r1, r1, tmul);

    mpz_fdiv_r(r0, r0, Q)
    mpz_fdiv_r(r1, r1, Q)
    mpz_fdiv_r(r2, r2, Q)
    mpz_fdiv_r(r3, r3, Q)
    mpz_fdiv_r(r4, r4, Q)
    mpz_fdiv_r(r5, r5, Q)

    fq6_t_set(rop, res)
    fq_t_release(_tmul)
    fq6_t_release(_res)


# fq12_t operations
cdef void fq12_t_init(fq12_t rop):
    for i in range(12):
        fq_t_init(&rop[i])


cdef void fq12_t_set(fq12_t rop, fq12_t a):
    for i in range(12):
        mpz_set(&rop[i], &a[i])


cdef void fq12_t_init_set_fq_t(fq12_t rop, fq_t a):
    fq_t_init_set(&rop[0], a)
    for i in range(1, 12):
        fq_t_init(&rop[i])


cdef void fq12_t_set_fq12(fq12_t rop, a):
    for i in range(12):
        mpz_set_pylong(&rop[i], a[i])


cdef tuple fq12_t_get_fq12(fq12_t a):
    res = []
    for i in range(12):
        res.append(mpz_get_pylong(&a[i]))
    return tuple(res)


cdef str fq12_t_get_pystr(fq12_t x, base):
    fq6_str_list = [fq6_t_get_pystr(&x[i*6], base) for i in range(2)]
    return f'Fq12({", ".join(fq6_str_list)})'


cdef unsigned int fq12_t_eq(fq12_t a, fq12_t b):
    eq_cnt = 0
    for i in range(12):
        if mpz_cmp(&a[i], &b[i]) == 0:
            eq_cnt += 1
    if eq_cnt == 12:
        return 1
    else:
        return 0


cdef void fq12_t_neg(fq12_t rop, fq12_t a_op):
    for i in range(12):
        fqt_neg(&rop[i], &a_op[i])


cdef void fq12_t_invert(fq12_t rop, fq12_t x_op):
    cdef mpz_ptr res
    cdef int _res
    res = fq12_t_get(&_res)
    # factor = inv(a*a - nor(b*b))
    fq12_t_set(res, x_op)
    fq6_t_mul(&res[0], &res[0], &res[0])
    fq6_t_mul(&res[6], &res[6], &res[6])
    fq6_t_mul_by_nonresidue(&res[6], &res[6])
    fq6_t_sub(&res[0], &res[0], &res[6])
    fq6_t_invert(&res[0], &res[0])
    # res1 = a*factor
    fq6_t_mul(&rop[0], &x_op[0], &res[0])
    # res2 = -b*factor
    fq6_t_neg(&res[0], &res[0])
    fq6_t_mul(&rop[6], &x_op[6], &res[0])
    fq12_t_release(_res)


cdef void fq12_t_floordiv(fq12_t rop, fq12_t a_op, fq12_t x_op):
    cdef mpz_ptr res
    cdef int _res
    res = fq12_t_get(&_res)
    fq12_t_invert(res, x_op)
    fq12_t_mul(rop, res, a_op)
    fq12_t_release(_res)


cdef void fq12_t_qi_pow(fq12_t rop, fq12_t x_op, unsigned int i):
    i %= 12
    if i == 0:
        fq12_t_set(rop, x_op)
        return
    fq6_t_qi_pow(&rop[0], &x_op[0], i)
    fq6_t_qi_pow(&rop[6], &x_op[6], i)
    fq6_t_mul(&rop[6], &rop[6], fc_12[i-1])


cdef void fq12_t_pow(fq12_t rop, fq12_t a_op, mpz_t e_op):
    cdef mpz_ptr res, tmul
    cdef mp_bitcnt_t bits_left, bit_n
    cdef int _res, _tmul
    res = fq12_t_get(&_res)
    tmul = fq12_t_get(&_tmul)
    fq12_t_set(res, fq12_t_one)
    fq12_t_set(tmul, a_op)

    bit_n = 0
    bits_left = mpz_popcount(e_op)
    while bits_left > 0:
        if mpz_tstbit(e_op, bit_n):
            fq12_t_mul(res, res, tmul)
            bits_left -= 1
        fq12_t_mul(tmul, tmul, tmul)
        bit_n += 1

    fq12_t_set(rop, res)
    fq12_t_release(_res)
    fq12_t_release(_tmul)


cdef void fq12_t_add(fq12_t rop, fq12_t a_op, fq12_t m_op):
    for i in range(12):
        fqt_add(&rop[i], &a_op[i], &m_op[i])


cdef void fq_t_sub_fq12_t(fq12_t rop, fq_t a_op, fq12_t m_op):
    fqt_sub(&rop[0], a_op, &m_op[0])
    for i in range(1, 12):
        fqt_neg(&rop[i], &m_op[i])


cdef void fq12_t_sub(fq12_t rop, fq12_t a_op, fq12_t m_op):
    for i in range(12):
        fqt_sub(&rop[i], &a_op[i], &m_op[i])


cdef void fq12_t_mul_fq_t(fq12_t rop, fq12_t a_op, fq_t m_op):
    for i in range(12):
        fqt_mul(&rop[i], &a_op[i], m_op)


cdef void fq12_t_mul(fq12_t rop, fq12_t a_op, fq12_t m_op):
    '''
    r0 = (am - bn + cq - cr - dr - dq + eo - ep - fp - fo +
          gw - gx - hx - hw + iu - iv - jv - ju + ks - kt - lt - ls)
    r1 = (an + bm + cq + cr - dr + dq + eo + ep - fp + fo +
          gw + gx - hx + hw + iu + iv - jv + ju + ks + kt - lt + ls)
    r2 = (ao - bp + cm - dn + eq - er - fr - fq + gs - ht +
          iw - ix - jx - jw + ku - kv - lv - lu)
    r3 = (ap + bo + cn + dm + eq + er - fr + fq + gt + hs +
          iw + ix - jx + jw + ku + kv - lv + lu)
    r4 = (aq - br + co - dp + em - fn + gu - hv + is - jt +
          kw - kx - lx - lw)
    r5 = (ar + bq + cp + do + en + fm + gv + hu + it + js +
          kw + kx - lx + lw)
    r6 = (as - bt + cw - cx - dx - dw + eu - ev - fv - fu +
          gm - hn + iq - ir - jr - jq + ko - kp - lp - lo)
    r7 = (at + bs + cw + cx - dx + dw + eu + ev - fv + fu +
          gn + hm + iq + ir - jr + jq + ko + kp - lp + lo)
    r8 = (au - bv + cs - dt + ew - ex - fx - fw + go - hp +
          im - jn + kq - kr - lr - lq)
    r9 = (av + bu + ct + ds + ew + ex - fx + fw + gp + ho +
          in + jm + kq + kr - lr + lq)
    r10 = aw - bx + cu - dv + es - ft + gq - hr + io - jp + km - ln
    r11 = ax + bw + cv + du + et + fs + gr + hq + ip + jo + kn + lm
    '''
    cdef mpz_ptr a, b, c, d, e, f, g, h, i, j, k, l, r0, r1, r2, r3, r4, r5
    cdef mpz_ptr m, n, o, p, q, r, s, t, u, v, w, x, r6, r7, r8, r9, r10, r11
    cdef mpz_ptr tmul, res
    cdef int _tmul, _res
    res = fq12_t_get(&_res)
    tmul = fq_t_get(&_tmul)

    a = &a_op[0]; b = &a_op[1]; c = &a_op[2];
    d = &a_op[3]; e = &a_op[4]; f = &a_op[5];
    g = &a_op[6]; h = &a_op[7]; i = &a_op[8];
    j = &a_op[9]; k = &a_op[10]; l = &a_op[11];
    m = &m_op[0]; n = &m_op[1]; o = &m_op[2];
    p = &m_op[3]; q = &m_op[4]; r = &m_op[5];
    s = &m_op[6]; t = &m_op[7]; u = &m_op[8];
    v = &m_op[9]; w = &m_op[10]; x = &m_op[11];
    r0 = &res[0]; r1 = &res[1]; r2 = &res[2];
    r3 = &res[3]; r4 = &res[4]; r5 = &res[5];
    r6 = &res[6]; r7 = &res[7]; r8 = &res[8];
    r9 = &res[9]; r10 = &res[10]; r11 = &res[11];

    mpz_mul(r0, a, m); mpz_mul(r1, a, n); mpz_mul(r2, a, o);
    mpz_mul(r3, a, p); mpz_mul(r4, a, q); mpz_mul(r5, a, r);
    mpz_mul(r6, a, s); mpz_mul(r7, a, t); mpz_mul(r8, a, u);
    mpz_mul(r9, a, v); mpz_mul(r10, a, w); mpz_mul(r11, a, x);

    mpz_mul(tmul, b, n); mpz_sub(r0, r0, tmul);
    mpz_mul(tmul, b, m); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, b, p); mpz_sub(r2, r2, tmul);
    mpz_mul(tmul, b, o); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, b, r); mpz_sub(r4, r4, tmul);
    mpz_mul(tmul, b, q); mpz_add(r5, r5, tmul);
    mpz_mul(tmul, b, t); mpz_sub(r6, r6, tmul);
    mpz_mul(tmul, b, s); mpz_add(r7, r7, tmul);
    mpz_mul(tmul, b, v); mpz_sub(r8, r8, tmul);
    mpz_mul(tmul, b, u); mpz_add(r9, r9, tmul);
    mpz_mul(tmul, b, x); mpz_sub(r10, r10, tmul);
    mpz_mul(tmul, b, w); mpz_add(r11, r11, tmul);

    mpz_mul(tmul, c, q); mpz_add(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, c, m); mpz_add(r2, r2, tmul);
    mpz_mul(tmul, c, n); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, c, o); mpz_add(r4, r4, tmul);
    mpz_mul(tmul, c, p); mpz_add(r5, r5, tmul);
    mpz_mul(tmul, c, w); mpz_add(r6, r6, tmul); mpz_add(r7, r7, tmul);
    mpz_mul(tmul, c, s); mpz_add(r8, r8, tmul);
    mpz_mul(tmul, c, t); mpz_add(r9, r9, tmul);
    mpz_mul(tmul, c, u); mpz_add(r10, r10, tmul);
    mpz_mul(tmul, c, v); mpz_add(r11, r11, tmul);

    mpz_mul(tmul, c, r); mpz_sub(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, d, n); mpz_sub(r2, r2, tmul);
    mpz_mul(tmul, d, m); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, d, p); mpz_sub(r4, r4, tmul);
    mpz_mul(tmul, d, o); mpz_add(r5, r5, tmul);
    mpz_mul(tmul, c, x); mpz_sub(r6, r6, tmul); mpz_add(r7, r7, tmul);
    mpz_mul(tmul, d, t); mpz_sub(r8, r8, tmul);
    mpz_mul(tmul, d, s); mpz_add(r9, r9, tmul);
    mpz_mul(tmul, d, v); mpz_sub(r10, r10, tmul);
    mpz_mul(tmul, d, u); mpz_add(r11, r11, tmul);

    mpz_mul(tmul, d, r); mpz_sub(r0, r0, tmul); mpz_sub(r1, r1, tmul);
    mpz_mul(tmul, e, q); mpz_add(r2, r2, tmul); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, e, m); mpz_add(r4, r4, tmul);
    mpz_mul(tmul, e, n); mpz_add(r5, r5, tmul);
    mpz_mul(tmul, d, x); mpz_sub(r6, r6, tmul); mpz_sub(r7, r7, tmul);
    mpz_mul(tmul, e, w); mpz_add(r8, r8, tmul); mpz_add(r9, r9, tmul);
    mpz_mul(tmul, e, s); mpz_add(r10, r10, tmul);
    mpz_mul(tmul, e, t); mpz_add(r11, r11, tmul);

    mpz_mul(tmul, d, q); mpz_sub(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, e, r); mpz_sub(r2, r2, tmul); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, f, n); mpz_sub(r4, r4, tmul);
    mpz_mul(tmul, f, m); mpz_add(r5, r5, tmul);
    mpz_mul(tmul, d, w); mpz_sub(r6, r6, tmul); mpz_add(r7, r7, tmul);
    mpz_mul(tmul, e, x); mpz_sub(r8, r8, tmul); mpz_add(r9, r9, tmul);
    mpz_mul(tmul, f, t); mpz_sub(r10, r10, tmul);
    mpz_mul(tmul, f, s); mpz_add(r11, r11, tmul);

    mpz_mul(tmul, e, o); mpz_add(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, f, r); mpz_sub(r2, r2, tmul); mpz_sub(r3, r3, tmul);
    mpz_mul(tmul, g, u); mpz_add(r4, r4, tmul);
    mpz_mul(tmul, g, v); mpz_add(r5, r5, tmul);
    mpz_mul(tmul, e, u); mpz_add(r6, r6, tmul); mpz_add(r7, r7, tmul);
    mpz_mul(tmul, f, x); mpz_sub(r8, r8, tmul); mpz_sub(r9, r9, tmul);
    mpz_mul(tmul, g, q); mpz_add(r10, r10, tmul);
    mpz_mul(tmul, g, r); mpz_add(r11, r11, tmul);

    mpz_mul(tmul, e, p); mpz_sub(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, f, q); mpz_sub(r2, r2, tmul); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, h, v); mpz_sub(r4, r4, tmul);
    mpz_mul(tmul, h, u); mpz_add(r5, r5, tmul);
    mpz_mul(tmul, e, v); mpz_sub(r6, r6, tmul); mpz_add(r7, r7, tmul);
    mpz_mul(tmul, f, w); mpz_sub(r8, r8, tmul); mpz_add(r9, r9, tmul);
    mpz_mul(tmul, h, r); mpz_sub(r10, r10, tmul);
    mpz_mul(tmul, h, q); mpz_add(r11, r11, tmul);

    mpz_mul(tmul, f, p); mpz_sub(r0, r0, tmul); mpz_sub(r1, r1, tmul);
    mpz_mul(tmul, g, s); mpz_add(r2, r2, tmul);
    mpz_mul(tmul, g, t); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, i, s); mpz_add(r4, r4, tmul);
    mpz_mul(tmul, i, t); mpz_add(r5, r5, tmul);
    mpz_mul(tmul, f, v); mpz_sub(r6, r6, tmul); mpz_sub(r7, r7, tmul);
    mpz_mul(tmul, g, o); mpz_add(r8, r8, tmul);
    mpz_mul(tmul, g, p); mpz_add(r9, r9, tmul);
    mpz_mul(tmul, i, o); mpz_add(r10, r10, tmul);
    mpz_mul(tmul, i, p); mpz_add(r11, r11, tmul);

    mpz_mul(tmul, f, o); mpz_sub(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, h, t); mpz_sub(r2, r2, tmul);
    mpz_mul(tmul, h, s); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, j, t); mpz_sub(r4, r4, tmul);
    mpz_mul(tmul, j, s); mpz_add(r5, r5, tmul);
    mpz_mul(tmul, f, u); mpz_sub(r6, r6, tmul); mpz_add(r7, r7, tmul);
    mpz_mul(tmul, h, p); mpz_sub(r8, r8, tmul);
    mpz_mul(tmul, h, o); mpz_add(r9, r9, tmul);
    mpz_mul(tmul, j, p); mpz_sub(r10, r10, tmul);
    mpz_mul(tmul, j, o); mpz_add(r11, r11, tmul);

    mpz_mul(tmul, g, w); mpz_add(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, i, w); mpz_add(r2, r2, tmul); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, k, w); mpz_add(r4, r4, tmul); mpz_add(r5, r5, tmul);
    mpz_mul(tmul, g, m); mpz_add(r6, r6, tmul);
    mpz_mul(tmul, g, n); mpz_add(r7, r7, tmul);
    mpz_mul(tmul, i, m); mpz_add(r8, r8, tmul);
    mpz_mul(tmul, i, n); mpz_add(r9, r9, tmul);
    mpz_mul(tmul, k, m); mpz_add(r10, r10, tmul);
    mpz_mul(tmul, k, n); mpz_add(r11, r11, tmul);

    mpz_mul(tmul, g, x); mpz_sub(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, i, x); mpz_sub(r2, r2, tmul); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, k, x); mpz_sub(r4, r4, tmul); mpz_add(r5, r5, tmul);
    mpz_mul(tmul, h, n); mpz_sub(r6, r6, tmul);
    mpz_mul(tmul, h, m); mpz_add(r7, r7, tmul);
    mpz_mul(tmul, j, n); mpz_sub(r8, r8, tmul);
    mpz_mul(tmul, j, m); mpz_add(r9, r9, tmul);
    mpz_mul(tmul, l, n); mpz_sub(r10, r10, tmul);
    mpz_mul(tmul, l, m); mpz_add(r11, r11, tmul);

    mpz_mul(tmul, h, x); mpz_sub(r0, r0, tmul); mpz_sub(r1, r1, tmul);
    mpz_mul(tmul, j, x); mpz_sub(r2, r2, tmul); mpz_sub(r3, r3, tmul);
    mpz_mul(tmul, l, x); mpz_sub(r4, r4, tmul); mpz_sub(r5, r5, tmul);
    mpz_mul(tmul, i, q); mpz_add(r6, r6, tmul); mpz_add(r7, r7, tmul);
    mpz_mul(tmul, k, q); mpz_add(r8, r8, tmul); mpz_add(r9, r9, tmul);

    mpz_mul(tmul, h, w); mpz_sub(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, j, w); mpz_sub(r2, r2, tmul); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, l, w); mpz_sub(r4, r4, tmul); mpz_add(r5, r5, tmul);
    mpz_mul(tmul, i, r); mpz_sub(r6, r6, tmul); mpz_add(r7, r7, tmul);
    mpz_mul(tmul, k, r); mpz_sub(r8, r8, tmul); mpz_add(r9, r9, tmul);

    mpz_mul(tmul, i, u); mpz_add(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, k, u); mpz_add(r2, r2, tmul); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, j, r); mpz_sub(r6, r6, tmul); mpz_sub(r7, r7, tmul);
    mpz_mul(tmul, l, r); mpz_sub(r8, r8, tmul); mpz_sub(r9, r9, tmul);

    mpz_mul(tmul, i, v); mpz_sub(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, k, v); mpz_sub(r2, r2, tmul); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, j, q); mpz_sub(r6, r6, tmul); mpz_add(r7, r7, tmul);
    mpz_mul(tmul, l, q); mpz_sub(r8, r8, tmul); mpz_add(r9, r9, tmul);

    mpz_mul(tmul, j, v); mpz_sub(r0, r0, tmul); mpz_sub(r1, r1, tmul);
    mpz_mul(tmul, l, v); mpz_sub(r2, r2, tmul); mpz_sub(r3, r3, tmul);
    mpz_mul(tmul, k, o); mpz_add(r6, r6, tmul); mpz_add(r7, r7, tmul);

    mpz_mul(tmul, j, u); mpz_sub(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, l, u); mpz_sub(r2, r2, tmul); mpz_add(r3, r3, tmul);
    mpz_mul(tmul, k, p); mpz_sub(r6, r6, tmul); mpz_add(r7, r7, tmul);

    mpz_mul(tmul, k, s); mpz_add(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, l, p); mpz_sub(r6, r6, tmul); mpz_sub(r7, r7, tmul);

    mpz_mul(tmul, k, t); mpz_sub(r0, r0, tmul); mpz_add(r1, r1, tmul);
    mpz_mul(tmul, l, o); mpz_sub(r6, r6, tmul); mpz_add(r7, r7, tmul);

    mpz_mul(tmul, l, t); mpz_sub(r0, r0, tmul); mpz_sub(r1, r1, tmul);

    mpz_mul(tmul, l, s); mpz_sub(r0, r0, tmul); mpz_add(r1, r1, tmul);

    mpz_fdiv_r(r0, r0, Q)
    mpz_fdiv_r(r1, r1, Q)
    mpz_fdiv_r(r2, r2, Q)
    mpz_fdiv_r(r3, r3, Q)
    mpz_fdiv_r(r4, r4, Q)
    mpz_fdiv_r(r5, r5, Q)
    mpz_fdiv_r(r6, r6, Q)
    mpz_fdiv_r(r7, r7, Q)
    mpz_fdiv_r(r8, r8, Q)
    mpz_fdiv_r(r9, r9, Q)
    mpz_fdiv_r(r10, r10, Q)
    mpz_fdiv_r(r11, r11, Q)

    fq12_t_set(rop, res)
    fq12_t_release(_res)
    fq_t_release(_tmul)


# ec.py, paring.py calc operations on fq_t, fq2_t, fq12_t operands
cdef void fq2_t_to_affine(fq2_t ropx, fq2_t ropy, int *ropinf,
                          fq2_t x, fq2_t y, fq2_t z, int inf):
    cdef mpz_ptr z_sq, z_cu
    cdef int _z_sq, _z_cu
    if inf:
        fq2_t_set(ropx, x)
        fq2_t_set(ropy, y)
        ropinf[0] = inf

    z_sq = fq2_t_get(&_z_sq)
    z_cu = fq2_t_get(&_z_cu)
    # x = x * invert(z^2)
    # y = y * invert(z^3)
    fq2_t_mul(z_sq, z, z)
    fq2_t_mul(z_cu, z_sq, z)
    fq2_t_invert(z_sq, z_sq)
    fq2_t_invert(z_cu, z_cu)
    fq2_t_mul(ropx, x, z_sq)
    fq2_t_mul(ropy, y, z_cu)
    ropinf[0] = inf

    fq2_t_release(_z_sq)
    fq2_t_release(_z_cu)


cdef void fq2_t_double_point(fq2_t ropx, fq2_t ropy, int *ropinf,
                             fq2_t px, fq2_t py, int *pinf):
    cdef mpz_ptr xr, yr, tmul
    cdef int _xr, _yr, _tmul
    xr = fq2_t_get(&_xr)
    yr = fq2_t_get(&_yr)
    tmul = fq2_t_get(&_tmul)

    # left = 3 * px * px
    fq2_t_pow(yr, px, mpz_n2)
    fq2_t_mul_fq_t(yr, yr, mpz_n3)

    # s = left * invert(2*py)
    fq2_t_mul_fq_t(xr, py, mpz_n2)
    fq2_t_invert(xr, xr)
    fq2_t_mul(yr, yr, xr)

    # xr = fq_pow(Q, s, 2) - 2*px
    fq2_t_pow(xr, yr, mpz_n2)
    fq2_t_sub(xr, xr, px)
    fq2_t_sub(xr, xr, px)

    # yr = s * (px - xr) - py
    fq2_t_sub(tmul, px, xr)
    fq2_t_mul(yr, yr, tmul)
    fq2_t_sub(ropy, yr, py)

    fq2_t_set(ropx, xr)
    ropinf[0] = 0

    fq2_t_release(_xr)
    fq2_t_release(_yr)
    fq2_t_release(_tmul)


cdef void fq2_t_add_points(fq2_t ropx, fq2_t ropy, int *ropinf,
                           fq2_t x1, fq2_t y1, int *inf1,
                           fq2_t x2, fq2_t y2, int *inf2):
    cdef mpz_ptr xr, yr, tmul
    cdef int _xr, _yr, _tmul

    if inf1[0]:
        fq2_t_set(ropx, x2)
        fq2_t_set(ropy, y2)
        ropinf[0] = inf2[0]
        return
    if inf1[0]:
        fq2_t_set(ropx, x1)
        fq2_t_set(ropy, y1)
        ropinf[0] = inf1[0]
        return
    if fq2_t_eq(x1, x2) and fq2_t_eq(y1, y2):
        fq2_t_double_point(ropx, ropy, ropinf, x1, y1, inf1)
        return
    if fq2_t_eq(x1, x2):
        fq2_t_set(ropx, fq2_t_zero)
        fq2_t_set(ropy, fq2_t_zero)
        ropinf[0] = 1
        return

    xr = fq2_t_get(&_xr)
    yr = fq2_t_get(&_yr)
    tmul = fq2_t_get(&_tmul)

    # s = (y2 - y1) * inv(x2 - x1)
    fq2_t_sub(xr, y2, y1)
    fq2_t_sub(yr, x2, x1)
    fq2_t_invert(yr, yr)
    fq2_t_mul(xr, xr, yr)

    # ropx = s*s - x1 - x2
    fq2_t_pow(yr, xr, mpz_n2)
    fq2_t_sub(yr, yr, x1)
    fq2_t_sub(yr, yr, x2)

    # ropy = s * (x1 - xr) - y1
    fq2_t_sub(tmul, x1, yr)
    fq2_t_mul(xr, xr, tmul)
    fq2_t_sub(ropy, xr, y1)

    fq2_t_set(ropx, yr)
    ropinf[0] = 0

    fq2_t_release(_xr)
    fq2_t_release(_yr)
    fq2_t_release(_tmul)


cdef void fq_t_double_point_jacobian(fq_t ropx, fq_t ropy, fq_t ropz,
                                     fq_t x_op, fq_t y_op, fq_t z_op):
    cdef mpz_ptr xres, yres, zres
    cdef int _xres, _yres, _zres
    xres = fq_t_get(&_xres)
    yres = fq_t_get(&_yres)
    zres = fq_t_get(&_zres)

    # M = 3*X^2 + a*Z^4
    # a is 0 on bls12-381
    mpz_pow_ui(yres, x_op, 2)
    mpz_mul_ui(yres, yres, 3)

    # S = 4*X*Y^2
    mpz_pow_ui(zres, y_op, 2)
    mpz_mul(zres, zres, x_op)
    mpz_mul_ui(zres, zres, 4)

    # X' = M^2 - 2*S
    mpz_mul(xres, yres, yres)
    mpz_sub(xres, xres, zres)
    mpz_sub(xres, xres, zres)

    # Y' = M*(S - X') - 8*Y^4
    mpz_sub(zres, zres, xres)
    mpz_mul(yres, yres, zres)
    mpz_pow_ui(zres, y_op, 4)
    mpz_mul_ui(zres, zres, 8)
    mpz_sub(yres, yres, zres)

    # Z' = 2*Y*Z
    mpz_mul(zres, y_op, z_op)
    mpz_mul_ui(zres, zres, 2)

    mpz_fdiv_r(ropx, xres, Q)
    mpz_fdiv_r(ropy, yres, Q)
    mpz_fdiv_r(ropz, zres, Q)

    fq_t_release(_xres)
    fq_t_release(_yres)
    fq_t_release(_yres)


cdef void fq2_t_double_point_jacobian(fq2_t ropx, fq2_t ropy, fq2_t ropz,
                                      fq2_t x_op, fq2_t y_op, fq2_t z_op):
    cdef mpz_ptr xres, yres, zres
    cdef int _xres, _yres, _zres
    xres = fq2_t_get(&_xres)
    yres = fq2_t_get(&_yres)
    zres = fq2_t_get(&_zres)

    # M = 3*X^2 + a*Z^4
    # a is 0 on bls12-381
    fq2_t_pow(yres, x_op, mpz_n2)
    fq2_t_mul_fq_t(yres, yres, mpz_n3)

    # S = 4*X*Y^2
    fq2_t_pow(zres, y_op, mpz_n2)
    fq2_t_mul(zres, zres, x_op)
    fq2_t_mul_fq_t(zres, zres, mpz_n4)

    # X' = M^2 - 2*S
    fq2_t_mul(xres, yres, yres)
    fq2_t_sub(xres, xres, zres)
    fq2_t_sub(xres, xres, zres)

    # Y' = M*(S - X') - 8*Y^4
    fq2_t_sub(zres, zres, xres)
    fq2_t_mul(yres, yres, zres)
    fq2_t_pow(zres, y_op, mpz_n4)
    fq2_t_mul_fq_t(zres, zres, mpz_n8)
    fq2_t_sub(yres, yres, zres)

    # Z' = 2*Y*Z
    fq2_t_mul(zres, y_op, z_op)
    fq2_t_mul_fq_t(zres, zres, mpz_n2)

    for i in range(2):
        mpz_fdiv_r(&ropx[i], &xres[i], Q)
        mpz_fdiv_r(&ropy[i], &yres[i], Q)
        mpz_fdiv_r(&ropz[i], &zres[i], Q)

    fq2_t_release(_xres)
    fq2_t_release(_yres)
    fq2_t_release(_yres)


cdef void fq12_t_double_point_jacobian(fq12_t ropx, fq12_t ropy, fq12_t ropz,
                                       fq12_t x_op, fq12_t y_op, fq12_t z_op):
    cdef mpz_ptr xres, yres, zres
    cdef int _xres, _yres, _zres
    xres = fq12_t_get(&_xres)
    yres = fq12_t_get(&_yres)
    zres = fq12_t_get(&_zres)

    # M = 3*X^2 + a*Z^4
    # a is 0 on bls12-381
    fq12_t_pow(yres, x_op, mpz_n2)
    fq12_t_mul_fq_t(yres, yres, mpz_n3)

    # S = 4*X*Y^2
    fq12_t_pow(zres, y_op, mpz_n2)
    fq12_t_mul(zres, zres, x_op)
    fq12_t_mul_fq_t(zres, zres, mpz_n4)

    # X' = M^2 - 2*S
    fq12_t_mul(xres, yres, yres)
    fq12_t_sub(xres, xres, zres)
    fq12_t_sub(xres, xres, zres)

    # Y' = M*(S - X') - 8*Y^4
    fq12_t_sub(zres, zres, xres)
    fq12_t_mul(yres, yres, zres)
    fq12_t_pow(zres, y_op, mpz_n4)
    fq12_t_mul_fq_t(zres, zres, mpz_n8)
    fq12_t_sub(yres, yres, zres)

    # Z' = 2*Y*Z
    fq12_t_mul(zres, y_op, z_op)
    fq12_t_mul_fq_t(zres, zres, mpz_n2)

    for i in range(12):
        mpz_fdiv_r(&ropx[i], &xres[i], Q)
        mpz_fdiv_r(&ropy[i], &yres[i], Q)
        mpz_fdiv_r(&ropz[i], &zres[i], Q)

    fq12_t_release(_xres)
    fq12_t_release(_yres)
    fq12_t_release(_yres)


cdef fq_t_add_points_jacobian(fq_t ropx, fq_t ropy, fq_t ropz, int *ropinf,
                              fq_t x1, fq_t y1, fq_t z1, int *inf1,
                              fq_t x2, fq_t y2, fq_t z2, int *inf2):
    cdef mpz_ptr u1, u2, s1, s2, h_sq
    cdef int _u1, _u2, _s1, _s2, _h_sq

    if inf1[0]:
        mpz_set(ropx, x2)
        mpz_set(ropy, y2)
        mpz_set(ropz, z2)
        ropinf[0] = inf2[0]
        return
    if inf2[0]:
        mpz_set(ropx, x1)
        mpz_set(ropy, y1)
        mpz_set(ropz, z1)
        ropinf[0] = inf1[0]
        return

    u1 = fq_t_get(&_u1)
    u2 = fq_t_get(&_u2)
    s1 = fq_t_get(&_s1)
    s2 = fq_t_get(&_s2)
    h_sq = fq_t_get(&_h_sq)

    # u1 = x1*z2^2
    mpz_pow_ui(u1, z2, 2)
    mpz_mul(u1, u1, x1)
    mpz_fdiv_r(u1, u1, Q)
    # u2 = x2*z1^2
    mpz_pow_ui(u2, z1, 2)
    mpz_mul(u2, u2, x2)
    mpz_fdiv_r(u2, u2, Q)
    # s1 = y1*z2^3
    mpz_pow_ui(s1, z2, 3)
    mpz_mul(s1, s1, y1)
    mpz_fdiv_r(s1, s1, Q)
    # s2 = y2*z1^3
    mpz_pow_ui(s2, z1, 3)
    mpz_mul(s2, s2, y2)
    mpz_fdiv_r(s2, s2, Q)

    if fq_t_eq(u1, u2):
        if not fq_t_eq(s1, s2):
            mpz_set(ropx, fq_t_one)
            mpz_set(ropy, fq_t_one)
            mpz_set(ropz, fq_t_zero)
            ropinf[0] = 1
        else:
            fq_t_double_point_jacobian(ropx, ropy, ropz, x1, y1, z1)
            ropinf[0] = 0
    else:
        # h = u2 - u1
        mpz_sub(u2, u2, u1)
        # r = s2 - s1
        mpz_sub(s2, s2, s1)
        # zr = h*z1*z2
        mpz_mul(ropz, z1, z2)
        mpz_mul(ropz, ropz, u2)
        mpz_fdiv_r(ropz, ropz, Q)
        # h^2
        mpz_mul(h_sq, u2, u2)
        # h^3
        mpz_mul(ropy, h_sq, u2)
        # xr = r^2 - h^3 - 2*u1*h^2
        mpz_mul(ropx, s2, s2)
        mpz_sub(ropx, ropx, ropy)
        mpz_mul(u1, u1, h_sq)
        mpz_sub(ropx, ropx, u1)
        mpz_sub(ropx, ropx, u1)
        mpz_fdiv_r(ropx, ropx, Q)
        # yr = r*(u1*h^2 - xr) - s1*h^3
        mpz_mul(s1, s1, ropy)
        mpz_sub(ropy, u1, ropx)
        mpz_mul(ropy, ropy, s2)
        mpz_sub(ropy, ropy, s1)
        mpz_fdiv_r(ropy, ropy, Q)
        ropinf[0] = 0

    fq_t_release(_u1)
    fq_t_release(_u2)
    fq_t_release(_s1)
    fq_t_release(_s2)
    fq_t_release(_h_sq)


cdef fq2_t_add_points_jacobian(fq2_t ropx, fq2_t ropy, fq2_t ropz, int *ropinf,
                               fq2_t x1, fq2_t y1, fq2_t z1, int *inf1,
                               fq2_t x2, fq2_t y2, fq2_t z2, int *inf2):
    cdef mpz_ptr u1, u2, s1, s2, h_sq
    cdef int _u1, _u2, _s1, _s2, _h_sq

    if inf1[0]:
        fq2_t_set(ropx, x2)
        fq2_t_set(ropy, y2)
        fq2_t_set(ropz, z2)
        ropinf[0] = inf2[0]
        return
    if inf2[0]:
        fq2_t_set(ropx, x1)
        fq2_t_set(ropy, y1)
        fq2_t_set(ropz, z1)
        ropinf[0] = inf1[0]
        return

    u1 = fq2_t_get(&_u1)
    u2 = fq2_t_get(&_u2)
    s1 = fq2_t_get(&_s1)
    s2 = fq2_t_get(&_s2)
    h_sq = fq2_t_get(&_h_sq)

    # u1 = x1*z2^2
    fq2_t_pow(u1, z2, mpz_n2)
    fq2_t_mul(u1, u1, x1)
    # u2 = x2*z1^2
    fq2_t_pow(u2, z1, mpz_n2)
    fq2_t_mul(u2, u2, x2)
    # s1 = y1*z2^3
    fq2_t_pow(s1, z2, mpz_n3)
    fq2_t_mul(s1, s1, y1)
    # s2 = y2*z1^3
    fq2_t_pow(s2, z1, mpz_n3)
    fq2_t_mul(s2, s2, y2)

    if fq2_t_eq(u1, u2):
        if not fq2_t_eq(s1, s2):
            fq2_t_set(ropx, fq2_t_one)
            fq2_t_set(ropy, fq2_t_one)
            fq2_t_set(ropz, fq2_t_zero)
            ropinf[0] = 1
        else:
            fq2_t_double_point_jacobian(ropx, ropy, ropz, x1, y1, z1)
            ropinf[0] = 0
    else:
        # h = u2 - u1
        fq2_t_sub(u2, u2, u1)
        # r = s2 - s1
        fq2_t_sub(s2, s2, s1)
        # zr = h*z1*z2
        fq2_t_mul(ropz, z1, z2)
        fq2_t_mul(ropz, ropz, u2)
        # h^2
        fq2_t_mul(h_sq, u2, u2)
        # h^3
        fq2_t_mul(ropy, h_sq, u2)
        # xr = r^2 - h^3 - 2*u1*h^2
        fq2_t_mul(ropx, s2, s2)
        fq2_t_sub(ropx, ropx, ropy)
        fq2_t_mul(u1, u1, h_sq)
        fq2_t_sub(ropx, ropx, u1)
        fq2_t_sub(ropx, ropx, u1)
        # yr = r*(u1*h^2 - xr) - s1*h^3
        fq2_t_mul(s1, s1, ropy)
        fq2_t_sub(ropy, u1, ropx)
        fq2_t_mul(ropy, ropy, s2)
        fq2_t_sub(ropy, ropy, s1)
        ropinf[0] = 0

    fq2_t_release(_u1)
    fq2_t_release(_u2)
    fq2_t_release(_s1)
    fq2_t_release(_s2)
    fq2_t_release(_h_sq)


cdef fq12_t_add_points_jacobian(fq12_t ropx, fq12_t ropy, fq12_t ropz,
                                int *ropinf,
                                fq12_t x1, fq12_t y1, fq12_t z1, int *inf1,
                                fq12_t x2, fq12_t y2, fq12_t z2, int *inf2):
    cdef mpz_ptr u1, u2, s1, s2, h_sq
    cdef int _u1, _u2, _s1, _s2, _h_sq

    if inf1[0]:
        fq12_t_set(ropx, x2)
        fq12_t_set(ropy, y2)
        fq12_t_set(ropz, z2)
        ropinf[0] = inf2[0]
        return
    if inf2[0]:
        fq12_t_set(ropx, x1)
        fq12_t_set(ropy, y1)
        fq12_t_set(ropz, z1)
        ropinf[0] = inf1[0]
        return

    u1 = fq12_t_get(&_u1)
    u2 = fq12_t_get(&_u2)
    s1 = fq12_t_get(&_s1)
    s2 = fq12_t_get(&_s2)
    h_sq = fq12_t_get(&_h_sq)

    # u1 = x1*z2^2
    fq12_t_pow(u1, z2, mpz_n2)
    fq12_t_mul(u1, u1, x1)
    # u2 = x2*z1^2
    fq12_t_pow(u2, z1, mpz_n2)
    fq12_t_mul(u2, u2, x2)
    # s1 = y1*z2^3
    fq12_t_pow(s1, z2, mpz_n3)
    fq12_t_mul(s1, s1, y1)
    # s2 = y2*z1^3
    fq12_t_pow(s2, z1, mpz_n3)
    fq12_t_mul(s2, s2, y2)

    if fq12_t_eq(u1, u2):
        if not fq12_t_eq(s1, s2):
            fq12_t_set(ropx, fq12_t_one)
            fq12_t_set(ropy, fq12_t_one)
            fq12_t_set(ropz, fq12_t_zero)
            ropinf[0] = 1
        else:
            fq12_t_double_point_jacobian(ropx, ropy, ropz, x1, y1, z1)
            ropinf[0] = 1
    else:
        # h = u2 - u1
        fq12_t_sub(u2, u2, u1)
        # r = s2 - s1
        fq12_t_sub(s2, s2, s1)
        # zr = h*z1*z2
        fq12_t_mul(ropz, z1, z2)
        fq12_t_mul(ropz, ropz, u2)
        # h^2
        fq12_t_mul(h_sq, u2, u2)
        # h^3
        fq12_t_mul(ropy, h_sq, u2)
        # xr = r^2 - h^3 - 2*u1*h^2
        fq12_t_mul(ropx, s2, s2)
        fq12_t_sub(ropx, ropx, ropy)
        fq12_t_mul(u1, u1, h_sq)
        fq12_t_sub(ropx, ropx, u1)
        fq12_t_sub(ropx, ropx, u1)
        # yr = r*(u1*h^2 - xr) - s1*h^3
        fq12_t_mul(s1, s1, ropy)
        fq12_t_sub(ropy, u1, ropx)
        fq12_t_mul(ropy, ropy, s2)
        fq12_t_sub(ropy, ropy, s1)
        ropinf[0] = 0

    fq12_t_release(_u1)
    fq12_t_release(_u2)
    fq12_t_release(_s1)
    fq12_t_release(_s2)
    fq12_t_release(_h_sq)


cdef void fq2_t_scalar_mult_jacobian(fq2_t ropx, fq2_t ropy, fq2_t ropz,
                                     int *ropinf, mpz_t c,
                                     fq2_t x, fq2_t y, fq2_t z, int *inf):
    cdef mp_bitcnt_t bits_left, bit_n
    cdef mpz_ptr tfq, xr, yr, zr
    cdef int _tfq, _xr, _yr, _zr, infr

    tfq = fq_t_get(&_tfq)

    mpz_fdiv_r(tfq, c, Q)
    if inf[0] or mpz_cmp(tfq, mpz_n0) == 0:
        fq2_t_set(ropx, fq2_t_one)
        fq2_t_set(ropy, fq2_t_one)
        fq2_t_set(ropz, fq2_t_zero)
        ropinf[0] = 1
        fq_t_release(_tfq)
        return

    xr = fq2_t_get(&_xr)
    yr = fq2_t_get(&_yr)
    zr = fq2_t_get(&_zr)
    fq2_t_set(xr, fq2_t_one)
    fq2_t_set(yr, fq2_t_one)
    fq2_t_set(zr, fq2_t_zero)
    infr = 1

    bit_n = 0
    bits_left = mpz_popcount(c)
    while bits_left > 0:
        if mpz_tstbit(c, bit_n):
            # result += addend
            fq2_t_add_points_jacobian(xr, yr, zr, &infr,
                                      xr, yr, zr, &infr, x, y, z, inf)
            bits_left -= 1
        # double point
        fq2_t_double_point_jacobian(x, y, z, x, y, z)
        bit_n += 1
    fq2_t_set(ropx, xr)
    fq2_t_set(ropy, yr)
    fq2_t_set(ropz, zr)
    ropinf[0] = infr

    fq_t_release(_tfq)
    fq2_t_release(_xr)
    fq2_t_release(_yr)
    fq2_t_release(_zr)


cdef void fq2_t_double_line_eval(fq12_t res,
                                 fq2_t rx_t, fq2_t ry_t, fq_t px, fq_t py):
    cdef mpz_ptr r12_x, r12_y, slope
    cdef int _r12_x, _r12_y, _slope
    r12_x = fq12_t_get(&_r12_x)
    r12_y = fq12_t_get(&_r12_y)
    slope = fq12_t_get(&_slope)

    #R12 = untwist(R)
    fq2_t_untwist(r12_x, r12_y, rx_t, ry_t)

    # slope = (3 * pow(R12.x, 2)) / (2 * R12.y)
    fq12_t_pow(res, r12_x, mpz_n2)
    fq12_t_mul_fq_t(res, res, mpz_n3)
    fq12_t_mul_fq_t(slope, r12_y, mpz_n2)
    fq12_t_invert(slope, slope)
    fq12_t_mul(slope, res, slope)

    # v = R12.y - slope * R12.x
    fq12_t_mul(res, slope, r12_x)
    fq12_t_sub(res, r12_y, res)

    # res = P.y - P.x * slope - v
    fq12_t_mul_fq_t(slope, slope, px)
    fq_t_sub_fq12_t(res, py, res)
    fq12_t_sub(res, res, slope)

    fq12_t_release(_r12_x)
    fq12_t_release(_r12_y)
    fq12_t_release(_slope)


cdef void fq2_t_add_line_eval(fq12_t res,
                              fq2_t rx, fq2_t ry, fq2_t qx, fq2_t qy,
                              fq_t px, fq_t py):
    cdef mpz_ptr r12x, r12y, q12x, q12y, nq12x, nq12y, slope, tmul
    cdef int _r12x, _r12y, _q12x, _q12y, _nq12x, _nq12y, _slope, _tmul
    r12x = fq12_t_get(&_r12x)
    r12y = fq12_t_get(&_r12y)
    q12x = fq12_t_get(&_q12x)
    q12y = fq12_t_get(&_q12y)
    nq12x = fq12_t_get(&_nq12x)
    nq12y = fq12_t_get(&_nq12y)
    slope = fq12_t_get(&_slope)
    tmul = fq12_t_get(&_tmul)

    #R12 = untwist(R)
    #Q12 = untwist(Q)
    fq2_t_untwist(r12x, r12y, rx, ry)
    fq2_t_untwist(q12x, q12y, qx, qy)

    # This is the case of a vertical line, where the denominator
    # will be 0.
    #if R12 == Q12.negate():
    #    return P.x - R12.x
    fq12_t_neg(nq12x, q12x)
    fq12_t_neg(nq12y, q12y)
    if fq12_t_eq(r12x, nq12x) and fq12_t_eq(r12y, nq12y):
        fq_t_sub_fq12_t(res, px, r12x)
        fq12_t_release(_r12x)
        fq12_t_release(_r12y)
        fq12_t_release(_q12x)
        fq12_t_release(_q12y)
        fq12_t_release(_nq12x)
        fq12_t_release(_nq12y)
        fq12_t_release(_slope)
        fq12_t_release(_tmul)
        return

    # slope = (Q12.y - R12.y) / (Q12.x - R12.x)
    fq12_t_sub(res, q12x, r12x)
    fq12_t_invert(res, res)
    fq12_t_sub(slope, q12y, r12y)
    fq12_t_mul(slope, slope, res)

    # v = (Q12.y * R12.x - R12.y * Q12.x) / (R12.x - Q12.x)
    fq12_t_mul(res, q12y, r12x)
    fq12_t_mul(tmul, r12y, q12x)
    fq12_t_sub(res, res, tmul)
    fq12_t_sub(tmul, r12x, q12x)
    fq12_t_invert(tmul, tmul)
    fq12_t_mul(tmul, res, tmul)

    #res = P.y - P.x * slope - v
    fq12_t_mul_fq_t(slope, slope, px)
    fq_t_sub_fq12_t(res, py, slope)
    fq12_t_sub(res, res, tmul)

    fq12_t_release(_r12x)
    fq12_t_release(_r12y)
    fq12_t_release(_q12x)
    fq12_t_release(_q12y)
    fq12_t_release(_nq12x)
    fq12_t_release(_nq12y)
    fq12_t_release(_slope)
    fq12_t_release(_tmul)


cdef void fq2_t_untwist(fq12_t ropx, fq12_t ropy, fq2_t x_t, fq2_t y_t):
    cdef mpz_ptr tmul
    cdef int _tmul
    tmul = fq_t_get(&_tmul)
    #m, n = x_t
    #new_x = (0, 0, 0, 0, (tw1*m - tw2*n) % Q, (tw1*n + tw2*m) % Q,
    #         0, 0, 0, 0, 0, 0)
    mpz_set(&ropx[0], mpz_n0)
    mpz_set(&ropx[1], mpz_n0)
    mpz_set(&ropx[2], mpz_n0)
    mpz_set(&ropx[3], mpz_n0)
    mpz_mul(&ropx[4], tw1, &x_t[0])
    mpz_mul(tmul, tw2, &x_t[1])
    mpz_sub(&ropx[4], &ropx[4], tmul)
    mpz_mul(&ropx[5], tw1, &x_t[1])
    mpz_mul(tmul, tw2, &x_t[0])
    mpz_add(&ropx[5], &ropx[5], tmul)
    mpz_set(&ropx[6], mpz_n0)
    mpz_set(&ropx[7], mpz_n0)
    mpz_set(&ropx[8], mpz_n0)
    mpz_set(&ropx[9], mpz_n0)
    mpz_set(&ropx[10], mpz_n0)
    mpz_set(&ropx[11], mpz_n0)

    #m, n = y_t
    #new_y = (0, 0, 0, 0, 0, 0,
    #         0, 0, (tw1*m - tw2*n) % Q, (tw1*n + tw2*m) % Q, 0, 0)
    mpz_set(&ropy[0], mpz_n0)
    mpz_set(&ropy[1], mpz_n0)
    mpz_set(&ropy[2], mpz_n0)
    mpz_set(&ropy[3], mpz_n0)
    mpz_set(&ropy[4], mpz_n0)
    mpz_set(&ropy[5], mpz_n0)
    mpz_set(&ropy[6], mpz_n0)
    mpz_set(&ropy[7], mpz_n0)
    mpz_mul(&ropy[8], tw1, &y_t[0])
    mpz_mul(tmul, tw2, &y_t[1])
    mpz_sub(&ropy[8], &ropy[8], tmul)
    mpz_mul(&ropy[9], tw1, &y_t[1])
    mpz_mul(tmul, tw2, &y_t[0])
    mpz_add(&ropy[9], &ropy[9], tmul)
    mpz_set(&ropy[10], mpz_n0)
    mpz_set(&ropy[11], mpz_n0)
    fq_t_release(_tmul)


cdef void fq_t_miller_loop(fq12_t res,
                           fq_t px, fq_t py, int pinf,
                           fq2_t qx, fq2_t qy, int qinf):
    cdef mp_bitcnt_t bits_left, bit_n
    cdef mpz_ptr rx, ry, tfq12
    cdef int _rx, _ry, _tfq12, rinf

    rx = fq2_t_get(&_rx)
    ry = fq2_t_get(&_ry)
    tfq12 = fq12_t_get(&_tfq12)

    fq2_t_set(rx, qx)
    fq2_t_set(ry, qy)
    rinf = qinf
    fq12_t_set(res, fq12_t_one)

    bit_n = 0
    bits_left = mpz_popcount(NX)
    while bits_left > 0:
        if mpz_tstbit(NX, bit_n):
            bits_left -= 1
        bit_n += 1

    bit_n -= 1
    while bit_n > 0:
        bit_n -= 1
        # Compute sloped line lrr
        fq2_t_double_line_eval(tfq12, rx, ry, px, py)
        fq12_t_pow(res, res, mpz_n2)
        fq12_t_mul(res, res, tfq12)
        # R = 2 * R
        fq2_t_double_point(rx, ry, &rinf, rx, ry, &rinf)
        if mpz_tstbit(NX, bit_n):
            # Compute sloped line lrq
            fq2_t_add_line_eval(tfq12, rx, ry, qx, qy, px, py)
            fq12_t_mul(res, res, tfq12)
            # R = R + Q
            fq2_t_add_points(rx, ry, &rinf, rx, ry, &rinf, qx, qy, &qinf)

    fq2_t_release(_rx)
    fq2_t_release(_ry)
    fq12_t_release(_tfq12)


cdef void fq12_t_final_exp(fq12_t rop, fq12_t a_op):
    cdef mpz_ptr res
    cdef int _res
    res = fq12_t_get(&_res)
    # ans = pow(e, mpz_final_exp_e)
    fq12_t_pow(res, a_op, mpz_final_exp_e)
    # ans = ans.qi_power(2) * ans
    fq12_t_qi_pow(rop, res, 2)
    fq12_t_mul(res, rop, res)
    # ans = ans.qi_power(6) / ans
    fq12_t_qi_pow(rop, res, 6)
    fq12_t_floordiv(rop, rop, res)
    fq12_t_release(_res)


cdef void fq_t_ate_pairing_multi(int n, fq12_t rop,
                                 fq_t *ps_x, fq_t *ps_y, int *ps_inf,
                                 fq2_t *qs_x, fq2_t *qs_y, int *qs_inf):
    cdef mpz_ptr ml_res
    cdef int _ml_res
    ml_res = fq12_t_get(&_ml_res)

    fq12_t_set(rop, fq12_t_one)
    for i in range(n):
        fq_t_miller_loop(ml_res,
                         ps_x[i], ps_y[i], ps_inf[i],
                         qs_x[i], qs_y[i], qs_inf[i])
        fq12_t_mul(rop, rop, ml_res)
    fq12_t_final_exp(rop, rop)

    fq12_t_release(_ml_res)


# fq operations
def fq_invert(P, A):
    cdef mpz_ptr z, a
    cdef int _z, _a
    a = fq_t_get(&_a)
    fq_t_set_pylong(a, A)
    if P == Q_int:
        fq_t_invert(a, a, Q)
    else:
        z = fq_t_get(&_z)
        fq_t_set_pylong(z, P)
        fq_t_invert(a, a, z)
        fq_t_release(_z)
    res = fq_t_get_pylong(a)
    fq_t_release(_a)
    return res


def fq_floordiv(P, A, X):
    cdef mpz_ptr z, a, x
    cdef int _z, _a, _x
    a = fq_t_get(&_a)
    x = fq_t_get(&_x)
    fq_t_set_pylong(a, A)
    fq_t_set_pylong(x, X)
    if P == Q_int:
        fq_t_floordiv(a, a, x, Q)
    else:
        z = fq_t_get(&_z)
        fq_t_set_pylong(z, P)
        fq_t_floordiv(a, a, x, z)
        fq_t_release(_z)
    res = fq_t_get_pylong(a)
    fq_t_release(_a)
    fq_t_release(_x)
    return res


def fq_pow(P, A, E):
    cdef mpz_ptr z, a, e
    cdef int _z, _a, _e
    a = fq_t_get(&_a)
    e = fq_t_get(&_e)
    fq_t_set_pylong(a, A)
    fq_t_set_pylong(e, E)
    if P == Q_int:
        fq_t_pow(a, a, e, Q)
    else:
        z = fq_t_get(&_z)
        fq_t_set_pylong(z, P)
        fq_t_pow(a, a, e, z)
        fq_t_release(_z)
    res = fq_t_get_pylong(a)
    fq_t_release(_a)
    fq_t_release(_e)
    return res


# fq2 operations
def fq2_invert(t_x):
    cdef mpz_ptr x
    cdef int _x
    x = fq2_t_get(&_x)
    fq2_t_set_fq2(x, t_x)
    fq2_t_invert(x, x)
    res = fq2_t_get_fq2(x)
    fq2_t_release(_x)
    return res


def fq2_floordiv(t_a, t_x):
    cdef mpz_ptr a, x
    cdef int _a, _x
    a = fq2_t_get(&_a)
    x = fq2_t_get(&_x)
    fq2_t_set_fq2(a, t_a)
    fq2_t_set_fq2(x, t_x)
    fq2_t_floordiv(a, a, x)
    res = fq2_t_get_fq2(a)
    fq2_t_release(_a)
    fq2_t_release(_x)
    return res


def fq2_qi_pow(t_x, i):
    cdef mpz_ptr x
    cdef int _x
    i %= 2
    if i == 0:
        return t_x
    x = fq2_t_get(&_x)
    fq2_t_set_fq2(x, t_x)
    fq2_t_qi_pow(x, x, i)
    res = fq2_t_get_fq2(x)
    fq2_t_release(_x)
    return res


def fq2_pow(t_a, E):
    cdef mpz_ptr a, e
    cdef int _e, _a
    a = fq2_t_get(&_a)
    e = fq_t_get(&_e)
    fq2_t_set_fq2(a, t_a)
    fq_t_set_pylong(e, E)
    fq2_t_pow(a, a, e)
    res = fq2_t_get_fq2(a)
    fq2_t_release(_a)
    fq_t_release(_e)
    return res


def fq2_add(t_a, t_m):
    cdef mpz_ptr a, b
    cdef int _a, _b
    a = fq2_t_get(&_a)
    b = fq2_t_get(&_b)
    fq2_t_set_fq2(a, t_a)
    fq2_t_set_fq2(b, t_m)
    fq2_t_add(a, a, b)
    res = fq2_t_get_fq2(a)
    fq2_t_release(_a)
    fq2_t_release(_b)
    return res


def fq2_mul(t_a, t_m):
    cdef mpz_ptr a, b
    cdef int _a, _b
    a = fq2_t_get(&_a)
    b = fq2_t_get(&_b)
    fq2_t_set_fq2(a, t_a)
    fq2_t_set_fq2(b, t_m)
    fq2_t_mul(a, a, b)
    res = fq2_t_get_fq2(a)
    fq2_t_release(_a)
    fq2_t_release(_b)
    return res


# fq6 operations
def fq6_invert(t_x):
    cdef mpz_ptr x
    cdef int _x
    x = fq6_t_get(&_x)
    fq6_t_set_fq6(x, t_x)
    fq6_t_invert(x, x)
    res = fq6_t_get_fq6(x)
    fq6_t_release(_x)
    return res


def fq6_floordiv(t_a, t_x):
    cdef mpz_ptr a, x
    cdef int _a, _x
    a = fq6_t_get(&_a)
    x = fq6_t_get(&_x)
    fq6_t_set_fq6(a, t_a)
    fq6_t_set_fq6(x, t_x)
    fq6_t_floordiv(a, a, x)
    res = fq6_t_get_fq6(a)
    fq6_t_release(_a)
    fq6_t_release(_x)
    return res


def fq6_qi_pow(t_x, i):
    cdef mpz_ptr x
    cdef int _x
    i %= 6
    if i == 0:
        return t_x
    x = fq6_t_get(&_x)
    fq6_t_set_fq6(x, t_x)
    fq6_t_qi_pow(x, x, i)
    res = fq6_t_get_fq6(x)
    fq6_t_release(_x)
    return res


def fq6_add(t_a, t_m):
    cdef mpz_ptr a, b
    cdef int _a, _b
    a = fq6_t_get(&_a)
    b = fq6_t_get(&_b)
    fq6_t_set_fq6(a, t_a)
    fq6_t_set_fq6(b, t_m)
    fq6_t_add(a, a, b)
    res = fq6_t_get_fq6(a)
    fq6_t_release(_a)
    fq6_t_release(_b)
    return res


def fq6_mul(t_a, t_m):
    cdef mpz_ptr a, b
    cdef int _a, _b
    a = fq6_t_get(&_a)
    b = fq6_t_get(&_b)
    fq6_t_set_fq6(a, t_a)
    fq6_t_set_fq6(b, t_m)
    fq6_t_mul(a, a, b)
    res = fq6_t_get_fq6(a)
    fq6_t_release(_a)
    fq6_t_release(_b)
    return res


# fq12 operations
def fq12_invert(t_x):
    cdef mpz_ptr x
    cdef int _x
    x = fq12_t_get(&_x)
    fq12_t_set_fq12(x, t_x)
    fq12_t_invert(x, x)
    res = fq12_t_get_fq12(x)
    fq12_t_release(_x)
    return res


def fq12_floordiv(t_a, t_x):
    cdef mpz_ptr a, x
    cdef int _a, _x
    a = fq12_t_get(&_a)
    x = fq12_t_get(&_x)
    fq12_t_set_fq12(a, t_a)
    fq12_t_set_fq12(x, t_x)
    fq12_t_floordiv(a, a, x)
    res = fq12_t_get_fq12(a)
    fq12_t_release(_a)
    fq12_t_release(_x)
    return res


def fq12_qi_pow(t_x, i):
    cdef mpz_ptr x
    cdef int _x
    i %= 12
    if i == 0:
        return t_x
    x = fq12_t_get(&_x)
    fq12_t_set_fq12(x, t_x)
    fq12_t_qi_pow(x, x, i)
    res = fq12_t_get_fq12(x)
    fq12_t_release(_x)
    return res


def fq12_pow(t_a, E):
    cdef mpz_ptr a, e
    cdef int _a, _e
    a = fq12_t_get(&_a)
    e = fq_t_get(&_e)
    fq12_t_set_fq12(a, t_a)
    fq_t_set_pylong(e, E)
    fq12_t_pow(a, a, e)
    res = fq12_t_get_fq12(a)
    fq12_t_release(_a)
    fq_t_release(_e)
    return res


def fq12_mul_fq(t_a, t_m):
    cdef mpz_ptr a, m
    cdef int _a, _m
    a = fq12_t_get(&_a)
    m = fq_t_get(&_m)
    fq12_t_set_fq12(a, t_a)
    fq_t_set_pylong(m, t_m)
    fq12_t_mul_fq_t(a, a, m)
    res = fq12_t_get_fq12(a)
    fq12_t_release(_a)
    fq_t_release(_m)
    return res


def fq12_add(t_a, t_m):
    cdef mpz_ptr a, m
    cdef int _a, _m
    a = fq12_t_get(&_a)
    m = fq12_t_get(&_m)
    fq12_t_set_fq12(a, t_a)
    fq12_t_set_fq12(m, t_m)
    fq12_t_add(a, a, m)
    res = fq12_t_get_fq12(a)
    fq12_t_release(_a)
    fq12_t_release(_m)
    return res


def fq12_mul(t_a, t_m):
    cdef mpz_ptr a, m
    cdef int _a, _m
    a = fq12_t_get(&_a)
    m = fq12_t_get(&_m)
    fq12_t_set_fq12(a, t_a)
    fq12_t_set_fq12(m, t_m)
    fq12_t_mul(a, a, m)
    res = fq12_t_get_fq12(a)
    fq12_t_release(_a)
    fq12_t_release(_m)
    return res


# ec.py, paring.py calc operations
def fq2_to_affine(X, Y, Z, INF):
    cdef mpz_ptr x, y, z
    cdef int _x, _y, _z, inf
    x = fq2_t_get(&_x)
    y = fq2_t_get(&_y)
    z = fq2_t_get(&_z)
    fq2_t_set_fq2(x, X)
    fq2_t_set_fq2(y, Y)
    fq2_t_set_fq2(z, Z)
    inf = 1 if INF else 0

    fq2_t_to_affine(x, y, &inf, x, y, z, inf)

    rx = fq2_t_get_fq2(x)
    ry = fq2_t_get_fq2(y)
    rinf = True if inf else False

    fq2_t_release(_x)
    fq2_t_release(_y)
    fq2_t_release(_z)
    return rx, ry, rinf


def fq2_double_point(X, Y, INF):
    cdef mpz_ptr x, y
    cdef int _x, _y, inf
    x = fq2_t_get(&_x)
    y = fq2_t_get(&_y)
    fq2_t_set_fq2(x, X)
    fq2_t_set_fq2(y, Y)
    inf = 1 if INF else 0

    fq2_t_double_point(x, y, &inf, x, y, &inf)
    xr = fq2_t_get_fq2(x)
    yr = fq2_t_get_fq2(y)
    infr = True if inf else False

    fq2_t_release(_x)
    fq2_t_release(_y)
    return xr, yr, infr


def fq2_add_points(X1, Y1, INF1, X2, Y2, INF2):
    cdef mpz_ptr x1, y1, x2, y2
    cdef int _x1, _y1, _x2, _y2, inf1, inf2
    x1 = fq2_t_get(&_x1)
    y1 = fq2_t_get(&_y1)
    x2 = fq2_t_get(&_x2)
    y2 = fq2_t_get(&_y2)
    fq2_t_set_fq2(x1, X1)
    fq2_t_set_fq2(y1, Y1)
    inf1 = 1 if INF1 else 0
    fq2_t_set_fq2(x2, X2)
    fq2_t_set_fq2(y2, Y2)
    inf2 = 1 if INF2 else 0

    fq2_t_add_points(x1, y1, &inf1, x1, y1, &inf1, x2, y2, &inf2)
    xr = fq2_t_get_fq2(x1)
    yr = fq2_t_get_fq2(y1)
    infr = True if inf1 else False

    fq2_t_release(_x1)
    fq2_t_release(_y1)
    fq2_t_release(_x2)
    fq2_t_release(_y2)
    return xr, yr, infr


def fq_double_point_jacobian(X, Y, Z):
    cdef mpz_ptr x, y, z
    cdef int _x, _y, _z
    x = fq_t_get(&_x)
    y = fq_t_get(&_y)
    z = fq_t_get(&_z)

    fq_t_set_pylong(x, X)
    fq_t_set_pylong(y, Y)
    fq_t_set_pylong(z, Z)

    fq_t_double_point_jacobian(x, y, z, x, y, z)
    xr = fq_t_get_pylong(x)
    yr = fq_t_get_pylong(y)
    zr = fq_t_get_pylong(z)

    fq_t_release(_x)
    fq_t_release(_y)
    fq_t_release(_z)
    return xr, yr, zr


def fq2_double_point_jacobian(X, Y, Z):
    cdef mpz_ptr x, y, z
    cdef int _x, _y, _z
    x = fq2_t_get(&_x)
    y = fq2_t_get(&_y)
    z = fq2_t_get(&_z)

    fq2_t_set_fq2(x, X)
    fq2_t_set_fq2(y, Y)
    fq2_t_set_fq2(z, Z)

    fq2_t_double_point_jacobian(x, y, z, x, y, z)
    xr = fq2_t_get_fq2(x)
    yr = fq2_t_get_fq2(y)
    zr = fq2_t_get_fq2(z)

    fq2_t_release(_x)
    fq2_t_release(_y)
    fq2_t_release(_z)
    return xr, yr, zr


def fq12_double_point_jacobian(X, Y, Z):
    cdef mpz_ptr x, y, z
    cdef int _x, _y, _z
    x = fq12_t_get(&_x)
    y = fq12_t_get(&_y)
    z = fq12_t_get(&_z)

    fq12_t_set_fq12(x, X)
    fq12_t_set_fq12(y, Y)
    fq12_t_set_fq12(z, Z)

    fq12_t_double_point_jacobian(x, y, z, x, y, z)
    xr = fq12_t_get_fq12(x)
    yr = fq12_t_get_fq12(y)
    zr = fq12_t_get_fq12(z)

    fq12_t_release(_x)
    fq12_t_release(_y)
    fq12_t_release(_z)
    return xr, yr, zr


def fq_add_points_jacobian(X1, Y1, Z1, INF1, X2, Y2, Z2, INF2):
    cdef mpz_ptr x1, y1, z1, x2, y2, z2
    cdef int _x1, _y1, _z1, _x2, _y2, _z2
    cdef int inf1, inf2
    inf1 = 1 if INF1 else 0
    inf2 = 1 if INF2 else 0
    if inf1:
        return X2, Y2, Z2, INF2
    if inf2:
        return X1, Y1, Z1, INF1
    x1 = fq_t_get(&_x1)
    y1 = fq_t_get(&_y1)
    z1 = fq_t_get(&_z1)
    x2 = fq_t_get(&_x2)
    y2 = fq_t_get(&_y2)
    z2 = fq_t_get(&_z2)

    fq_t_set_pylong(x1, X1)
    fq_t_set_pylong(y1, Y1)
    fq_t_set_pylong(z1, Z1)
    fq_t_set_pylong(x2, X2)
    fq_t_set_pylong(y2, Y2)
    fq_t_set_pylong(z2, Z2)

    fq_t_add_points_jacobian(x1, y1, z1, &inf1,
                             x1, y1, z1, &inf1, x2, y2, z2, &inf2)
    xr = fq_t_get_pylong(x1)
    yr = fq_t_get_pylong(y1)
    zr = fq_t_get_pylong(z1)
    infr = True if inf1 else False

    fq_t_release(_x1)
    fq_t_release(_y1)
    fq_t_release(_z1)
    fq_t_release(_x2)
    fq_t_release(_y2)
    fq_t_release(_z2)
    return xr, yr, zr, infr


def fq2_add_points_jacobian(X1, Y1, Z1, INF1, X2, Y2, Z2, INF2):
    cdef mpz_ptr x1, y1, z1, x2, y2, z2
    cdef int _x1, _y1, _z1, _x2, _y2, _z2
    cdef int inf1, inf2
    inf1 = 1 if INF1 else 0
    inf2 = 1 if INF2 else 0
    if inf1:
        return X2, Y2, Z2, INF2
    if inf2:
        return X1, Y1, Z1, INF1
    x1 = fq2_t_get(&_x1)
    y1 = fq2_t_get(&_y1)
    z1 = fq2_t_get(&_z1)
    x2 = fq2_t_get(&_x2)
    y2 = fq2_t_get(&_y2)
    z2 = fq2_t_get(&_z2)

    fq2_t_set_fq2(x1, X1)
    fq2_t_set_fq2(y1, Y1)
    fq2_t_set_fq2(z1, Z1)
    fq2_t_set_fq2(x2, X2)
    fq2_t_set_fq2(y2, Y2)
    fq2_t_set_fq2(z2, Z2)

    fq2_t_add_points_jacobian(x1, y1, z1, &inf1,
                              x1, y1, z1, &inf1, x2, y2, z2, &inf2)
    xr = fq2_t_get_fq2(x1)
    yr = fq2_t_get_fq2(y1)
    zr = fq2_t_get_fq2(z1)
    infr = True if inf1 else False

    fq2_t_release(_x1)
    fq2_t_release(_y1)
    fq2_t_release(_z1)
    fq2_t_release(_x2)
    fq2_t_release(_y2)
    fq2_t_release(_z2)
    return xr, yr, zr, infr


def fq12_add_points_jacobian(X1, Y1, Z1, INF1, X2, Y2, Z2, INF2):
    cdef mpz_ptr x1, y1, z1, x2, y2, z2
    cdef int _x1, _y1, _z1, _x2, _y2, _z2
    cdef int inf1, inf2
    inf1 = 1 if INF1 else 0
    inf2 = 1 if INF2 else 0
    if inf1:
        return X2, Y2, Z2, INF2
    if inf2:
        return X1, Y1, Z1, INF1
    x1 = fq12_t_get(&_x1)
    y1 = fq12_t_get(&_y1)
    z1 = fq12_t_get(&_z1)
    x2 = fq12_t_get(&_x2)
    y2 = fq12_t_get(&_y2)
    z2 = fq12_t_get(&_z2)

    fq12_t_set_fq12(x1, X1)
    fq12_t_set_fq12(y1, Y1)
    fq12_t_set_fq12(z1, Z1)
    fq12_t_set_fq12(x2, X2)
    fq12_t_set_fq12(y2, Y2)
    fq12_t_set_fq12(z2, Z2)

    fq12_t_add_points_jacobian(x1, y1, z1, &inf1,
                               x1, y1, z1, &inf1, x2, y2, z2, &inf2)
    xr = fq12_t_get_fq12(x1)
    yr = fq12_t_get_fq12(y1)
    zr = fq12_t_get_fq12(z1)
    infr = True if inf1 else False

    fq12_t_release(_x1)
    fq12_t_release(_y1)
    fq12_t_release(_z1)
    fq12_t_release(_x2)
    fq12_t_release(_y2)
    fq12_t_release(_z2)
    return xr, yr, zr, infr


def fq2_scalar_mult_jacobian(C, X, Y, Z, INF):
    cdef mpz_ptr c, x, y, z
    cdef int _c, _x, _y, _z, inf
    c = fq_t_get(&_c)
    x = fq2_t_get(&_x)
    y = fq2_t_get(&_y)
    z = fq2_t_get(&_z)
    fq_t_set_pylong(c, C)
    fq2_t_set_fq2(x, X)
    fq2_t_set_fq2(y, Y)
    fq2_t_set_fq2(z, Z)
    inf = 1 if INF else 0

    fq2_t_scalar_mult_jacobian(x, y, z, &inf, c, x, y, z, &inf)
    xr = fq2_t_get_fq2(x)
    yr = fq2_t_get_fq2(y)
    zr = fq2_t_get_fq2(z)
    infr = True if inf else False

    fq_t_release(_c)
    fq2_t_release(_x)
    fq2_t_release(_y)
    fq2_t_release(_z)
    return xr, yr, zr, infr


def fq2_double_line_eval(RX, RY, PX, PY):
    cdef mpz_ptr rx, ry, px, py, rop
    cdef int _rx, _ry, _px, _py, _rop
    rx = fq2_t_get(&_rx)
    ry = fq2_t_get(&_ry)
    px = fq_t_get(&_px)
    py = fq_t_get(&_py)
    rop = fq12_t_get(&_rop)

    fq2_t_set_fq2(rx, RX)
    fq2_t_set_fq2(ry, RY)
    fq_t_set_pylong(px, PX)
    fq_t_set_pylong(py, PY)
    fq2_t_double_line_eval(rop, rx, ry, px, py)
    res = fq12_t_get_fq12(rop)

    fq2_t_release(_rx)
    fq2_t_release(_ry)
    fq_t_release(_px)
    fq_t_release(_py)
    fq12_t_release(_rop)
    return res


def fq2_add_line_eval(RX, RY, QX, QY, PX, PY):
    cdef mpz_ptr rx, ry, qx, qy, px, py, rop
    cdef int _rx, _ry, _qx, _qy, _px, _py, _rop
    rx = fq2_t_get(&_rx)
    ry = fq2_t_get(&_ry)
    qx = fq2_t_get(&_qx)
    qy = fq2_t_get(&_qy)
    px = fq_t_get(&_px)
    py = fq_t_get(&_py)
    rop = fq12_t_get(&_rop)

    fq2_t_set_fq2(rx, RX)
    fq2_t_set_fq2(ry, RY)
    fq2_t_set_fq2(qx, QX)
    fq2_t_set_fq2(qy, QY)
    fq_t_set_pylong(px, PX)
    fq_t_set_pylong(py, PY)
    fq2_t_add_line_eval(rop, rx, ry, qx, qy, px, py)
    res = fq12_t_get_fq12(rop)

    fq2_t_release(_rx)
    fq2_t_release(_ry)
    fq2_t_release(_qx)
    fq2_t_release(_qy)
    fq_t_release(_px)
    fq_t_release(_py)
    fq12_t_release(_rop)
    return res


def fq2_untwist(X, Y):
    cdef mpz_ptr x, y, ropx, ropy
    cdef int _x, _y, _ropx, _ropy
    x = fq2_t_get(&_x)
    y = fq2_t_get(&_y)
    ropx = fq12_t_get(&_ropx)
    ropy = fq12_t_get(&_ropy)
    fq2_t_set_fq2(x, X)
    fq2_t_set_fq2(y, Y)
    fq2_t_untwist(ropx, ropy, x, y)
    resx = fq12_t_get_fq12(ropx)
    resy = fq12_t_get_fq12(ropy)
    fq2_t_release(_x)
    fq2_t_release(_y)
    fq12_t_release(_ropx)
    fq12_t_release(_ropy)
    return resx, resy


def fq_miller_loop(PX, PY, PINF, QX, QY, QINF):
    cdef mpz_ptr px, py, qx, qy, rop
    cdef int _px, _py, _qx, _qy, _rop, pinf, qinf
    px = fq_t_get(&_px)
    py = fq_t_get(&_py)
    qx = fq2_t_get(&_qx)
    qy = fq2_t_get(&_qy)
    rop = fq12_t_get(&_rop)

    fq_t_set_pylong(px, PX)
    fq_t_set_pylong(py, PY)
    pinf = 1 if PINF else 0
    fq2_t_set_fq2(qx, QX)
    fq2_t_set_fq2(qy, QY)
    qinf = 1 if QINF else 0

    fq_t_miller_loop(rop, px, py, pinf, qx, qy, qinf)
    res = fq12_t_get_fq12(rop)

    fq_t_release(_px)
    fq_t_release(_py)
    fq2_t_release(_qx)
    fq2_t_release(_qy)
    fq12_t_release(_rop)
    return res


def fq12_final_exp(t_x):
    cdef mpz_ptr e
    cdef int _e
    e = fq12_t_get(&_e)
    fq12_t_set_fq12(e, t_x)
    fq12_t_final_exp(e, e)
    res = fq12_t_get_fq12(e)
    fq12_t_release(_e)
    return res


def fq_ate_pairing_multi(Ps, Qs):
    cdef int n
    cdef fq_t *ps_x
    cdef fq_t *ps_y
    cdef int *ps_inf
    cdef fq2_t *qs_x
    cdef fq2_t *qs_y
    cdef int *qs_inf

    cdef mpz_ptr prod
    cdef int _prod
    prod = fq12_t_get(&_prod)

    n = len(Qs)

    ps_x = <fq_t *>malloc(n * sizeof(fq_t))
    ps_y = <fq_t *>malloc(n * sizeof(fq_t))
    ps_inf = <int *>malloc(n * sizeof(int))

    qs_x = <fq2_t *>malloc(n * sizeof(fq2_t))
    qs_y = <fq2_t *>malloc(n * sizeof(fq2_t))
    qs_inf = <int *>malloc(n * sizeof(int))

    for i in range(n):
        fq_t_init(ps_x[i])
        fq_t_init(ps_y[i])
        fq2_t_init(qs_x[i])
        fq2_t_init(qs_y[i])

        mpz_set_pylong(ps_x[i], Ps[i][0])
        mpz_set_pylong(ps_y[i], Ps[i][1])
        ps_inf[i] = 1 if Ps[i][2] else 0

        fq2_t_set_fq2(qs_x[i], Qs[i][0])
        fq2_t_set_fq2(qs_y[i], Qs[i][1])
        qs_inf[i] = 1 if Qs[i][2] else 0

    fq_t_ate_pairing_multi(n, prod,
                           ps_x, ps_y, ps_inf,
                           qs_x, qs_y, qs_inf)
    res = fq12_t_get_fq12(prod)

    for i in range(n):
        mpz_clear(ps_x[i])
        mpz_clear(ps_y[i])
        mpz_clear(&qs_x[i][0])
        mpz_clear(&qs_x[i][1])
        mpz_clear(&qs_y[i][0])
        mpz_clear(&qs_y[i][1])

    free(ps_x)
    free(ps_y)
    free(ps_inf)
    free(qs_x)
    free(qs_y)
    free(qs_inf)

    fq12_t_release(_prod)
    return res


# fq queues vars and classes
DEF Q_BLOCK_SIZE = 100
DEF Q2_BLOCK_SIZE = 100
DEF Q6_BLOCK_SIZE = 100
DEF Q12_BLOCK_SIZE = 100


cdef fq_t fq_q[Q_BLOCK_SIZE]
cdef int fq_qi[Q_BLOCK_SIZE]
cdef int fq_qt  # tail
cdef fq2_t fq2_q[Q2_BLOCK_SIZE]
cdef int fq2_qi[Q2_BLOCK_SIZE]
cdef int fq2_qt  # tail
cdef fq6_t fq6_q[Q6_BLOCK_SIZE]
cdef int fq6_qi[Q6_BLOCK_SIZE]
cdef int fq6_qt  # tail
cdef fq12_t fq12_q[Q12_BLOCK_SIZE]
cdef int fq12_qi[Q12_BLOCK_SIZE]
cdef int fq12_qt


cdef mpz_ptr fq_t_get(int *pi):
    global fq_qt
    cdef mpz_ptr p
    if fq_qt <= 0:
        pi[0] = -1
        raise Exception('fq_q queue empty')
    pi[0] = fq_qi[fq_qt]
    p = fq_q[fq_qt]
    fq_qt -= 1
    return p


cdef void fq_t_release(int x):
    global fq_qt
    if x >= 0:
        fq_qt += 1
        fq_qi[fq_qt] = x


cdef void fq_t_alloc():
    global fq_qt
    cdef int i
    for i in range(Q_BLOCK_SIZE):
        fq_qt = i
        fq_qi[i] = i
        fq_t_init(fq_q[i])


cdef mpz_ptr fq2_t_get(int *pi):
    global fq2_qt
    cdef mpz_ptr p
    if fq2_qt <= 0:
        pi[0] = -1
        raise Exception('fq2_q queue empty')
    pi[0] = fq2_qi[fq2_qt]
    p = fq2_q[fq2_qt]
    fq2_qt -= 1
    return p


cdef void fq2_t_release(int x):
    global fq2_qt
    if x >= 0:
        fq2_qt += 1
        fq2_qi[fq2_qt] = x


cdef void fq2_t_alloc():
    global fq2_qt
    cdef int i
    for i in range(Q2_BLOCK_SIZE):
        fq2_qt = i
        fq2_qi[i] = i
        fq2_t_init(fq2_q[i])


cdef mpz_ptr fq6_t_get(int *pi):
    global fq6_qt
    cdef mpz_ptr p
    if fq6_qt <= 0:
        pi[0] = -1
        raise Exception('fq6_q queue empty')
    pi[0] = fq6_qi[fq6_qt]
    p = fq6_q[fq6_qt]
    fq6_qt -= 1
    return p


cdef void fq6_t_release(int x):
    global fq6_qt
    if x >= 0:
        fq6_qt += 1
        fq6_qi[fq6_qt] = x


cdef void fq6_t_alloc():
    global fq6_qt
    cdef int i
    for i in range(Q6_BLOCK_SIZE):
        fq6_qt = i
        fq6_qi[i] = i
        fq6_t_init(fq6_q[i])


cdef mpz_ptr fq12_t_get(int *pi):
    global fq12_qt
    cdef mpz_ptr p
    if fq12_qt <= 0:
        pi[0] = -1
        raise Exception('fq12_q queue empty')
    pi[0] = fq12_qi[fq12_qt]
    p = fq12_q[fq12_qt]
    fq12_qt -= 1
    return p


cdef void fq12_t_release(int x):
    global fq12_qt
    if x >= 0:
        fq12_qt += 1
        fq12_qi[fq12_qt] = x


cdef void fq12_t_alloc():
    global fq12_qt
    cdef int i
    for i in range(Q12_BLOCK_SIZE):
        fq12_qt = i
        fq12_qi[i] = i
        fq12_t_init(fq12_q[i])


# module init
def fq_queues_init():
    fq_t_alloc()
    fq2_t_alloc()
    fq6_t_alloc()
    fq12_t_alloc()


def params_init():
    global Q_int
    Q_int = int('0x1a0111ea397fe69a4b1ba7b6434bacd764774b'
                '84f38512bf6730d2a0f6b0f6241eabfffeb153ff'
                'ffb9feffffffffaaab', 16)
    mpz_init2(Q, INIT_BITS)
    mpz_set_pylong(Q, Q_int)
    mpz_init2(B, INIT_BITS)
    mpz_set_pylong(B, 0x4)
    mpz_init2(N, INIT_BITS)
    mpz_set_pylong(N, int('0x73eda753299d7d483339d80809a1d80553bda4'
                          '02fffe5bfeffffffff00000001', 16))
    mpz_init2(NX, INIT_BITS)
    mpz_set_pylong(NX, 0xd201000000010000)

    # twist/untwist inverse constants
    mpz_init2(tw1, INIT_BITS)
    mpz_set_pylong(tw1, int('0xd0088f51cbff34d258dd3db21a5d66'
                            'bb23ba5c279c2895fb39869507b587b1'
                            '20f55ffff58a9ffffdcff7fffffffd556', 16))
    mpz_init2(tw2, INIT_BITS)
    mpz_set_pylong(tw2, int('0xd0088f51cbff34d258dd3db21a5d66'
                            'bb23ba5c279c2895fb39869507b587b1'
                            '20f55ffff58a9ffffdcff7fffffffd555', 16))

    fq_t_init_set_pylong(fq_t_zero, 0)
    fq_t_init_set_pylong(fq_t_one, 1)
    fq_t_init_set_pylong_mod(fq2_t_root, Q, -1)
    fq2_t_init_set_fq_t(fq2_t_zero, fq_t_zero)
    fq2_t_init_set_fq_t(fq2_t_one, fq_t_one)
    fq2_t_init_set_fq2(fq6_t_root, (1, 1))
    fq6_t_init_set_fq_t(fq6_t_zero, fq_t_zero)
    fq6_t_init_set_fq_t(fq6_t_one, fq_t_one)
    fq6_t_init_set_fq6(fq12_t_root, (0, 0, 1, 0, 0, 0))
    fq12_t_init_set_fq_t(fq12_t_zero, fq_t_zero)
    fq12_t_init_set_fq_t(fq12_t_one, fq_t_one)

    fq_t_init(mpz_final_exp_e)
    mpz_pow_ui(mpz_final_exp_e, Q, 4)
    fq_t_init(mpz_n2)
    mpz_pow_ui(mpz_n2, Q, 2)
    mpz_sub(mpz_final_exp_e, mpz_final_exp_e, mpz_n2)
    mpz_add_ui(mpz_final_exp_e, mpz_final_exp_e, 1)
    mpz_fdiv_q(mpz_final_exp_e, mpz_final_exp_e, N)

    fq_t_set_pylong(mpz_n2, 2)
    fq_t_init_set_pylong(mpz_n0, 0)
    fq_t_init_set_pylong(mpz_n3, 3)
    fq_t_init_set_pylong(mpz_n4, 4)
    fq_t_init_set_pylong(mpz_n8, 8)


def frob_coeffs_init():
    fq2_t_init_set_fq2(fc_6[0][0],
                       (0,
                       int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                           '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                           'fd8bfd00000000aaac', 16)))
    fq2_t_init_set_fq2(fc_6[0][1],
                       (int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                            '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                            'fd8bfd00000000aaad', 16), 0))
    fq2_t_init_set_fq2(fc_6[1][0],
                       (int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                            '3be6f89688de17d813620a00022e01fffffffeff'
                            'fe', 16), 0))
    fq2_t_init_set_fq2(fc_6[1][1],
                       (int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                            '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                            'fd8bfd00000000aaac', 16), 0))
    fq2_t_init_set_fq2(fc_6[2][0], (0, 1))
    fq2_t_init_set_fq2(fc_6[2][1],
                       (int('0x1a0111ea397fe69a4b1ba7b6434bacd764774b'
                            '84f38512bf6730d2a0f6b0f6241eabfffeb153ff'
                            'ffb9feffffffffaaaa', 16), 0))
    fq2_t_init_set_fq2(fc_6[3][0],
                       (int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                            '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                            'fd8bfd00000000aaac', 16), 0))
    fq2_t_init_set_fq2(fc_6[3][1],
                       (int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                            '3be6f89688de17d813620a00022e01fffffffeff'
                            'fe', 16), 0))
    fq2_t_init_set_fq2(fc_6[4][0],
                       (0,
                       int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                           '3be6f89688de17d813620a00022e01fffffffeff'
                           'fe', 16)))
    fq2_t_init_set_fq2(fc_6[4][1],
                       (int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                            '3be6f89688de17d813620a00022e01fffffffeff'
                            'ff', 16), 0))


    fq6_t_init_set_fq6(fc_12[0],
                       (int('0x1904d3bf02bb0667c231beb4202c0d1f0fd603'
                            'fd3cbd5f4f7b2443d784bab9c4f67ea53d63e781'
                            '3d8d0775ed92235fb8', 16),
                        int('0xfc3e2b36c4e03288e9e902231f9fb854a14787'
                            'b6c7b36fec0c8ec971f63c5f282d5ac14d6c7ec2'
                            '2cf78a126ddc4af3', 16), 0, 0, 0, 0))
    fq6_t_init_set_fq6(fc_12[1],
                       (int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                            '3be6f89688de17d813620a00022e01fffffffeff'
                            'ff', 16), 0, 0, 0, 0, 0))
    fq6_t_init_set_fq6(fc_12[2],
                       (int('0x135203e60180a68ee2e9c448d77a2cd91c3ded'
                            'd930b1cf60ef396489f61eb45e304466cf3e67fa'
                            '0af1ee7b04121bdea2', 16),
                        int('0x6af0e0437ff400b6831e36d6bd17ffe48395da'
                            'bc2d3435e77f76e17009241c5ee67992f72ec05f'
                            '4c81084fbede3cc09', 16), 0, 0, 0, 0))
    fq6_t_init_set_fq6(fc_12[3],
                       (int('0x5f19672fdf76ce51ba69c6076a0f77eaddb3a9'
                            '3be6f89688de17d813620a00022e01fffffffeff'
                            'fe', 16), 0, 0, 0, 0, 0))
    fq6_t_init_set_fq6(fc_12[4],
                       (int('0x144e4211384586c16bd3ad4afa99cc9170df35'
                            '60e77982d0db45f3536814f0bd5871c1908bd478'
                            'cd1ee605167ff82995', 16),
                        int('0x5b2cfd9013a5fd8df47fa6b48b1e045f398162'
                            '40c0b8fee8beadf4d8e9c0566c63a3e6e257f873'
                            '29b18fae980078116', 16), 0, 0, 0, 0))
    fq6_t_init_set_fq6(fc_12[5],
                       (int('0x1a0111ea397fe69a4b1ba7b6434bacd764774b'
                            '84f38512bf6730d2a0f6b0f6241eabfffeb153ff'
                            'ffb9feffffffffaaaa', 16), 0, 0, 0, 0, 0))
    fq6_t_init_set_fq6(fc_12[6],
                       (int('0xfc3e2b36c4e03288e9e902231f9fb854a14787'
                            'b6c7b36fec0c8ec971f63c5f282d5ac14d6c7ec2'
                            '2cf78a126ddc4af3', 16),
                        int('0x1904d3bf02bb0667c231beb4202c0d1f0fd603'
                            'fd3cbd5f4f7b2443d784bab9c4f67ea53d63e781'
                            '3d8d0775ed92235fb8', 16), 0, 0, 0, 0))
    fq6_t_init_set_fq6(fc_12[7],
                       (int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                            '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                            'fd8bfd00000000aaac', 16), 0, 0, 0, 0, 0))
    fq6_t_init_set_fq6(fc_12[8],
                       (int('0x6af0e0437ff400b6831e36d6bd17ffe48395da'
                            'bc2d3435e77f76e17009241c5ee67992f72ec05f'
                            '4c81084fbede3cc09', 16),
                        int('0x135203e60180a68ee2e9c448d77a2cd91c3ded'
                            'd930b1cf60ef396489f61eb45e304466cf3e67fa'
                            '0af1ee7b04121bdea2', 16), 0, 0, 0, 0))
    fq6_t_init_set_fq6(fc_12[9],
                       (int('0x1a0111ea397fe699ec02408663d4de85aa0d85'
                            '7d89759ad4897d29650fb85f9b409427eb4f49ff'
                            'fd8bfd00000000aaad', 16), 0, 0, 0, 0, 0))
    fq6_t_init_set_fq6(fc_12[10],
                       (int('0x5b2cfd9013a5fd8df47fa6b48b1e045f398162'
                            '40c0b8fee8beadf4d8e9c0566c63a3e6e257f873'
                            '29b18fae980078116', 16),
                        int('0x144e4211384586c16bd3ad4afa99cc9170df35'
                            '60e77982d0db45f3536814f0bd5871c1908bd478'
                            'cd1ee605167ff82995', 16), 0, 0, 0, 0))


def module_init():
    fq_queues_init()
    params_init()
    frob_coeffs_init()


module_init()
