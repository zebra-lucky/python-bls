# cython: language_level=3
# -*- coding: utf-8 -*-

from .types cimport *


cdef mp_bitcnt_t INIT_BITS = 320


ctypedef mpz_t fq_t
ctypedef __mpz_struct fq2_t[2]
ctypedef __mpz_struct fq6_t[6]
ctypedef __mpz_struct fq12_t[12]


# fq queues
cdef mpz_ptr fq_t_get(int *pi)
cdef void fq_t_release(int x)
cdef void fq_t_alloc()

cdef mpz_ptr fq2_t_get(int *pi)
cdef void fq2_t_release(int x)
cdef void fq2_t_alloc()

cdef mpz_ptr fq6_t_get(int *pi)
cdef void fq6_t_release(int x)
cdef void fq6_t_alloc()

cdef mpz_ptr fq12_t_get(int *pi)
cdef void fq12_t_release(int x)
cdef void fq12_t_alloc()


# curve q, n, and other constants
cdef mpz_t Q, B, N, NX, tw1, tw2
cdef mpz_t mpz_final_exp_e
cdef fq_t fq_t_zero, fq_t_one, fq2_t_root
cdef fq2_t fq2_t_zero, fq2_t_one, fq6_t_root
cdef fq6_t fq6_t_zero, fq6_t_one, fq12_t_root
cdef fq12_t fq12_t_zero, fq12_t_one

# mpz_t 2, 3, 4, 8
cdef mpz_t mpz_n0, mpz_n2, mpz_n3, mpz_n4, mpz_n8

# frob coeffs
cdef fq2_t fc_6[5][2]
cdef fq6_t fc_12[11]


cdef str fq_t_get_pystr(fq_t x, base)
cdef str fq2_t_get_pystr(fq2_t x, base)
cdef str fq6_t_get_pystr(fq6_t x, base)
cdef str fq12_t_get_pystr(fq12_t x, base)
