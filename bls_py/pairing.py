import logging
from collections import namedtuple

from . import bls12381
from .fields import Fq, Fq2, Fq12


# Struct for elliptic curve parameters
EC = namedtuple("EC", "q a b gx gy g2x g2y n h x k sqrt_n3 sqrt_n3m1o2")

# use secp256k1 as default
default_ec = EC(*bls12381.parameters())
default_ec_twist = EC(*bls12381.parameters())


def double_line_eval(R, P, ec=default_ec):
    """
    Creates an equation for a line tangent to R,
    and evaluates this at the point P. f(x) = y - sv - v.
    f(P).
    """
    Q = ec.q
    rx = R.x
    ry = R.y
    px = P.x
    py = P.y
    if type(rx) != Fq2 or type(px) != Fq:
        raise Exception("invalid elements")
    return Fq12(Q, fq2_double_line_eval(rx.ZT, ry.ZT, px.Z, py.Z))


def add_line_eval(R, Q, P, ec=default_ec):
    """
    Creates an equation for a line between R and Q,
    and evaluates this at the point P. f(x) = y - sv - v.
    f(P).
    """
    rx = R.x
    ry = R.y
    qx = Q.x
    qy = Q.y
    px = P.x
    py = P.y
    if type(rx) != Fq2 or type(qx) != Fq2 or type(px) != Fq:
        raise Exception("invalid elements")
    Q = ec.q
    return Fq12(Q, fq2_add_line_eval(rx.ZT, ry.ZT, qx.ZT, qy.ZT,
                                     px.Z, py.Z))


def miller_loop(P, Q, ec=default_ec):
    """
    Performs a double and add algorithm for the ate pairing. This algorithm
    is taken from Craig Costello's "Pairing for Beginners".
    """
    px = P.x
    py = P.y
    pinf = P.infinity
    qx = Q.x
    qy = Q.y
    qinf = Q.infinity
    if type(px) != Fq or type(qx) != Fq2:
        raise Exception("invalid elements")
    Q = ec.q
    return Fq12(Q, fq_miller_loop(px.Z, py.Z, pinf, qx.ZT, qy.ZT, qinf))


def final_exponentiation(element, ec):
    """
    Performs a final exponentiation on bls12-381 curve to map
    the result of the miller loop to a unique element of Fq12.
    """
    return Fq12(ec.q, fq12_final_exp(element.ZT))


def ate_pairing(P, Q, ec=default_ec):
    """
    Performs one ate pairing.
    """
    element = miller_loop(P, Q, ec)
    return final_exponentiation(element, ec)


def ate_pairing_multi(Ps, Qs, ec=default_ec):
    """
    Computes multiple pairings at once. This is more efficient,
    since we can multiply all the results of the miller loops,
    and perform just one final exponentiation.
    """
    Ps = tuple((ps.x.Z, ps.y.Z, ps.infinity) for ps in Ps)
    Qs = tuple((qs.x.ZT, qs.y.ZT, qs.infinity) for qs in Qs)
    return Fq12(ec.q, fq_ate_pairing_multi(Ps, Qs))


try:
    from .fields_t import (fq12_final_exp, fq2_double_line_eval,
                           fq2_add_line_eval, fq_miller_loop,
                           fq_ate_pairing_multi)
except ImportError as e:
    logging.error(f'Can not import from fields_t_c: {e}')


"""
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
"""
