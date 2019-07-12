import logging
from collections import namedtuple
from copy import deepcopy

from . import bls12381
from .fields import Fq, Fq2, Fq12
from .util import hash256, hash512


# Struct for elliptic curve parameters
EC = namedtuple("EC", "q a b gx gy g2x g2y n h x k sqrt_n3 sqrt_n3m1o2")

# use secp256k1 as default
default_ec = EC(*bls12381.parameters())
default_ec_twist = EC(*bls12381.parameters_twist())


class AffinePoint:
    """
    Elliptic curve point, can represent only bls12-381 curve (a=0 or a=0,0
    for twisted curve), and use Fq, Fq2 or Fq12 coordinates.
    """
    def __init__(self, x, y, infinity, ec=default_ec):
        type_x = type(x)
        type_y = type(y)
        if (type_x != type_y
                or not type_x in (Fq, Fq2, Fq12)
                or not type_y in (Fq, Fq2, Fq12)):
            raise Exception("x,y should be field elements")
        self.FE = type(x)
        self.x = x
        self.y = y
        self.infinity = infinity
        self.ec = ec

    def is_on_curve(self):
        """
        Check that y^2 = x^3 + ax + b
        """
        return is_on_curve_affine(self, self.ec, self.FE)

    def __add__(self, other):
        if other == 0:
            return self
        if not type(other) == AffinePoint:
            raise Exception("Incorrect object")

        return add_points(self, other, self.ec, self.FE)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(other.negate())

    def __rsub__(self, other):
        return self.negate().__add__(other)

    def __str__(self):
        return ("AffinePoint(x=" + self.x.__str__() +
                ", y=" + self.y.__str__() +
                ", i=" + str(self.infinity) + ")\n")

    def __repr__(self):
        return ("AffinePoint(x=" + self.x.__repr__() +
                ", y=" + self.y.__repr__() +
                ", i=" + str(self.infinity) + ")\n")

    def __eq__(self, other):
        if not type(other) == AffinePoint:
            return False
        return (self.x == other.x and
                self.y == other.y and
                self.infinity == other.infinity)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __mul__(self, c):
        if not type(c) in (Fq, int):
            raise ValueError("Error, must be int or Fq")
        return scalar_mult_jacobian(c, self.to_jacobian(), self.ec).to_affine()

    def negate(self):
        return AffinePoint(self.x, -self.y, self.infinity, self.ec)

    def __rmul__(self, c):
        return self.__mul__(c)

    def to_jacobian(self):
        return to_jacobian(self, self.ec, self.FE)

    # Lexicographically greater than the negation
    def lex_gt_neg(self):
        if self.FE == Fq:
            if self.y > (self.ec.q // 2):
                return True
        elif self.FE == Fq2:
            if self.y[1] > (self.ec.q // 2):
                return True
        return False

    def serialize(self):
        output = bytearray(self.x.serialize())

        # If the y coordinate is the bigger one of the two, set the first
        # bit to 1.
        if self.lex_gt_neg():
            output[0] |= 0x80

        return bytes(output)

    def __deepcopy__(self, memo):
        return AffinePoint(deepcopy(self.x, memo),
                           deepcopy(self.y, memo),
                           self.infinity,
                           self.ec)


class JacobianPoint:
    """
    Elliptic curve point, can represent only bls12-381 curve (a=0 or a=0,0
    for twisted curve), and use Fq, Fq2 or Fq12 coordinates. Uses Jacobian
    coordinates so that point addition does not require slow inversion.
    """
    def __init__(self, x, y, z, infinity, ec=default_ec):
        type_x = type(x)
        type_y = type(y)
        type_z = type(z)
        if (type_x not in (Fq, Fq2, Fq12)
                or type_y not in (Fq, Fq2, Fq12)
                or type_z not in (Fq, Fq2, Fq12)):
            raise Exception("x,y should be field elements")
        self.FE = type(x)
        self.x = x
        self.y = y
        self.z = z
        self.infinity = infinity
        self.ec = ec

    def is_on_curve(self):
        return is_on_curve_jacobian(self, self.ec, self.FE)

    def __add__(self, other):
        if other == 0:
            return self
        if not type(other) == JacobianPoint:
            raise ValueError("Incorrect object")

        return add_points_jacobian(self, other, self.ec, self.FE)

    def __radd__(self, other):
        return self.__add__(other)

    def __str__(self):
        return ("JacobianPoint(x=" + self.x.__str__() +
                ", y=" + self.y.__str__() +
                "z=" + self.z.__str__() +
                ", i=" + str(self.infinity) + ")\n")

    def __eq__(self, other):
        if not type(other) == JacobianPoint:
            return False
        return self.to_affine() == other.to_affine()

    def __ne__(self, other):
        return not self.__eq__(other)

    def __mul__(self, c):
        if not type(c) in (Fq, int):
            raise ValueError("Error, must be int or Fq")
        return scalar_mult_jacobian(c, self, self.ec)

    def __rmul__(self, c):
        return self.__mul__(c)

    def to_affine(self):
        return to_affine(self, self.ec, self.FE)

    def serialize(self):
        return self.to_affine().serialize()

    def __deepcopy__(self, memo):
        return JacobianPoint(deepcopy(self.x, memo),
                             deepcopy(self.y, memo),
                             deepcopy(self.z, memo),
                             self.infinity,
                             self.ec)


def is_on_curve_affine(p, ec, FE):
    Q = ec.q
    px = p.x
    py = p.y
    pinf = p.infinity
    if FE == Fq:
        return fq_is_on_curve_affine(px.Z, py.Z, pinf)
    elif FE == Fq2:
        return fq2_is_on_curve_affine(px.ZT, py.ZT, pinf)
    elif FE == Fq12:
        return fq12_is_on_curve_affine(px.ZT, py.ZT, pinf)
    else:
        raise ValueError('FE must be Fq, Fq2 or Fq12')


def is_on_curve_jacobian(p, ec, FE):
    Q = ec.q
    px = p.x
    py = p.y
    pz = p.z
    pinf = p.infinity
    if FE == Fq:
        return fq_is_on_curve_jacobian(px.Z, py.Z, pz.Z, pinf)
    elif FE == Fq2:
        return fq2_is_on_curve_jacobian(px.ZT, py.ZT, pz.ZT, pinf)
    elif FE == Fq12:
        return fq12_is_on_curve_jacobian(px.ZT, py.ZT, pz.ZT, pinf)
    else:
        raise ValueError('FE must be Fq, Fq2 or Fq12')


def to_affine(p, ec, FE):
    Q = ec.q
    px = p.x
    py = p.y
    pz = p.z
    pinf = p.infinity
    if FE == Fq:
        xr, yr, infr = fq_to_affine(px.Z, py.Z, pz.Z, pinf)
    elif FE == Fq2:
        xr, yr, infr = fq2_to_affine(px.ZT, py.ZT, pz.ZT, pinf)
    elif FE == Fq12:
        xr, yr, infr = fq12_to_affine(px.ZT, py.ZT, pz.ZT, pinf)
    else:
        raise ValueError('FE must be Fq, Fq2 or Fq12')
    return AffinePoint(FE(Q, xr), FE(Q, yr), infr, ec)


def to_jacobian(p, ec, FE):
    Q = ec.q
    px = p.x
    py = p.y
    pinf = p.infinity
    if FE == Fq:
        xr, yr, zr, infr = fq_to_jacobian(px.Z, py.Z, pinf)
    elif FE == Fq2:
        xr, yr, zr, infr = fq2_to_jacobian(px.ZT, py.ZT, pinf)
    elif FE == Fq12:
        xr, yr, zr, infr = fq12_to_jacobian(px.ZT, py.ZT, pinf)
    else:
        raise ValueError('FE must be Fq, Fq2 or Fq12')
    return JacobianPoint(FE(Q, xr), FE(Q, yr), FE(Q, zr), infr, ec)


def y_for_x(x, ec=default_ec, FE=Fq):
    """
    Solves y = sqrt(x^3 + ax + b) for both valid ys
    """
    if not type(x) == FE:
        x = FE(ec.q, x)

    # u = x * x * x + ec.a * x + ec.b
    # a is 0 for bls12-381
    u = x * x * x + ec.b

    y = u.modsqrt()
    if y == 0 or not AffinePoint(x, y, False, ec).is_on_curve():
        raise ValueError("No y for point x")
    return [y, ec.q - y]


def double_point(p, ec, FE):
    """
    Basic elliptic curve point doubling
    """
    Q = ec.q
    px = p.x
    py = p.y
    pinf = p.infinity
    if FE == Fq:
        xr, yr, infr = fq_double_point(px.Z, py.Z, pinf)
    elif FE == Fq2:
        xr, yr, infr = fq2_double_point(px.ZT, py.ZT, pinf)
    elif FE == Fq12:
        xr, yr, infr = fq12_double_point(px.ZT, py.ZT, pinf)
    else:
        raise ValueError('FE must be Fq, Fq2 or Fq12', FE)
    return AffinePoint(FE(Q, xr), FE(Q, yr), infr, ec)


def add_points(p1, p2, ec=default_ec, FE=Fq):
    """
    Basic elliptic curve point addition
    """
    if p1.infinity:
        ec = p2.ec
    Q = ec.q
    x1 = p1.x
    y1 = p1.y
    inf1 = p1.infinity
    x2 = p2.x
    y2 = p2.y
    inf2 = p2.infinity
    if FE == Fq:
        xr, yr, infr = fq_add_points(x1.Z, y1.Z, inf1, x2.Z, y2.Z, inf2)
    elif FE == Fq2:
        xr, yr, infr = fq2_add_points(x1.ZT, y1.ZT, inf1, x2.ZT, y2.ZT, inf2)
    elif FE == Fq12:
        xr, yr, infr = fq12_add_points(x1.ZT, y1.ZT, inf1, x2.ZT, y2.ZT, inf2)
    else:
        raise ValueError('FE must be Fq, Fq2 or Fq12', FE)
    return AffinePoint(FE(Q, xr), FE(Q, yr), infr, ec)


def add_points_jacobian(p1, p2, ec=default_ec, FE=Fq):
    """
    Jacobian elliptic curve point addition
    http://www.hyperelliptic.org/EFD/oldefd/jacobian.html
    """
    if p1.infinity:
        ec = p2.ec
    Q = ec.q
    x1 = p1.x
    y1 = p1.y
    z1 = p1.z
    inf1 = p1.infinity
    x2 = p2.x
    y2 = p2.y
    z2 = p2.z
    inf2 = p2.infinity
    if FE == Fq:
        xr, yr, zr, infr = fq_add_points_jacobian(x1.Z, y1.Z, z1.Z, inf1,
                                                  x2.Z, y2.Z, z2.Z, inf2)
    elif FE == Fq2:
        xr, yr, zr, infr = fq2_add_points_jacobian(x1.ZT, y1.ZT, z1.ZT, inf1,
                                                   x2.ZT, y2.ZT, z2.ZT, inf2)
    elif FE == Fq12:
        xr, yr, zr, infr = fq12_add_points_jacobian(x1.ZT, y1.ZT, z1.ZT, inf1,
                                                    x2.ZT, y2.ZT, z2.ZT, inf2)
    else:
        raise ValueError('FE must be Fq, Fq2 or Fq12')
    return JacobianPoint(FE(Q, xr), FE(Q, yr), FE(Q, zr), infr, ec)


def scalar_mult(c, p1, ec=default_ec, FE=Fq):
    """
    Double and add:
    https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication
    """
    if p1.infinity or c % ec.q == 0:
        return AffinePoint(FE.zero(ec.q), FE.zero(ec.q), ec)
    result = AffinePoint(FE.zero(ec.q), FE.zero(ec.q), True, ec)
    addend = p1
    while c > 0:
        if c & 1:
            result += addend

        # double point
        addend += addend
        c = c >> 1

    return result


def scalar_mult_jacobian(c, p1, ec=default_ec, FE=Fq):
    """
    Double and add:
    https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication
    """
    Q = ec.q
    if type(c) == FE:
        c = c.Z
    px = p1.x
    py = p1.y
    pz = p1.z
    pinf = p1.infinity
    type_p = type(px)
    if type_p == Fq:
        xr, yr, zr, infr = fq_scalar_mult_jacobian(c,
                                                   px.Z, py.Z, pz.Z, pinf)
        return JacobianPoint(Fq(Q, xr), Fq(Q, yr), Fq(Q, zr), infr, ec)
    elif type_p == Fq2:
        xr, yr, zr, infr = fq2_scalar_mult_jacobian(c,
                                                    px.ZT, py.ZT, pz.ZT, pinf)
        return JacobianPoint(Fq2(Q, xr), Fq2(Q, yr), Fq2(Q, zr), infr, ec)
    elif type_p == Fq12:
        xr, yr, zr, infr = fq12_scalar_mult_jacobian(c,
                                                     px.ZT, py.ZT, pz.ZT, pinf)
        return JacobianPoint(Fq12(Q, xr), Fq12(Q, yr), Fq12(Q, zr), infr, ec)
    else:
        raise Exception("point should be Fq, Fq2 or Fq12 elements")


def generator_Fq(ec=default_ec):
    return AffinePoint(ec.gx, ec.gy, False, ec)


def generator_Fq2(ec=default_ec_twist):
    return AffinePoint(ec.g2x, ec.g2y, False, ec)


def untwist(point, ec=default_ec):
    """
    Given a point on G2 on the twisted curve, this converts it's
    coordinates back from Fq2 to Fq12. See Craig Costello book, look
    up twists.
    """
    Q = ec.q
    x = point.x
    y = point.y
    type_x = type(x)
    if type_x == Fq2:
        x_t, y_t = fq2_untwist(x.ZT, y.ZT)
    elif type_x == Fq12:
        x_t, y_t = fq12_untwist(x.ZT, y.ZT)
    else:
        raise Exception("point should be Fq2 or Fq12 elements")
    return AffinePoint(Fq12(Q, x_t), Fq12(Q, y_t), False, ec)


def twist(point, ec=default_ec_twist):
    """
    Given an untwisted point, this converts it's
    coordinates to a point on the twisted curve. See Craig Costello
    book, look up twists.
    """
    Q = ec.q
    x = point.x
    y = point.y
    type_x = type(x)
    if type_x == Fq2:
        x_t, y_t = fq2_twist(x.ZT, y.ZT)
    elif type_x == Fq12:
        x_t, y_t = fq12_twist(x.ZT, y.ZT)
    else:
        raise Exception("point should be Fq2 or Fq12 elements")
    return AffinePoint(Fq12(Q, x_t), Fq12(Q, y_t), False, ec)


def psi(P, ec=default_ec):
    ut = untwist(P, ec)
    t = AffinePoint(ut.x.qi_power(1), ut.y.qi_power(1), False, ec)
    t2 = twist(t, ec)
    return AffinePoint(t2.x[0][0], t2.y[0][0], False, ec)


# Rust pairings hashing https://github.com/zkcrypto/pairing
# https://www.di.ens.fr/~fouque/pub/latincrypt12.pdf
def sw_encode(t, ec=default_ec, FE=Fq):
    if t == 0:  # Maps t=0 to the point at infinity
        return JacobianPoint(FE.one(ec.q), FE.one(ec.q),
                             FE.zero(ec.q), True, ec)
    # Parity will ensure that sw_encode(t) = -sw_encode(-t)
    parity = False
    if FE == Fq:
        if t > -t:
            parity = True
    elif FE == Fq2:
        if t[1] > (-t)[1]:
            parity = True

    # w = t^2 + b + 1
    w = t * t + ec.b + 1

    # Map w=0 to generator and its negation
    if w == 0:
        if FE == Fq:
            ret = generator_Fq(ec)
        else:
            return generator_Fq2(ec)
        if parity:
            ret = ret.negate()
        return ret

    w = ~w * ec.sqrt_n3 * t

    # At least 1 of x1, x2, x3 is guaranteed to be a valid x
    # x1 = -wt + sqrt(-3)
    x1 = -w * t + ec.sqrt_n3m1o2

    # x2 = -x1 - 1
    x2 = FE.from_fq(ec.q, Fq(ec.q, -1)) - x1

    # x3 = 1/w^2 + 1
    x3 = ~(w * w) + 1

    # Constant time algorithm for finding the correct x, from
    # FT paper.
    Xx1 = 1
    Xx2 = 1
    try:
        y_for_x(x1, ec, FE)
    except:  # noqa: E772
        Xx1 = -1
    try:
        y_for_x(x2, ec, FE)
    except:  # noqa: E772
        Xx2 = -1

    index = (((Xx1 - 1) * Xx2) % 3)

    xs = [x1, x2, x3]
    ys = y_for_x(xs[index], ec, FE)
    ret = AffinePoint(xs[index], ys[0], False, ec)
    if ret.lex_gt_neg() is not parity:
        ret = ret.negate()
    return ret


# Performs a Fouque-Tibouchi hash
def hash_to_point_prehashed_Fq(m, ec=default_ec):
    if type(m) != bytes:
        m = m.encode("utf-8")
    t0 = Fq(ec.q, int.from_bytes(hash512(m + b"G1_0"), "big"))
    t1 = Fq(ec.q, int.from_bytes(hash512(m + b"G1_1"), "big"))

    P = sw_encode(t0, ec, Fq) + sw_encode(t1, ec, Fq)

    return P * ec.h  # Scaling by cofactor


def hash_to_point_Fq(m, ec=default_ec):
    h = hash256(m)
    return hash_to_point_prehashed_Fq(h, ec)


# Performs a Fouque-Tibouchi hash
def hash_to_point_prehashed_Fq2(m, ec=default_ec_twist):
    if type(m) != bytes:
        m = m.encode("utf-8")
    t0_0 = Fq(ec.q, int.from_bytes(hash512(m + b"G2_0_c0"), "big"))
    t0_1 = Fq(ec.q, int.from_bytes(hash512(m + b"G2_0_c1"), "big"))
    t1_0 = Fq(ec.q, int.from_bytes(hash512(m + b"G2_1_c0"), "big"))
    t1_1 = Fq(ec.q, int.from_bytes(hash512(m + b"G2_1_c1"), "big"))

    t0 = Fq2(ec.q, t0_0, t0_1)
    t1 = Fq2(ec.q, t1_0, t1_1)

    P = sw_encode(t0, ec, Fq2) + sw_encode(t1, ec, Fq2)

    # This is the cofactor multiplication, done in a more
    # efficient way. See page 11 of "Efficient hash maps
    # to G2 on BLS curves" by Budrioni and Pintore.
    x = -ec.x
    psi2P = psi(psi(2*P, ec), ec)
    t0 = x*P
    t1 = x*t0
    t2 = (t1 + t0) - P
    t3 = psi((x+1) * P, ec)
    return t2 - t3 + psi2P


def hash_to_point_Fq2(m, ec=default_ec_twist):
    h = hash256(m)
    return hash_to_point_prehashed_Fq2(h, ec)


try:
    from .fields_t import (
        fq_is_on_curve_affine, fq2_is_on_curve_affine, fq12_is_on_curve_affine,
        fq_is_on_curve_jacobian, fq2_is_on_curve_jacobian,
        fq12_is_on_curve_jacobian, fq_to_jacobian, fq2_to_jacobian,
        fq12_to_jacobian, fq_to_affine, fq2_to_affine, fq12_to_affine,
        fq_double_point, fq2_double_point, fq12_double_point,
        fq_add_points, fq2_add_points, fq12_add_points,
        fq2_twist, fq12_twist, fq2_untwist, fq12_untwist,
        fq_scalar_mult_jacobian, fq2_scalar_mult_jacobian,
        fq12_scalar_mult_jacobian, fq_add_points_jacobian,
        fq2_add_points_jacobian, fq12_add_points_jacobian,
        fq_double_point_jacobian, fq2_double_point_jacobian,
        fq12_double_point_jacobian
    )
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
