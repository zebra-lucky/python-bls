from collections import namedtuple
from copy import deepcopy

from . import bls12381
from .fields import (FieldExtBase, Fq, Fq2, Fq6, Fq12,
                     fq2_t_add_fq_int, fq2_t_sub_fq_int,
                     fq2_t_mul_fq2_t, fq2_t_mul_fq_int, fq_int_sub_fq2_t,
                     fq2_t_add_fq2_t, fq2_t_sub_fq2_t, fq6_t_mul_fq6_t,
                     fq6_t_mul_fq2_t, fq6_t_mul_fq_int, fq6_t_add_fq6_t,
                     fq6_t_sub_fq6_t, fq12_t_mul_fq12_t, fq12_t_mul_fq2_t,
                     fq12_t_mul_fq_int, fq12_t_add_fq12_t, fq12_t_sub_fq12_t)
from .util import hash256, hash512


# Struct for elliptic curve parameters
EC = namedtuple("EC", "q a b gx gy g2x g2y n h x k sqrt_n3 sqrt_n3m1o2")

# use secp256k1 as default
default_ec = EC(*bls12381.parameters())
default_ec_twist = EC(*bls12381.parameters_twist())


class AffinePoint:
    """
    Elliptic curve point, can represent any curve, and use Fq or Fq2
    coordinates.
    """
    def __init__(self, x, y, infinity, ec=default_ec):
        if (not isinstance(x, Fq) and not isinstance(x, FieldExtBase) or
           (not isinstance(y, Fq) and not isinstance(y, FieldExtBase)) or
           type(x) != type(y)):
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
        if self.infinity:
            return True
        left = self.y * self.y
        right = self.x * self.x * self.x + self.ec.a * self.x + self.ec.b

        return left == right

    def __add__(self, other):
        if other == 0:
            return self
        if not isinstance(other, AffinePoint):
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
        if not isinstance(other, AffinePoint):
            return False
        return (self.x == other.x and
                self.y == other.y and
                self.infinity == other.infinity)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __mul__(self, c):
        if not isinstance(c, (Fq, int)):
            raise ValueError("Error, must be int or Fq")
        return scalar_mult_jacobian(c, self.to_jacobian(), self.ec).to_affine()

    def negate(self):
        return AffinePoint(self.x, -self.y, self.infinity, self.ec)

    def __rmul__(self, c):
        return self.__mul__(c)

    def to_jacobian(self):
        return JacobianPoint(self.x, self.y, self.FE.one(self.ec.q),
                             self.infinity, self.ec)

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
    Elliptic curve point, can represent any curve, and use Fq or Fq2
    coordinates. Uses Jacobian coordinates so that point addition
    does not require slow inversion.
    """
    def __init__(self, x, y, z, infinity, ec=default_ec):
        if (not isinstance(x, Fq) and not isinstance(x, FieldExtBase) or
           (not isinstance(y, Fq) and not isinstance(y, FieldExtBase)) or
           (not isinstance(z, Fq) and not isinstance(z, FieldExtBase))):
            raise Exception("x,y should be field elements")
        self.FE = type(x)
        self.x = x
        self.y = y
        self.z = z
        self.infinity = infinity
        self.ec = ec

    def is_on_curve(self):
        if self.infinity:
            return True
        return self.to_affine().is_on_curve()

    def __add__(self, other):
        if other == 0:
            return self
        if not isinstance(other, JacobianPoint):
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
        if not isinstance(other, JacobianPoint):
            return False
        return self.to_affine() == other.to_affine()

    def __ne__(self, other):
        return not self.__eq__(other)

    def __mul__(self, c):
        if not isinstance(c, (Fq, int)):
            raise ValueError("Error, must be int or Fq")
        return scalar_mult_jacobian(c, self, self.ec)

    def __rmul__(self, c):
        return self.__mul__(c)

    def to_affine(self):
        if self.infinity:
            return AffinePoint(Fq.zero(self.ec.q), Fq.zero(self.ec.q),
                               self.infinity, self.ec)
        new_x = self.x / pow(self.z, 2)
        new_y = self.y / pow(self.z, 3)
        return AffinePoint(new_x, new_y, self.infinity, self.ec)

    def serialize(self):
        return self.to_affine().serialize()

    def __deepcopy__(self, memo):
        return JacobianPoint(deepcopy(self.x, memo),
                             deepcopy(self.y, memo),
                             deepcopy(self.z, memo),
                             self.infinity,
                             self.ec)


def y_for_x(x, ec=default_ec, FE=Fq):
    """
    Solves y = sqrt(x^3 + ax + b) for both valid ys
    """
    if not isinstance(x, FE):
        x = FE(ec.q, x)

    u = x * x * x + ec.a * x + ec.b

    y = u.modsqrt()
    if y == 0 or not AffinePoint(x, y, False, ec).is_on_curve():
        raise ValueError("No y for point x")
    return [y, ec.q - y]


def double_point(p1, ec=default_ec, FE=Fq):
    """
    Basic elliptic curve point doubling
    """
    x, y = p1.x, p1.y
    left = 3 * x * x
    left = left + ec.a
    s = left / (2 * y)
    new_x = s * s - x - x
    new_y = s * (x - new_x) - y
    return AffinePoint(new_x, new_y, False, ec)


def add_points(p1, p2, ec=default_ec, FE=Fq):
    """
    Basic elliptic curve point addition
    """
    assert(p1.is_on_curve())
    assert(p2.is_on_curve())
    if p1.infinity:
        return p2
    if p2.infinity:
        return p1
    if p1 == p2:
        return double_point(p1, ec, FE)
    if p1.x == p2.x:
        return AffinePoint(FE.zero(), FE.zero(), True, ec)

    x1, y1 = p1.x, p1.y
    x2, y2 = p2.x, p2.y
    s = (y2 - y1) / (x2 - x1)
    new_x = s * s - x1 - x2
    new_y = s * (x1 - new_x) - y1
    return AffinePoint(new_x, new_y, False, ec)


def double_point_jacobian(p1, ec=default_ec, FE=Fq):
    """
    Jacobian elliptic curve point doubling
    http://www.hyperelliptic.org/EFD/oldefd/jacobian.html
    """
    X, Y, Z = p1.x, p1.y, p1.z
    if Y == FE.zero(ec.q) or p1.infinity:
        return JacobianPoint(FE.one(ec.q), FE.one(ec.q),
                             FE.zero(ec.q), True, ec)
    Q = ec.q
    if FE == Fq2:
        xr, yr, zr = double_point_jacobian_fq2_t(X.ZT, Y.ZT, Z.ZT, ec)
        return JacobianPoint(Fq2(Q, xr), Fq2(Q, yr), Fq2(Q, zr), False, ec)
    elif FE == Fq and ec == default_ec:
        xr, yr, zr = double_point_jacobian_fq_int(X.Z, Y.Z, Z.Z, ec)
        return JacobianPoint(Fq(Q, xr), Fq(Q, yr), Fq(Q, zr), False, ec)
    elif FE == Fq and ec == default_ec_twist:
        xr, yr, zr = double_point_jacobian_fq_int_twist(X.Z, Y.Z, Z.Z, ec)
        return JacobianPoint(Fq2(Q, xr), Fq2(Q, yr), Fq2(Q, zr), False, ec)
    elif FE == Fq12:
        xr, yr, zr = double_point_jacobian_fq12_t(X.ZT, Y.ZT, Z.ZT, ec)
        return JacobianPoint(Fq12(Q, xr), Fq12(Q, yr), Fq12(Q, zr), False, ec)
    else:
        raise ValueError('FE must be Fq, Fq2 or Fq12')


def double_point_jacobian_fq_int(X, Y, Z, ec):
    '''dobule point with fq int X, Y, Z, returning tuple'''
    P = ec.q
    ec_a = ec.a.Z
    # S = 4*X*Y^2
    S = 4*X*Y*Y%P

    Z_sq = Z*Z%P
    Z_4th = Z_sq*Z_sq%P
    Y_sq = Y*Y%P
    Y_4th = Y_sq*Y_sq%P

    # M = 3*X^2 + a*Z^4
    M = (3*X*X%P + ec_a*Z_4th%P)%P

    # X' = M^2 - 2*S
    X_p = (M*M%P - 2*S%P)%P
    # Y' = M*(S - X') - 8*Y^4
    Y_p = (M*((S - X_p)%P)%P - 8*Y_4th%P)%P
    # Z' = 2*Y*Z
    Z_p = 2*Y*Z%P
    return X_p, Y_p, Z_p


def double_point_jacobian_fq_int_twist(X, Y, Z, ec):
    '''dobule point with fq int X, Y, Z, returning tuple'''
    P = ec.q
    ec_a = ec.a.ZT
    addi_f = fq2_t_add_fq_int
    subi_f = fq2_t_sub_fq_int
    muli_f = fq2_t_mul_fq_int
    mul_f = fq2_t_mul_fq2_t

    # S = 4*X*Y^2
    S = 4*X*Y*Y%P

    Z_sq = Z*Z%P
    Z_4th = Z_sq*Z_sq%P
    Y_sq = Y*Y%P
    Y_4th = Y_sq*Y_sq%P

    # M = 3*X^2 + a*Z^4
    M = addi_f(P, muli_f(P, ec_a, Z_4th), 3*X*X%P)

    # X' = M^2 - 2*S
    X_p = subi_f(P, mul_f(P, M, M), 2*S%P)
    # Y' = M*(S - X') - 8*Y^4
    Y_p = subi_f(P, mul_f(P, M, fq_int_sub_fq2_t(P, S, X_p)) , 8*Y_4th%P)
    # Z' = 2*Y*Z
    Z_p = addi_f(P, (0,0), 2*Y*Z%P)
    return X_p, Y_p, Z_p


def double_point_jacobian_fq2_t(X, Y, Z, ec):
    '''dobule point with fq2 tuples X, Y, Z, returning tuple of tuples'''
    if ec == default_ec_twist:
        mul_ec_a = fq2_t_mul_fq2_t
        ec_a = ec.a.ZT
    else:
        mul_ec_a = fq2_t_mul_fq_int
        ec_a = ec.a.Z
    func_t = (fq2_t_mul_fq2_t, fq2_t_mul_fq_int, mul_ec_a,
               fq2_t_add_fq2_t, fq2_t_sub_fq2_t)
    return double_point_jacobian_fqx_t(func_t, X, Y, Z, ec.q, ec_a)


def double_point_jacobian_fq12_t(X, Y, Z, ec):
    '''dobule point with fq12 tuples X, Y, Z, returning tuple of tuples'''
    if ec == default_ec_twist:
        mul_ec_a = fq12_t_mul_fq2_t
        ec_a = ec.a.ZT
    else:
        mul_ec_a = fq12_t_mul_fq_int
        ec_a = ec.a.Z
    func_t = (fq12_t_mul_fq12_t, fq12_t_mul_fq_int, mul_ec_a,
               fq12_t_add_fq12_t, fq12_t_sub_fq12_t)
    return double_point_jacobian_fqx_t(func_t, X, Y, Z, ec.q, ec_a)


def double_point_jacobian_fqx_t(func_t, X, Y, Z, P, ec_a):
    '''dobule point with tuples X, Y, Z, returning tuple of tuples,
    using func_t functions tuple for operations'''
    mul_f, mul_i_f, mul_ec_a, add_f, sub_f = func_t
    # S = 4*X*Y^2
    S = mul_f(P, mul_f(P, mul_i_f(P, X, 4), Y), Y)

    Z_sq = mul_f(P, Z, Z)
    Z_4th = mul_f(P, Z_sq, Z_sq)
    Y_sq = mul_f(P, Y, Y)
    Y_4th = mul_f(P, Y_sq, Y_sq)

    # M = 3*X^2 + a*Z^4
    M = mul_f(P, mul_i_f(P, X, 3), X)
    M = add_f(P, mul_ec_a(P, Z_4th, ec_a), M)

    # X' = M^2 - 2*S
    X_p = sub_f(P, mul_f(P, M, M), mul_i_f(P, S, 2))
    # Y' = M*(S - X') - 8*Y^4
    Y_p = sub_f(P, mul_f(P, M, sub_f(P, S, X_p)), mul_i_f(P, Y_4th, 8))
    # Z' = 2*Y*Z
    Z_p = mul_f(P, mul_i_f(P, Y, 2), Z)
    return X_p, Y_p, Z_p


def add_points_jacobian(p1, p2, ec=default_ec, FE=Fq):
    """
    Jacobian elliptic curve point addition
    http://www.hyperelliptic.org/EFD/oldefd/jacobian.html
    """
    if p1.infinity:
        return p2
    if p2.infinity:
        return p1

    if FE == Fq:
        U1, U2, S1, S2 = calc_u1_u2_s1_s2_fq_int(p1.x.Z, p1.y.Z, p1.z.Z,
                                                 p2.x.Z, p2.y.Z, p2.z.Z,
                                                 ec)
    elif FE == Fq2:
        U1, U2, S1, S2 = calc_u1_u2_s1_s2_fqx_t(fq2_t_mul_fq2_t,
                                                p1.x.ZT, p1.y.ZT, p1.z.ZT,
                                                p2.x.ZT, p2.y.ZT, p2.z.ZT,
                                                ec)
    elif FE == Fq12:
        U1, U2, S1, S2 = calc_u1_u2_s1_s2_fqx_t(fq12_t_mul_fq12_t,
                                                p1.x.ZT, p1.y.ZT, p1.z.ZT,
                                                p2.x.ZT, p2.y.ZT, p2.z.ZT,
                                                ec)
    else:
        raise ValueError('FE must be Fq, Fq2 or Fq12')

    if U1 == U2:
        if S1 != S2:
            return JacobianPoint(FE.one(ec.q), FE.one(ec.q),
                                 FE.zero(ec.q), True, ec)
        else:
            return double_point_jacobian(p1, ec, FE)

    type_u1 = type(U1)
    if type_u1 == int:
        return calc_jp_on_fq_int_us(U1, U2, S1, S2, p1.z.Z, p2.z.Z, ec)
    elif type_u1 == tuple and len(U1) == 2:
        func_t = (fq2_t_mul_fq2_t, fq2_t_sub_fq2_t, fq2_t_mul_fq_int)
        return calc_jp_on_fqx_t_us(func_t, U1, U2, S1, S2,
                                   p1.z.ZT, p2.z.ZT, ec)
    elif type_u1 == tuple and len(U1) == 12:
        func_t = (fq12_t_mul_fq12_t, fq12_t_sub_fq12_t, fq12_t_mul_fq_int)
        return calc_jp_on_fqx_t_us(func_t, U1, U2, S1, S2,
                                   p1.z.ZT, p2.z.ZT, ec)
    else:
        raise ValueError('FE must be Fq, Fq2 or Fq12')


def calc_u1_u2_s1_s2_fq_int(x1, y1, z1, x2, y2, z2, ec):
    '''x, y, z inputs of type int, returning tuple of int'''
    P = ec.q
    # U1 = X1*Z2^2
    U1 = x1*z2*z2%P
    # U2 = X2*Z1^2
    U2 = x2*z1*z1%P
    # S1 = Y1*Z2^3
    S1 = y1*z2*z2*z2%P
    # S2 = Y2*Z1^3
    S2 = y2*z1*z1*z1%P
    return(U1, U2, S1, S2)


def calc_u1_u2_s1_s2_fqx_t(mul_f, x1_t, y1_t, z1_t, x2_t, y2_t, z2_t, ec):
    '''x, y, z inputs of type fq2_t, returning tuple of fq2_t tuples'''
    P = ec.q
    # U1 = X1*Z2^2
    U1 = mul_f(P, mul_f(P, x1_t, z2_t), z2_t)
    # U2 = X2*Z1^2
    U2 = mul_f(P, mul_f(P, x2_t, z1_t), z1_t)
    # S1 = Y1*Z2^3
    S1 = mul_f(P, mul_f(P, mul_f(P, y1_t, z2_t), z2_t), z2_t)
    # S2 = Y2*Z1^3
    S2 = mul_f(P, mul_f(P, mul_f(P, y2_t, z1_t), z1_t), z1_t)
    return(U1, U2, S1, S2)


def calc_jp_on_fq_int_us(U1, U2, S1, S2, Z1, Z2, ec):
    '''calc jacobian point with int U1, U2, S1, S2, Z1, Z2'''
    P = ec.q
    # H = U2 - U1
    H = (U2-U1)%P
    # R = S2 - S1
    R = (S2-S1)%P
    H_sq = H*H%P
    H_cu = H*H_sq%P
    # X3 = R^2 - H^3 - 2*U1*H^2
    X3 = (R*R%P - H_cu - 2*U1*H_sq%P)%P
    # Y3 = R*(U1*H^2 - X3) - S1*H^3
    Y3 = (R*(U1*H_sq%P - X3)%P - S1*H_cu%P)%P
    # Z3 = H*Z1*Z2
    Z3 = H*Z1*Z2%P
    return JacobianPoint(Fq(P, X3), Fq(P, Y3), Fq(P, Z3), False, ec)


def calc_jp_on_fqx_t_us(func_t, U1, U2, S1, S2, Z1, Z2, ec):
    '''calc jacobian point with tuples U1, U2, S1, S2, Z1, Z2,
    using func_t functions tuple for operations'''
    P = ec.q
    mul_f, sub_f, muli_f = func_t
    # H = U2 - U1
    H = sub_f(P, U2, U1)
    # R = S2 - S1
    R = sub_f(P, S2, S1)
    H_sq = mul_f(P, H, H)
    H_cu = mul_f(P, H, H_sq)
    # X3 = R^2 - H^3 - 2*U1*H^2
    X3 = sub_f(P,
               sub_f(P, mul_f(P, R, R), H_cu),
               mul_f(P, H_sq, muli_f(P, U1, 2)))
    # Y3 = R*(U1*H^2 - X3) - S1*H^3
    Y3 = sub_f(P,
               mul_f(P, R, sub_f(P, mul_f(P, U1, H_sq), X3)),
               mul_f(P, S1, H_cu))
    # Z3 = H*Z1*Z2
    Z3 = mul_f(P, mul_f(P, Z1, Z2), H)
    if len(U1) == 2:
        return JacobianPoint(Fq2(P, X3), Fq2(P, Y3), Fq2(P, Z3), False, ec)
    else:
        return JacobianPoint(Fq12(P, X3), Fq12(P, Y3), Fq12(P, Z3), False, ec)


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
    if p1.infinity or c % ec.q == 0:
        return JacobianPoint(FE.one(ec.q), FE.one(ec.q),
                             FE.zero(ec.q), True, ec)
    if isinstance(c, FE):
        c = c.Z
    result = JacobianPoint(FE.one(ec.q), FE.one(ec.q),
                           FE.zero(ec.q), True, ec)
    addend = p1
    while c > 0:
        if c & 1:
            result += addend
        # double point
        addend += addend
        c = c >> 1
    return result


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
    f = Fq12.one(ec.q)
    wsq = Fq12(ec.q, f.root, Fq6.zero(ec.q))
    wcu = Fq12(ec.q, Fq6.zero(ec.q), f.root)
    return AffinePoint(point.x / wsq, point.y / wcu, False, ec)


def twist(point, ec=default_ec_twist):
    """
    Given an untwisted point, this converts it's
    coordinates to a point on the twisted curve. See Craig Costello
    book, look up twists.
    """
    f = Fq12.one(ec.q)
    wsq = Fq12(ec.q, f.root, Fq6.zero(ec.q))
    wcu = Fq12(ec.q, Fq6.zero(ec.q), f.root)
    new_x = (point.x * wsq)
    new_y = (point.y * wcu)
    return AffinePoint(new_x, new_y, False, ec)


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
