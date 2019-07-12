# -*- coding: utf-8 -*-

from .fields_t import (fq_invert, fq_floordiv, fq_pow,
                       fq2_neg, fq2_invert, fq2_floordiv,
                       fq2_pow, fq2_qi_pow, fq2_mul_by_nonresidue,
                       fq2_add, fq2_mul,
                       fq2_add_fq, fq2_sub, fq2_sub_fq, fq_sub_fq2,
                       fq2_mul_fq, fq12_mul_fq2, fq12_invert,
                       fq6_add, fq6_mul, fq6_floordiv,
                       fq6_mul_fq2, fq6_neg, fq6_invert, fq6_pow, fq6_qi_pow,
                       fq6_mul_by_nonresidue, fq6_add_fq2,
                       fq6_add_fq, fq6_sub, fq6_sub_fq2, fq6_sub_fq,
                       fq2_sub_fq6, fq_sub_fq6, fq6_mul_fq, fq12_floordiv,
                       fq12_pow, fq12_mul_fq, fq12_add, fq12_mul,
                       fq12_mul_fq6, fq12_neg, fq12_qi_pow,
                       fq12_add_fq, fq12_add_fq6, fq12_add_fq2,
                       fq12_sub, fq12_sub_fq, fq12_sub_fq6, fq12_sub_fq2,
                       fq_sub_fq12, fq6_sub_fq12, fq2_sub_fq12)


bls12381_q = int('0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf'
                 '6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab', 16)


FQ2_ONE_TUPLE = (1, 0)
FQ2_ZERO_TUPLE = (0, 0)
FQ6_ROOT_TUPLE = (1, 1)
FQ6_ONE_TUPLE = (1, 0, 0, 0, 0, 0)
FQ6_ZERO_TUPLE = (0, 0, 0, 0, 0, 0)
FQ12_ROOT_TUPLE = (0, 0, 1, 0, 0, 0)
FQ12_ONE_TUPLE = (1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
FQ12_ZERO_TUPLE = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


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
        return Fq(self.Q, fq_invert(self.Q, self.Z))

    def __pow__(self, X):
        return Fq(self.Q, fq_pow(self.Q, self.Z, X))

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
            return Fq(self.Q, fq_floordiv(self.Q, self.Z, X.Z))
        elif tx == int:
            return Fq(self.Q, fq_floordiv(self.Q, self.Z, X))
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
        return Fq2(self.Q, fq2_neg(self.ZT))

    def __invert__(self):
        return Fq2(self.Q, fq2_invert(self.ZT))

    def __pow__(self, e):
        return Fq2(self.Q, fq2_pow(self.ZT, e))

    def qi_power(self, i):
        global bls12381_q
        if self.Q != bls12381_q:
            raise NotImplementedError
        i %= 2
        if i == 0:
            return self
        return Fq2(self.Q, fq2_qi_pow(self.ZT, i))

    def mul_by_nonresidue(self):
        # multiply by u + 1
        return Fq2(self.Q, fq2_mul_by_nonresidue(self.ZT))

    def __add__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq2:
            return Fq2(Q, fq2_add(self.ZT, other.ZT))
        elif tx == int:
            return Fq2(Q, fq2_add_fq(self.ZT, other))
        elif tx == Fq:
            return Fq2(Q, fq2_add_fq(self.ZT, other.Z))
        else:
            return NotImplemented

    def __sub__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq2:
            return Fq2(Q, fq2_sub(self.ZT, other.ZT))
        elif tx == int:
            return Fq2(Q, fq2_sub_fq(self.ZT, other))
        elif tx == Fq:
            return Fq2(Q, fq2_sub_fq(self.ZT, other.Z))
        else:
            return NotImplemented

    def __rsub__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq2:
            return Fq2(Q, fq2_sub(other.ZT, self.ZT))
        elif tx == int:
            return Fq2(Q, fq_sub_fq2(other, self.ZT))
        elif tx == Fq:
            return Fq2(Q, fq_sub_fq2(other.Z, self.ZT))
        else:
            return NotImplemented

    def __mul__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq2:
            return Fq2(Q, fq2_mul(self.ZT, other.ZT))
        elif tx == int:
            return Fq2(Q, fq2_mul_fq(self.ZT, other))
        elif tx == Fq:
            return Fq2(Q, fq2_mul_fq(self.ZT, other.Z))
        else:
            return NotImplemented

    def __floordiv__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq12:
            return Fq12(Q, fq12_mul_fq2(fq12_invert(other.ZT), self.ZT))
        elif tx == Fq6:
            return Fq6(Q, fq6_mul_fq2(fq2_invert(other.ZT), self.ZT))
        elif tx == Fq2:
            return Fq2(Q, fq2_floordiv(self.ZT, other.ZT))
        elif tx == int:
            return Fq2(Q, fq2_mul_fq(self.ZT, fq_invert(Q, other)))
        elif tx == Fq:
            return Fq2(Q, fq2_mul_fq(self.ZT, fq_invert(Q, other.Z)))
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
        if type(a) != Fq2 or type(b) != Fq2 or type(c) != Fq2:
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
        return Fq6(self.Q, fq6_neg(self.ZT))

    def __invert__(self):
        return Fq6(self.Q, fq6_invert(self.ZT))

    def __pow__(self, e):
        return Fq6(self.Q, fq6_pow(self.ZT, e))

    def qi_power(self, i):
        global bls12381_q
        if self.Q != bls12381_q:
            raise NotImplementedError
        i %= 6
        if i == 0:
            return self
        return Fq6(self.Q, fq6_qi_pow(self.ZT, i))

    def mul_by_nonresidue(self):
        # multiply by v
        return Fq6(self.Q, fq6_mul_by_nonresidue(self.ZT))

    def __add__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq6:
            return Fq6(Q, fq6_add(self.ZT, other.ZT))
        elif tx == Fq2:
            return Fq6(Q, fq6_add_fq2(self.ZT, other.ZT))
        elif tx == Fq:
            return Fq6(Q, fq6_add_fq(self.ZT, other.Z))
        elif tx == int:
            return Fq6(Q, fq6_add_fq(self.ZT, other))
        else:
            return NotImplemented

    def __sub__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq6:
            return Fq6(Q, fq6_sub(self.ZT, other.ZT))
        elif tx == Fq2:
            return Fq6(Q, fq6_sub_fq2(self.ZT, other.ZT))
        elif tx == Fq:
            return Fq6(Q, fq6_sub_fq(self.ZT, other.Z))
        elif tx == int:
            return Fq6(Q, fq6_sub_fq(self.ZT, other))
        else:
            return NotImplemented

    def __rsub__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq6:
            return Fq6(Q, fq6_sub(other.ZT, self.ZT))
        elif tx == Fq2:
            return Fq6(Q, fq2_sub_fq6(other.ZT, self.ZT))
        elif tx == Fq:
            return Fq6(Q, fq_sub_fq6(other.Z, self.ZT))
        elif tx == int:
            return Fq6(Q, fq_sub_fq6(other, self.ZT))
        else:
            return NotImplemented

    def __mul__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq6:
            return Fq6(Q, fq6_mul(self.ZT, other.ZT))
        elif tx == Fq2:
            return Fq6(Q, fq6_mul_fq2(self.ZT, other.ZT))
        elif tx == Fq:
            return Fq6(Q, fq6_mul_fq(self.ZT, other.Z))
        elif tx == int:
            return Fq6(Q, fq6_mul_fq(self.ZT, other))
        else:
            return NotImplemented

    def __floordiv__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq12:
            return Fq12(Q, fq12_mul_fq6(fq12_invert(other.ZT), self.ZT))
        elif tx == Fq6:
            return Fq6(Q, fq6_floordiv(self.ZT, other.ZT))
        elif tx == Fq2:
            return Fq6(Q, fq6_mul_fq2(self.ZT, fq2_invert(other.ZT)))
        elif tx == Fq:
            return Fq6(Q, fq6_mul_fq(self.ZT, fq_invert(Q, other.Z)))
        elif tx == int:
            return Fq6(Q, fq6_mul_fq(self.ZT, fq_invert(Q, other)))
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
        return Fq12(self.Q, fq12_neg(self.ZT))

    def __invert__(self):
        return Fq12(self.Q, fq12_invert(self.ZT))

    def __pow__(self, e):
        return Fq12(self.Q, fq12_pow(self.ZT, e))

    def qi_power(self, i):
        global bls12381_q
        if self.Q != bls12381_q:
            raise NotImplementedError
        i %= 12
        if i == 0:
            return self
        return Fq12(self.Q, fq12_qi_pow(self.ZT, i))

    def __add__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq12:
            return Fq12(Q, fq12_add(self.ZT, other.ZT))
        elif tx == Fq:
            return Fq12(Q, fq12_add_fq(self.ZT, other.Z))
        elif tx == Fq6:
            return Fq12(Q, fq12_add_fq6(self.ZT, other.ZT))
        elif tx == Fq2:
            return Fq12(Q, fq12_add_fq2(self.ZT, other.ZT))
        elif tx == int:
            return Fq12(Q, fq12_add_fq(self.ZT, other))
        else:
            return NotImplemented

    def __sub__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq12:
            return Fq12(Q, fq12_sub(self.ZT, other.ZT))
        elif tx == Fq:
            return Fq12(Q, fq12_sub_fq(self.ZT, other.Z))
        elif tx == Fq6:
            return Fq12(Q, fq12_sub_fq6(self.ZT, other.ZT))
        elif tx == Fq2:
            return Fq12(Q, fq12_sub_fq2(self.ZT, other.ZT))
        elif tx == int:
            return Fq12(Q, fq12_sub_fq(self.ZT, other))
        else:
            return NotImplemented

    def __rsub__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq12:
            return Fq12(Q, fq12_sub(other.ZT, self.ZT))
        elif tx == Fq:
            return Fq12(Q, fq_sub_fq12(other.Z, self.ZT))
        elif tx == Fq6:
            return Fq12(Q, fq6_sub_fq12(other.ZT, self.ZT))
        elif tx == Fq2:
            return Fq12(Q, fq2_sub_fq12(other.ZT, self.ZT))
        elif tx == int:
            return Fq12(Q, fq_sub_fq12(other, self.ZT))
        else:
            return NotImplemented

    def __mul__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq12:
            return Fq12(Q, fq12_mul(self.ZT, other.ZT))
        elif tx == int:
            return Fq12(Q, fq12_mul_fq(self.ZT, other))
        elif tx == Fq:
            return Fq12(Q, fq12_mul_fq(self.ZT, other.Z))
        elif tx == Fq2:
            return Fq12(Q, fq12_mul_fq2(self.ZT, other.ZT))
        elif tx == Fq6:
            return Fq12(Q, fq12_mul_fq6(self.ZT, other.ZT))
        else:
            return NotImplemented

    def __floordiv__(self, other):
        Q = self.Q
        tx = type(other)
        if tx == Fq12:
            return Fq12(Q, fq12_floordiv(self.ZT, other.ZT))
        elif tx == int:
            return Fq12(Q, fq12_mul(self.ZT, fq_invert(Q, other)))
        elif tx == Fq:
            return Fq12(Q, fq12_mul(self.ZT, fq_invert(Q, other.Z)))
        elif tx == Fq2:
            return Fq12(Q, fq12_mul_fq2(self.ZT, fq2_invert(other.ZT)))
        elif tx == Fq6:
            return Fq12(Q, fq12_mul_fq6(self.ZT, fq6_invert(other.ZT)))
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
