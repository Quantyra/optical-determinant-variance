"""Exact real arithmetic and optical transformations used by the two checks."""


from dataclasses import dataclass
from fractions import Fraction as F
import math


@dataclass(frozen=True)
class Q2:
    a: F = F(0)
    b: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, "a", F(self.a))
        object.__setattr__(self, "b", F(self.b))

    @staticmethod
    def coerce(v):
        return v if isinstance(v, Q2) else Q2(v)

    def __add__(self, v):
        v = self.coerce(v)
        return Q2(self.a + v.a, self.b + v.b)

    __radd__ = __add__

    def __neg__(self):
        return Q2(-self.a, -self.b)

    def __sub__(self, v):
        return self + -self.coerce(v)

    def __rsub__(self, v):
        return self.coerce(v) + -self

    def __mul__(self, v):
        v = self.coerce(v)
        return Q2(self.a * v.a + 2 * self.b * v.b,
                  self.a * v.b + self.b * v.a)

    __rmul__ = __mul__

    def __truediv__(self, v):
        v = self.coerce(v)
        d = v.a * v.a - 2 * v.b * v.b
        if not d:
            raise ZeroDivisionError
        return self * Q2(v.a / d, -v.b / d)

    def sign(self):
        if not self.a:
            return (self.b > 0) - (self.b < 0)
        if not self.b or self.a * self.b > 0:
            return (self.a > 0) - (self.a < 0)
        dominant = self.a if self.a * self.a > 2 * self.b * self.b else self.b
        return (dominant > 0) - (dominant < 0)

    def __abs__(self):
        return self * self.sign()

    def __float__(self):
        return float(self.a) + math.sqrt(2) * float(self.b)


ZERO, ONE, INV_SQRT2 = Q2(), Q2(1), Q2(0, F(1, 2))
SQRT2 = 2 * INV_SQRT2
A = 1 - SQRT2


def beamsplitter(fields, i, j):
    x, y = fields[i], fields[j]
    fields[i], fields[j] = INV_SQRT2 * (x + y), INV_SQRT2 * (x - y)


def rotate(fields, i, j, sine):
    x, y = fields[i], fields[j]
    cosine = F(12, 13)
    fields[i], fields[j] = cosine * x - sine * y, sine * x + cosine * y


def optical_column(occupied, source_column, phase=-1):
    # Rails 0,1 are the target. Rails 2,3 and 4,5 are the control pairs.
    # The four successive controlled signs consume heralds (6,7), ..., (12,13).
    fields = [ZERO] * 14
    fields[occupied[source_column]] = ONE

    def controlled_sign(leaf_one, herald):
        beamsplitter(fields, 1, leaf_one)
        for rail, detector in ((1, herald), (leaf_one, herald + 1)):
            z, u = fields[rail], fields[detector]
            fields[rail], fields[detector] = A * z + u, SQRT2 * z + u
        beamsplitter(fields, 1, leaf_one)

    rotate(fields, 4, 5, F(5, 13))
    controlled_sign(3, 6)
    beamsplitter(fields, 0, 1)
    controlled_sign(5, 8)
    beamsplitter(fields, 0, 1)
    fields[1] = phase * fields[1]
    beamsplitter(fields, 0, 1)
    controlled_sign(5, 10)
    beamsplitter(fields, 0, 1)
    controlled_sign(3, 12)
    rotate(fields, 4, 5, F(-5, 13))
    return [fields[i] for i in occupied]


def components(matrix):
    supports = [{j for j, entry in enumerate(row) if entry != ZERO} for row in matrix]
    pending = set(range(len(matrix)))
    while pending:
        rows = {pending.pop()}
        columns = set()
        while True:
            columns.update(*(supports[i] for i in rows))
            more = {i for i in pending if supports[i] & columns}
            if not more:
                break
            pending -= more
            rows |= more
        assert len(rows) == len(columns)
        yield sorted(rows), sorted(columns)


def multiply(x, y):
    # Integer coordinates in Q(sqrt(2)) after common denominators are removed.
    a, b = x
    c, d = y
    return a*c + 2*b*d, a*d + b*c


def permanent(matrix):
    # Independent subset recurrence, with no signs or random estimator.
    totals = {0: ONE}
    for row in matrix:
        following = {}
        for used, value in totals.items():
            for j, entry in enumerate(row):
                bit = 1 << j
                if entry != ZERO and not used & bit:
                    key = used | bit
                    following[key] = following.get(key, ZERO) + value * entry
        totals = following
    return totals[(1 << len(matrix))-1]
