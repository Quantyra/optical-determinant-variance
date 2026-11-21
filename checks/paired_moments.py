"""Exact shared-sign second moments. Integer pairs represent Q(sqrt(2)) numerators."""

from collections import defaultdict
from functools import cache
from itertools import combinations, product
from math import lcm

from optical import F, ONE, Q2, ZERO, optical_column, permanent


def fixture(b, z):
    ports = [b, 3, 4, *range(6, 14)]
    return list(zip(*(optical_column(ports, j, z) for j in range(len(ports)))))


def row_weights(a, d):
    # The four copies are P, P, Q(x), Q(y). Only even edge multiplicities survive.
    mass = [abs(u) if (abs(u)-abs(v)).sign() >= 0 else abs(v) for u,v in zip(a,d)]
    support = [j for j, value in enumerate(mass) if value != ZERO]
    raw = {}
    for j in support:
        for k in support:
            pairs = [(j,k),(k,j)] if j != k else [(v,v) for v in support]
            for u,v in pairs:
                # sqrt(m_j m_k / m_u m_v) is one for distinct P columns.
                scale = ONE if j != k else mass[j] / mass[u]
                raw[j,k,u,v] = [scale*t for t in
                                  (a[u]*a[v], d[u]*a[v], a[u]*d[v], d[u]*d[v])]
    den = lcm(*(z for coeff in raw.values() for h in coeff
                for z in (h.a.denominator, h.b.denominator)))
    return den, {key: tuple(x for h in coeff for x in (int(h.a*den), int(h.b*den)))
                 for key, coeff in raw.items()}


def add_product(out, state, left, right, sign):
    # Multiply in Z[sqrt(2),x,y]/(x^2,y^2). Store 1,x,y,xy in that order.
    dest = out[state]
    for i in range(4):
        a,b = left[2*i:2*i+2]
        if a == b == 0:
            continue
        for j in range(4):
            if i & j:
                continue
            c,d = right[2*j:2*j+2]
            k = 2*(i|j)
            dest[k] += sign*(a*c+2*b*d)
            dest[k+1] += sign*(a*d+b*c)


def second_moment(base, slope, order):
    n = len(base)
    full = (1<<n)-1
    rows = [row_weights(base[i], slope[i]) for i in order]
    support = [sorted({key[0] for key in w}) for _,w in rows]

    @cache
    def completable(k, used):
        return k == n or any(not used & (1<<j) and completable(k+1, used | (1<<j))
                             for j in support[k])

    @cache
    def moves(k, used):
        # Inversion sign for appending column j to a partial permutation.
        return {j: (used | (1<<j), (-1)**((used >> (j+1)).bit_count()))
                for j in support[k] if not used & (1<<j) and completable(k+1, used | (1<<j))}

    den = 1
    states = {(0,0,0,0): (1,0,0,0,0,0,0,0)}
    for i, (row_den, weights) in enumerate(rows):
        nxt = defaultdict(lambda: [0]*8)
        for masks, jet in states.items():
            p1,p2,q1,q2 = (moves(i,s) for s in masks)
            for j,(mj,sj) in p1.items():
                for k,(mk,sk) in p2.items():
                    options = ((j,k),(k,j)) if j != k else ((v,v) for v in q1)
                    for u,v in options:
                        if u in q1 and v in q2:
                            mu,su = q1[u]
                            mv,sv = q2[v]
                            add_product(nxt, (mj,mk,mu,mv), jet, weights[j,k,u,v], sj*sk*su*sv)
        states = nxt
        den *= row_den
    jet = states[(full,)*4]
    return [Q2(jet[i], jet[i+1])/den for i in range(0,8,2)]


def verify_fixture(b):
    base, at_zero = fixture(b, 1), fixture(b, 0)
    slope = [[u-v for u,v in zip(a,m)] for a,m in zip(base,at_zero)]
    p = Q2(F((25,144)[b],169))
    for z in (-2,-1,0,1,2):
        assert permanent(fixture(b,z)) == ONE+(z-1)*p
    assert any(x != ZERO for row in slope for x in row)
    for i,j in combinations(range(11),2):
        for k,l in combinations(range(11),2):
            assert slope[i][k]*slope[j][l] == slope[i][l]*slope[j][k]
    # Changing the row order in all four determinants preserves their product.
    jet = second_moment(base, slope, (1,0,9,10,3,4,2,7,8,5,6))
    assert jet[1] == jet[2]
    var = jet[3]-p*p
    exact = (Q2(F(1826730516171,51181312), F(-1276512114551,51181312)),
             Q2(F(1706119275154625,51181312), F(-46399821321775,1968512)))
    assert var == exact[b]
    verify_marked_cycle(base, slope)
    print(f"b={b}, p={p.a}, variance={var.a} + ({var.b}) sqrt(2)", flush=True)
    print(f"variance approximately {float(var):.9f}", flush=True)
    return var


def verify_marked_cycle(base, slope):
    # Check one nonzero four-edge Walsh coefficient with a nonzero cycle derivative.
    # Multiply by sqrt(product of masses on the four edges) to remove the radicals.
    n = len(base)
    mass = [[abs(a) if (abs(a)-abs(d)).sign() >= 0 else abs(d)
             for a,d in zip(row,dr)] for row,dr in zip(base,slope)]
    for (i,j),(k,l) in product(combinations(range(n),2), repeat=2):
        if any(mass[r][c] == ZERO for r,c in product((i,j),(k,l))):
            continue
        def cycle_term(u,v):
            w = mass[i][u]*mass[j][v]
            return (w*base[i][v]*base[j][u],
                    w*(base[i][v]*slope[j][u]+slope[i][v]*base[j][u]))
        t0,t1 = cycle_term(k,l)
        v0,v1 = cycle_term(l,k)
        if t1+v1 == ZERO:
            continue
        rest_rows = [r for r in range(n) if r not in (i,j)]
        rest_cols = [c for c in range(n) if c not in (k,l)]
        @cache
        def minor(h,mask):
            if h == len(rest_rows):
                return ONE,ZERO
            z,dz = ZERO,ZERO
            r = rest_rows[h]
            for c in rest_cols:
                if mask & (1<<c) or mass[r][c] == ZERO:
                    continue
                v,d = minor(h+1,mask|(1<<c))
                z += base[r][c]*v
                dz += slope[r][c]*v+base[r][c]*d
            return z,dz
        z,dz = minor(0,0)
        expected = -(t1+v1)*z-(t0+v0)*dz
        if expected == ZERO or z == ZERO:
            continue
        # Direct two-determinant sign expectation at this mask. Track the two
        # partial permutations and their inversion signs independently.
        @cache
        def paired(r,m1,m2):
            if r == n:
                return ONE,ZERO
            exceptional = r in (i,j)
            opts = ((k,l),(l,k)) if exceptional else ((c,c) for c in range(n))
            total, deriv = ZERO,ZERO
            for u,v in opts:
                if m1&(1<<u) or m2&(1<<v) or mass[r][u]==ZERO or mass[r][v]==ZERO:
                    continue
                sign = (-1)**((m1>>(u+1)).bit_count()+(m2>>(v+1)).bit_count())
                w = sign*mass[r][u] if exceptional else sign
                tail, dtail = paired(r+1,m1|(1<<u),m2|(1<<v))
                total += w*base[r][v]*tail
                deriv += w*(base[r][v]*dtail+slope[r][v]*tail)
            return total,deriv
        assert paired(0,0,0)[1] == expected
        print(f"marked cycle rows {i,j}, columns {k,l}: exact coefficient verified")
        return
    raise AssertionError("no marked test cycle")


def check_direct_signs():
    base = [[ONE,ONE],[-ONE,ONE]]
    slope = [[ONE,ZERO],[ONE,ZERO]]
    total = [0]*4
    for h,j,k,l in product((-1,1),repeat=4):
        dp = h*l-j*k
        dq = h*l+j*k
        dv = h*l-j*k
        x,y = dp*dq,dp*dv
        total = [s+t for s,t in zip(total,(x*x,x*y,y*x,y*y))]
    assert second_moment(base,slope,(0,1)) == [Q2(F(s,16)) for s in total]


if __name__ == "__main__":
    check_direct_signs()
    for b in (0,1):
        verify_fixture(b)
