"""Independent optical matrix check on a triple of live controlled signs."""

from itertools import product
from math import lcm

from optical import (
    A, F, INV_SQRT2, ONE, Q2, SQRT2, ZERO, beamsplitter, components, multiply,
)


C, T = F(77,85), F(36,85)
LEAF_C, LEAF_T = F(399,401), F(40,401)
Q = LEAF_T**2
LAMBDA = (1+A)/2


def rotation(fields, i, c, t):
    x,y = fields[i:i+2]
    fields[i:i+2] = [c*x-t*y,t*x+c*y]


def optical_matrix(n,z,cat=False):
    herald = 2*(n+1)
    gates = n + (n-1 if cat else 0)
    ports = [*range(0,herald,2), *range(herald,herald+4*gates)]
    columns = []
    for port in ports:
        f = [ZERO]*(herald+4*gates)
        f[port] = ONE
        gate = 0
        def cz(first,second):
            nonlocal gate
            beamsplitter(f,first,second)
            for arm,i in enumerate((first,second)):
                j = herald+2*gate+arm
                f[i],f[j] = A*f[i]+f[j],SQRT2*f[i]+f[j]
            beamsplitter(f,first,second)
            gate += 1
        def fanout(i):
            beamsplitter(f,2*i,2*i+1)
            cz(3,2*i+1)
            beamsplitter(f,2*i,2*i+1)
        c,t = (INV_SQRT2,INV_SQRT2) if cat else (C,T)
        rotation(f,0,c,t)
        if cat:
            rotation(f,2,F(5,13),F(12,13))
            for i in range(2,n+1):
                fanout(i)
        else:
            for i in range(1,n+1):
                rotation(f,2*i,LEAF_C,LEAF_T)
        for i in range(1,n+1):
            cz(1,2*i+1)
        rotation(f,0,c,t)
        f[1] *= z
        rotation(f,0,c,-t)
        for i in range(n,0,-1):
            cz(1,2*i+1)
        if cat:
            for i in range(n,1,-1):
                fanout(i)
            rotation(f,2,F(5,13),F(-12,13))
        else:
            for i in range(1,n+1):
                rotation(f,2*i,LEAF_C,-LEAF_T)
        rotation(f,0,c,-t)
        columns.append([f[i] for i in ports])
    return list(zip(*columns))


def permanent_jet(base, slope, rows, cols):
    # Subset recurrence with one derivative. Reorder sparse rows first.
    rows = sorted(rows,key=lambda i: sum(base[i][j]!=ZERO or slope[i][j]!=ZERO for j in cols))
    states = {0:(1,0,0,0)}
    denominator = 1
    for i in rows:
        terms = [(1<<k,base[i][j],slope[i][j]) for k,j in enumerate(cols)
                 if base[i][j]!=ZERO or slope[i][j]!=ZERO]
        den = lcm(*(v for _,a,d in terms for h in (a,d) for v in (h.a.denominator,h.b.denominator)))
        denominator *= den
        weights = [(bit,(int(a.a*den),int(a.b*den)),(int(d.a*den),int(d.b*den)))
                   for bit,a,d in terms]
        nxt = {}
        for mask,(x,y,dx,dy) in states.items():
            for bit,a,d in weights:
                if mask & bit:
                    continue
                u,v = multiply(a,(x,y))
                du,dv = multiply(a,(dx,dy))
                su,sv = multiply(d,(x,y))
                dest = mask|bit
                h,k,dh,dk = nxt.get(dest,(0,0,0,0))
                nxt[dest] = h+u,k+v,dh+du+su,dk+dv+sv
        states = nxt
    h,k,dh,dk = states.get((1<<len(cols))-1,(0,0,0,0))
    return Q2(h,k)/denominator,Q2(dh,dk)/denominator


def power(x,n):
    result = ONE
    for _ in range(n):
        result *= x
    return result


def check_triple(cat=False):
    n = 3
    base, at_zero = optical_matrix(n,1,cat), optical_matrix(n,0,cat)
    slope = [[a-b for a,b in zip(r,s)] for r,s in zip(base,at_zero)]
    assert len(list(components(base))) == 1
    ports = range(len(base))
    k = 2*C*C*T*T
    acceptance = Q2(F(25,169)) if cat else Q2(k*(1+(1-2*Q)**n))
    assert permanent_jet(base,slope,ports,ports) == (ONE,acceptance)
    derivatives = []
    for arms in product((0,1),repeat=n):
        first = n+1+(2*(n-1) if cat else 0)
        i,j,h = [first+2*m+arm for m,arm in enumerate(arms)]
        for r,c in ((i,0),(j,i),(h,j),(h,0),(i,i),(j,j)):
            assert slope[r][c]==ZERO
        direct = base[h][0]*base[i][i]*base[j][j]
        chain = base[i][0]*base[j][i]*base[h][j]
        u = INV_SQRT2 if cat else T
        assert direct == u*A*A and chain == u*INV_SQRT2*INV_SQRT2
        # Exercise the all-weightings proof at an asymmetric rational choice.
        # The orientation factors are reciprocal on the same six edges.
        def weight(row, col):
            return Q2(F(1+(row+2)*(col+3), 2+(row+col)%5))
        ratio = (weight(h,0)*weight(i,i)*weight(j,j)
                 / (weight(i,0)*weight(j,i)*weight(h,j)))
        x, y = ratio*chain, direct/ratio
        assert x*y == direct*chain and (x*y).sign() > 0
        assert (x+y)*(x+y)-4*direct*chain == (x-y)*(x-y)
        rows = [v for v in ports if v not in (i,j,h)]
        cols = [v for v in ports if v not in (0,i,j)]
        derivatives.append(permanent_jet(base,slope,rows,cols)[1])
    mean = sum(derivatives,ZERO)/len(derivatives)
    plus = power(1-Q+Q*LAMBDA,n-1)*(1-Q+Q*A)
    minus = power(1-Q-Q*LAMBDA,n-1)*(1-Q-Q*A)
    expected = Q2(F(25,338)) if cat else INV_SQRT2*T*C*C*(plus+minus)
    assert mean == expected
    rate = 8*A*A if cat else 8*A*A*power(1-Q+Q*LAMBDA,6)
    assert (rate-F(131,100)).sign() > 0
    assert k*(1+(1-2*Q)**150) < F(1,3) and 2*k > F(1,2)
    print(f"cat={cat}, one connected component, eight cycle complements, mean slope {float(mean):.12f}",flush=True)
    print(f"exponential base per triple {float(rate):.12f}")


def check_three_triples():
    n = 9
    base, at_zero = optical_matrix(n,1,True), optical_matrix(n,0,True)
    assert len(base) == 9*n-3 and len(list(components(base))) == 1
    first = n+1+2*(n-1)
    seen = set()
    for arms in product((0,1),repeat=n):
        selected = [first+2*m+arm for m,arm in enumerate(arms)]
        start = 0
        used_rows, used_columns, mask = set(), set(), set()
        for offset in range(0,n,3):
            i,j,h = selected[offset:offset+3]
            rows,cols = {i,j,h},{start,i,j}
            assert not used_rows & rows and not used_columns & cols
            used_rows |= rows
            used_columns |= cols
            edges = {(h,start),(i,i),(j,j),(i,start),(j,i),(h,j)}
            assert all(base[row][col]==at_zero[row][col] for row,col in edges)
            mask |= edges
            direct = base[h][start]*base[i][i]*base[j][j]
            chain = base[i][start]*base[j][i]*base[h][j]
            if offset == 0:
                u = INV_SQRT2
                assert direct == u*A*A and chain == u/2
            else:
                assert direct == A*A/SQRT2 and chain == ONE/(2*SQRT2)
            start = h
        seen.add(frozenset(mask))
    assert len(seen)==2**n
    print("three triples: 78 occupied ports, connected support, 512 distinct masks,"
          " all inter-triple endpoints and cycle products verified")


if __name__ == "__main__":
    check_triple()
    check_triple(cat=True)
    check_three_triples()
