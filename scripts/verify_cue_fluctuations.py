#!/usr/bin/env python3
"""Exact Weyl and connected-partition checks of Gram fluctuations.

The all-order proof is in notes/cue_gram_fluctuations.md. No simulation or
NumPy is used. Weyl values use the column Gram, independently of cumulants.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations, product
import math


@lru_cache(None)
def partitions(size):
    if not size:
        return ((),)
    output = []
    def visit(labels):
        if len(labels)==size:
            output.append(tuple(tuple(i for i,a in enumerate(labels) if a==j)
                                for j in range(max(labels)+1)))
            return
        for a in range(max(labels)+2):
            visit(labels+(a,))
    visit((0,))
    return tuple(output)


@lru_cache(None)
def cyclic_groups(size):
    # Set partitions followed by anchored cyclic orders, rather than cuts.
    return tuple((blocks[0],)+tail for blocks in partitions(size)
                 for tail in permutations(blocks[1:]))


def component_count(pi,sizes):
    owners = tuple(j for j,b in enumerate(sizes) for _ in range(b))
    parent = list(range(len(sizes)))
    def root(a):
        while parent[a]!=a:
            a = parent[a]
        return a
    for block in pi:
        for i in block[1:]:
            parent[root(owners[i])] = root(owners[block[0]])
    return len({root(a) for a in range(len(sizes))})


@lru_cache(None)
def trace_cumulant(frequencies,n):
    if sum(frequencies):
        return 0
    value = 0
    for groups in cyclic_groups(len(frequencies)):
        position,positions = 0,[0]
        for group in groups[:-1]:
            position += sum(frequencies[i] for i in group)
            positions.append(position)
        value += (-1)**(len(groups)-1)*max(n-max(positions)+min(positions),0)
    return value


def connected_cumulant(sizes,n):
    length = sum(sizes)
    connected = [pi for pi in partitions(length) if component_count(pi,sizes)==1]
    successors = []
    start = 0
    for size in sizes:
        successors.extend(start+(i+1)%size for i in range(size))
        start += size
    numerator = 0
    for x in product(range(n),repeat=length):
        frequencies = tuple(x[j]-x[i] for i,j in enumerate(successors))
        for pi in connected:
            value = 1
            for block in pi:
                # Sorting makes the cache symmetric; the trace cumulant is
                # symmetric in its slots, but their labels stay distinct.
                value *= trace_cumulant(tuple(sorted(frequencies[i] for i in block)),n)
                if not value:
                    break
            numerator += value
    return F(numerator,n**(length+len(sizes)))


def poly_multiply(a,b):
    result = defaultdict(int)
    for e,c in a.items():
        for f,d in b.items():
            result[tuple(x+y for x,y in zip(e,f))] += c*d
    return {e:c for e,c in result.items() if c}


def permutation_sign(p):
    return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))


class WeylJoint:
    def __init__(self,n,max_order=3,rows=None):
        self.n = n
        self.rows = n if rows is None else rows
        if not 1<=self.rows<=n:
            raise ValueError('row count must lie between 1 and n')
        zero = (0,)*(n-1)
        density = defaultdict(int)
        for p,q in product(permutations(range(n)),repeat=2):
            density[tuple(p[j]-q[j] for j in range(n-1))] += permutation_sign(p)*permutation_sign(q)
        self.density = dict(density)
        assert self.density[zero]==math.factorial(n)
        # Column Gram n V*V: sum_a (z_k/z_j)^a. Its nonzero
        # eigenvalues equal those of the row Gram, including rectangular V.
        entries = [[defaultdict(int) for _ in range(n)] for _ in range(n)]
        for j,k,a in product(range(n),range(n),range(self.rows)):
            exponent = tuple(a*((i==k)-(i==j)) for i in range(n-1))
            entries[j][k][exponent] += 1
        power = [[{zero:int(j==k)} for k in range(n)] for j in range(n)]
        self.traces = {0:{zero:n}}
        for order in range(1,max_order+1):
            new = [[defaultdict(int) for _ in range(n)] for _ in range(n)]
            for j,k,i in product(range(n),repeat=3):
                for e,c in poly_multiply(power[j][i],entries[i][k]).items():
                    new[j][k][e] += c
            power = [[dict(p) for p in row] for row in new]
            trace = defaultdict(int)
            for j in range(n):
                for e,c in power[j][j].items():
                    trace[e] += c
            self.traces[order] = dict(trace)

    @lru_cache(None)
    def moment(self,sizes):
        polynomial = {(0,)*(self.n-1):1}
        for size in sizes:
            polynomial = poly_multiply(polynomial,self.traces[size])
        value = sum(c*self.density.get(tuple(-x for x in e),0) for e,c in polynomial.items())
        return F(value,math.factorial(self.n)*self.n**sum(sizes)*self.rows**len(sizes))

    def cumulant(self,sizes):
        value = F(0)
        for rho in partitions(len(sizes)):
            term = (-1)**(len(rho)-1)*math.factorial(len(rho)-1)
            for group in rho:
                term *= self.moment(tuple(sorted(sizes[i] for i in group)))
            value += term
        return value

    def log_vandermonde_expectation(self):
        # log |Vandermonde(z)|^2 has Fourier coefficient -1/|k| on
        # each nonzero pair frequency. The Weyl density has finite support.
        value = F(0)
        for i in range(self.n):
            for j in range(i+1,self.n):
                for k in range(1,self.n):
                    exponent = tuple(k*((a==i)-(a==j)) for a in range(self.n-1))
                    reverse = tuple(-a for a in exponent)
                    value -= F(self.density.get(exponent,0)+self.density.get(reverse,0),
                               k*math.factorial(self.n))
        return value


def check_incidence(pi,cycles,inner_groups):
    length = sum(cycles)
    owners,inner_arcs = {},[]
    nodes = 0
    for block,groups in zip(pi,inner_groups):
        ids = list(range(nodes,nodes+len(groups)))
        nodes += len(groups)
        for node,group in zip(ids,groups):
            for i in group:
                owners[block[i]] = node
        inner_arcs.extend((ids[j],ids[(j+1)%len(ids)]) for j in range(len(ids)))
    arcs = list(inner_arcs)
    start = 0
    for size in cycles:
        arcs.extend((owners[start+i],owners[start+(i+1)%size]) for i in range(size))
        start += size
    assert len(arcs)==length+nodes
    divergence = [0]*nodes
    parent = list(range(nodes))
    def root(a):
        while parent[a]!=a:
            a = parent[a]
        return a
    for a,b in arcs:
        divergence[a] -= 1
        divergence[b] += 1
        parent[root(a)] = root(b)
    assert not any(divergence)  # A*1=0.
    components = len({root(a) for a in range(nodes)})
    assert components==component_count(pi,cycles)
    if components==1:
        assert length+nodes-(nodes-1)==length+1


def variance_second(n):
    return F(1,10*n)+F(1,6*n**3)-F(4,15*n**5)


CONNECTED_POLYNOMIALS = {
    (2,2): (F(0),F(-4,15),F(0),F(1,6),F(0),F(1,10)),
    (2,3): (F(0),F(0),F(-2,3),F(0),F(1,3),F(0),F(1,3)),
    (3,3): (F(0),F(6,35),F(0),F(-17,10),F(0),F(2,5),F(0),F(79,70)),
    (2,2,2): (F(0),F(-16,15),F(0),F(53,45),F(0),F(-7,45),F(0),F(2,45)),
}


def evaluate(polynomial,n):
    return sum(c*n**i for i,c in enumerate(polynomial))


def main():
    graph_cases = 0
    for sizes in ((2,2),(2,3),(2,2,2),(3,3)):
        for pi in partitions(sum(sizes)):
            t = component_count(pi,sizes)
            owners = tuple(j for j,b in enumerate(sizes) for _ in range(b))
            coefficient = 0
            for rho in partitions(len(sizes)):
                group_of = {i:a for a,group in enumerate(rho) for i in group}
                if all(len({group_of[owners[i]] for i in block})==1 for block in pi):
                    coefficient += (-1)**(len(rho)-1)*math.factorial(len(rho)-1)
            assert coefficient==int(t==1)
            for groups in product(*(cyclic_groups(len(block)) for block in pi)):
                check_incidence(pi,sizes,groups)
                graph_cases += 1
    # Prove the 26-term range table on the entire cone 0<=k<=l.
    vectors = ((1,0),(-1,0),(0,1),(0,-1))
    table = defaultdict(int)
    for groups in cyclic_groups(4):
        pos,positions = (0,0),[(0,0)]
        for group in groups[:-1]:
            pos = tuple(pos[j]+sum(vectors[i][j] for i in group) for j in range(2))
            positions.append(pos)
        upper = max(positions,key=lambda v:2*v[0]+3*v[1])
        lower = min(positions,key=lambda v:2*v[0]+3*v[1])
        for v in positions:
            # Nonnegative on both rays (1,1),(0,1) proves cone dominance.
            assert upper[1]-v[1]>=0 and sum(upper)-sum(v)>=0
            assert v[1]-lower[1]>=0 and sum(v)-sum(lower)>=0
        table[tuple(a-b for a,b in zip(upper,lower))] += (-1)**(len(groups)-1)
    table = {e:c for e,c in table.items() if c}
    assert table=={(-1,1):-1,(0,1):2,(1,1):-1},table
    assert len(cyclic_groups(4))==26
    cases = 0
    for sizes,p in CONNECTED_POLYNOMIALS.items():
        length = sum(sizes)
        assert len(p)==length+2 and p[0]==evaluate(p,1)==evaluate(p,-1)==0
        assert all(not c for i,c in enumerate(p) if i%2!=(length+1)%2)
        assert (math.factorial(length+1)*p[-1]).denominator==1
    for n in (1,2,3,4,5):
        weyl = WeylJoint(n)
        harmonic = sum(F(1,k) for k in range(1,n+1))
        assert weyl.log_vandermonde_expectation()==n*(harmonic-1)
        assert weyl.moment((2,))==F(4,3)-F(1,3*n*n)
        assert weyl.cumulant((2,2))==variance_second(n)
        for sizes in ((2,2),(2,3),(3,3),(2,2,2)):
            value = weyl.cumulant(sizes)
            assert value==evaluate(CONNECTED_POLYNOMIALS[sizes],n)/n**(sum(sizes)+len(sizes))
            if n<=3:
                direct = connected_cumulant(sizes,n)
                assert value==direct,(n,sizes,value,direct)
            print(f"n={n}, orders={sizes}: exact joint cumulant {value}")
            cases += 1
    for n in range(1,25):
        total = sum((n-k)*(n-l)*((k*k if k==l else 0)-max(k+l-n,0))
                    for k,l in product(range(1,n),repeat=2))
        assert F(4*total,n**6)==variance_second(n)
    assert [n**6*variance_second(n) for n in (2,3,4)]==[4,28,112]
    assert F(1,10)*F(79,70)-F(1,3)**2==F(11,6300)>0
    print(f"Cycle-component cancellation and incidence graph cases: {graph_cases}")
    print(f"Independent Weyl joint-cumulant checks: {cases}; variance double sums: 24 sizes")
    print("Exact variance: 1/(10n)+1/(6n^3)-4/(15n^5); limiting fluctuation variance 1/10.")
    print("Joint covariance for orders 2,3: [[1/10,1/3],[1/3,79/70]]; determinant 11/6300.")
    print("Independent logarithmic Weyl checks: 5 sizes; E log det(G_n)/n=H_n-1-log(n).")
    print("All-order concentration and Gaussian limits follow from the connected-network proof.")


if __name__=="__main__":
    main()
