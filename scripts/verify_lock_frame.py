#!/usr/bin/env python3
"""Independent exact audits of the resolved lock Fourier frame."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
import math


def add(a,b):
    return (a[0]+b[0],a[1]+b[1])


def mul(a,b):
    return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])


def scale(a,c):
    return (a[0]*c,a[1]*c)


def norm(a):
    return a[0]**2+a[1]**2


def gaussian_product(values):
    out=(Q(1),Q(0))
    for value in values:
        out=mul(out,value)
    return out


def gaussian_sum(values):
    out=(Q(0),Q(0))
    for value in values:
        out=add(out,value)
    return out


def character(k):
    return ((Q(1),Q(0)),(Q(0),Q(1)),(Q(-1),Q(0)),(Q(0),Q(-1)))[k%4]


def transform(a,k):
    assert len(a)==4
    return gaussian_sum(scale(character(-k*x),v) for x,v in enumerate(a))


def rectangle(f,g,s,t):
    n=len(f)
    return sum(f[x]*f[(x-s)%n]*g[(x+t)%n]*g[(x+t-s)%n] for x in range(n))


def full_locks(functions):
    n=len(functions[0])
    return {v:sum(functions[0][x]*math.prod(f[(x+h)%n] for f,h in zip(functions[1:],v))
                  for x in range(n))
            for v in product(range(n),repeat=len(functions)-1)}


def check_transforms():
    cases=0
    for kind in range(8):
        f=[Q((x+kind)%5-2,x%2+1) for x in range(4)]
        g=[Q((2*x+kind)%4,x%3+1) for x in range(4)]
        fh=[transform(f,a) for a in range(4)]
        gh=[transform(g,a) for a in range(4)]
        locks={(s,t):rectangle(f,g,s,t) for s,t in product(range(4),repeat=2)}
        hats={}
        for a,b in product(range(4),repeat=2):
            direct=gaussian_sum(scale(character(-a*s-b*t),value) for (s,t),value in locks.items())
            convoluted=scale(gaussian_sum(gaussian_product((fh[(a-b-c)%4],fh[(c-a)%4],
                                                             gh[(c+b)%4],gh[-c%4]))
                                              for c in range(4)),Q(1,4))
            assert direct==convoluted
            hats[a,b]=direct
            cases+=1
        assert sum(v*v for v in locks.values())==sum(norm(v) for v in hats.values())/16
        for s in range(4):
            paired_f=[f[x]*f[(x-s)%4] for x in range(4)]
            paired_g=[g[x]*g[(x-s)%4] for x in range(4)]
            for b in range(4):
                direct=gaussian_sum(scale(character(-b*t),locks[s,t]) for t in range(4))
                assert direct==mul(transform(paired_f,-b),transform(paired_g,b))
            assert sum(locks[s,t]**2 for t in range(4))==sum(
                norm(mul(transform(paired_f,-b),transform(paired_g,b))) for b in range(4))/4
    return cases


def check_full_locks():
    cases=0
    for slots,kind in product(range(2,5),range(3)):
        functions=[[Q((x+i+kind)%5-1,2) for x in range(4)] for i in range(slots)]
        models=[[v+Q((-1)**(x+i),3) for x,v in enumerate(f)] for i,f in enumerate(functions)]
        data=full_locks(functions)
        model_data=full_locks(models)
        spectra=[[transform(f,a) for a in range(4)] for f in functions]
        for frequencies in product(range(4),repeat=slots-1):
            direct=gaussian_sum(scale(character(-sum(a*h for a,h in zip(frequencies,v))),value)
                                for v,value in data.items())
            factorized=gaussian_product([spectra[0][-sum(frequencies)%4]]+
                                         [spectra[i][a] for i,a in enumerate(frequencies,1)])
            assert direct==factorized
            cases+=1
        cross=sum(data[v]*model_data[v] for v in data)
        direct=sum(math.prod(sum(f[x]*u[(x+s)%4] for x in range(4))
                              for f,u in zip(functions,models)) for s in range(4))
        assert cross==direct
        ds=[max(norm(transform([a-b for a,b in zip(f,u)],k)) for k in range(4))
            for f,u in zip(functions,models)]
        energies=[max(sum(v*v for v in f),sum(v*v for v in u)) for f,u in zip(functions,models)]
        bound=slots*sum(ds[i]*math.prod(energies[l] for l in range(slots) if l!=i) for i in range(slots))
        assert sum((data[v]-model_data[v])**2 for v in data)<=bound
        for i in range(slots):
            hybrid=functions[:i]+models[i:]
            replaced=functions[:i+1]+models[i+1:]
            left,right=full_locks(hybrid),full_locks(replaced)
            actual=sum((left[v]-right[v])**2 for v in left)
            other_energies=[sum(v*v for v in f) for l,f in enumerate(hybrid) if l!=i]
            assert actual<=ds[i]*math.prod(other_energies)
    return cases


def integer_corr(f):
    result=defaultdict(Q)
    for x,y in product(f,repeat=2):
        result[x-y]+=f[x]*f[y]
    return dict(result)


def integer_locks(f,g,b1=1,b2=1):
    result=defaultdict(Q)
    common=math.lcm(b1,b2)
    for m,m1,n,n1 in product(f,f,g,g):
        displacement=b2*(m-m1)
        if displacement==b1*(n-n1):
            assert displacement%common==0
            result[displacement//common,b1*n-b2*m]+=f[m]*f[m1]*g[n]*g[n1]
    return dict(result)


def check_counterexamples_and_dilations():
    a={0:Q(1),1:Q(1)}
    locks=integer_locks(a,a)
    ac=integer_corr(a)
    rho=defaultdict(Q)
    for s,t in product(ac,repeat=2):
        rho[s-t]+=ac[s]*ac[t]
    assert sum(v*v for v in locks.values())==8
    assert sum(v*v for v in rho.values())==70
    assert sum(locks.values())==rho[0]==6
    f=dict(enumerate(map(Q,(1,4,4,1,2))))
    u=dict(enumerate(map(Q,(2,5,2,2,1))))
    assert integer_corr(f)==integer_corr(u)
    lf,lu=integer_locks(f,a),integer_locks(u,a)
    assert sum(v*v for (k,j),v in lf.items() if k==1)==292
    assert sum(v*v for (k,j),v in lu.items() if k==1)==220
    cases=0
    for b1,b2 in product(range(1,5),repeat=2):
        locks=integer_locks(f,u,b1,b2)
        pf={b2*x:v for x,v in f.items()}
        pg={b1*x:v for x,v in u.items()}
        cf,cg=integer_corr(pf),integer_corr(pg)
        assert sum(locks.values())==sum(cf.get(s,Q(0))*cg.get(s,Q(0)) for s in set(cf)|set(cg))
        for k in {k for k,j in locks}:
            shift=math.lcm(b1,b2)*k
            assert sum(v for (kk,j),v in locks.items() if kk==k)==cf.get(shift,0)*cg.get(shift,0)
        cases+=1
    return cases


def check_cubes():
    cases=0
    for kind in range(8):
        h=[Q((x+kind)%4-1,3) for x in range(4)]
        v=[Q((2*x+kind)%5-2,4) for x in range(4)]
        e=[a-b for a,b in zip(h,v)]
        energy=sum(rectangle(h,h,s,t)**2 for s,t in product(range(4),repeat=2))
        def cube(function):
            return sum(math.prod(function[(x+aa*a+bb*b+cc*c)%4]
                                  for aa,bb,cc in product(range(2),repeat=3))
                       for x,a,b,c in product(range(4),repeat=4))
        assert energy==cube(h)
        error=sum((rectangle(h,h,s,t)-rectangle(v,v,s,t))**2 for s,t in product(range(4),repeat=2))
        bound_size=max(abs(z) for z in h+v)
        # Raise the U^3 stability inequality to the fourth power: all rational.
        assert (error/(16*4**4*bound_size**6))**4<=cube(e)/4**4
        cases+=1
    return cases


def cyclic_convolution(a,b):
    n=len(a)
    return [sum(a[j]*b[(i-j)%n] for j in range(n)) for i in range(n)]


def reduce_prime_character(coefficients):
    # Phi_p(z)=1+...+z^(p-1): subtract the last coefficient.
    return tuple(v-coefficients[-1] for v in coefficients[:-1])


def check_quadratic_characters():
    gauss_cases,cancellation_cases=0,0
    for p in (3,5,7,11,13):
        for a in range(p):
            coeff=[0]*p
            for x in range(p):
                coeff[(x*x-a*x)%p]+=1
            absolute_square=cyclic_convolution(coeff,[coeff[-i%p] for i in range(p)])
            assert reduce_prime_character(absolute_square)==(p,)+(0,)*(p-2)
            gauss_cases+=1
        for x,s,t in product(range(p),repeat=3):
            assert (x*x-(x-s)**2-(x+t)**2+(x+t-s)**2)%p==(-2*s*t)%p
        if p<5:
            continue
        for s,t in product(range(p),repeat=2):
            if not s or not t or s==t or (s+t)%p==0:
                continue
            offsets=(0,-s,t,t-s)
            surviving=[]
            for signs in product((-1,0,1),repeat=4):
                quadratic=sum(signs)%p
                linear=2*sum(a*d for a,d in zip(signs,offsets))%p
                if quadratic==linear==0:
                    surviving.append(signs)
            assert set(surviving)=={(0,0,0,0),(1,-1,-1,1),(-1,1,1,-1)}
            cancellation_cases+=1
        # The double character average in the cosine-square identity equals 1/p.
        coefficients=[0]*p
        for s,t in product(range(p),repeat=2):
            coefficients[(4*s*t)%p]+=1
        assert reduce_prime_character(coefficients)==(p,)+(0,)*(p-2)
    return gauss_cases,cancellation_cases


def main():
    print('Exact resolved two-dimensional transform cases:',check_transforms())
    print('Full-lock transform/product/stability frequency cases:',check_full_locks())
    print('Nonperiodic counterexamples passed; dilation cases:',check_counterexamples_and_dilations())
    print('Cubic identity and raised rational stability cases:',check_cubes())
    print('Quadratic Gauss/cancellation cases:',check_quadratic_characters())
    print('Resolved frame audits passed; no prime-specific variance theorem is inferred.')


if __name__=='__main__':
    main()
