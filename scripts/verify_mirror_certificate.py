#!/usr/bin/env python3
"""Exact mirror-certificate identities, pole isolation and residue signs."""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
from itertools import permutations
import json
from math import gcd, isqrt, lcm
from pathlib import Path
from christoffel_exact import det, solve_linear
from model_moments import MODEL_MOMENTS

Q = [F(2519), F(-8232), F(7368), F(-1932)]
D = list(map(F, (6345361, 104885808, 86095872, 3732624)))
N = list(map(F, (41472816, 131040168, 28469952)))
P = [a*(-1)**j for j, a in enumerate(D)]
BRACKETS = [(F(1,16), F(1,15)), (F(1), F(2)), (F(20), F(24))]


def trim(p):
    p = list(p)
    while len(p)>1 and p[-1]==0:
        p.pop()
    return p


def primitive(p):
    """Scale by a positive rational to primitive integer coefficients."""
    p = trim(p)
    common = lcm(*(x.denominator for x in p))
    integers = [int(x*common) for x in p]
    content = gcd(*integers)
    return [F(x//content) for x in integers] if content else [F(0)]


def divide(a,b):
    """Exact Euclidean division with an independently checked identity."""
    a,b = trim(a),trim(b)
    assert b!=[0]
    r = list(a)
    q = [F(0)]*max(1,len(a)-len(b)+1)
    while r!=[0] and len(r)>=len(b):
        k = len(r)-len(b)
        c = r[-1]/b[-1]
        q[k] += c
        for j,x in enumerate(b):
            r[j+k] -= c*x
        r = trim(r)
    assert trim(add(mul(q,b),r))==a
    assert r==[0] or len(r)<len(b)
    return trim(q),r


def sturm(p):
    """Signed remainder chain, retaining sign under positive rescaling."""
    p = primitive(p)
    assert len(p)>1
    out = [p,primitive([F(j)*p[j] for j in range(1,len(p))])]
    while out[-1]!=[0] and len(out[-1])>1:
        _,r = divide(out[-2],out[-1])
        if r==[0]:
            break
        out.append(primitive(scale(r,F(-1))))
    signs = [[(-1 if q[-1]<0 else 1)*
              ((-1)**(len(q)-1) if side<0 else 1)
              for q in out] for side in (-1,1)]
    variations = [sum(a!=b for a,b in zip(s,s[1:])) for s in signs]
    return out,signs,variations


def sturm_self_checks():
    # Include repeated roots, mixed real/nonreal roots, sign reversal and
    # no-root polynomials; the expected counts follow from explicit factors.
    cases = 0
    for roots in ((),(F(0),),(F(-2),F(1)),(F(-2),F(-2),F(1)),
                  (F(-3),F(-1),F(0),F(2),F(2))):
        for nonreal in (1,2,3):
            p = [F(1)]
            for root in roots:
                p = mul(p,[-root,F(1)])
            for _ in range(nonreal):
                p = mul(p,[F(nonreal),F(0),F(1)])
            for sign in (F(-1),F(1)):
                _,_,variations = sturm(scale(p,sign))
                assert variations[0]-variations[1]==len(set(roots))
                cases += 1
    return cases


def add(a, b):
    return [sum((p[j] if j < len(p) else F(0)) for p in (a,b))
            for j in range(max(len(a),len(b)))]


def mul(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for j, x in enumerate(a):
        for k, y in enumerate(b):
            out[j+k] += x*y
    return out


def scale(a, c):
    return [c*x for x in a]


def evaluate(a, x):
    out = F(0)
    for c in reversed(a):
        out = out*x+c
    return out


def iadd(a,b):
    return a[0]+b[0], a[1]+b[1]


def imul(a,b):
    v = [x*y for x in a for y in b]
    return min(v), max(v)


def idiv(a,b):
    assert b[0]*b[1] > 0
    return imul(a,(1/b[1],1/b[0]))


def ievaluate(a,x):
    out = F(0),F(0)
    for c in reversed(a):
        out = iadd(imul(out,x),(c,c))
    return out


def sqrt_interval(v):
    scale_int = 10**18
    lo = isqrt(v[0].numerator*scale_int**2//v[0].denominator)
    hi = isqrt(v[1].numerator*scale_int**2//v[1].denominator)+1
    assert F(lo,scale_int)**2 <= v[0] <= v[1] <= F(hi,scale_int)**2
    return F(lo,scale_int),F(hi,scale_int)


def poles():
    records = []
    intervals = []
    for lo,hi in BRACKETS:
        original = lo,hi
        signs = evaluate(P,lo),evaluate(P,hi)
        assert signs[0]*signs[1] < 0
        for _ in range(80):
            mid = (lo+hi)/2
            if evaluate(P,lo)*evaluate(P,mid) < 0:
                hi = mid
            else:
                lo = mid
        v = lo,hi
        assert evaluate(P,lo)*evaluate(P,hi) < 0
        numerator = ievaluate([N[0],-N[1],N[2]],v)
        derivative = ievaluate([D[1],-2*D[2],3*D[3]],v)
        weight = idiv(numerator,derivative)
        assert weight[0] > 0
        r = sqrt_interval(v)
        intervals.append((r,weight))
        records.append(dict(initial_bracket=list(map(str,original)),
                            endpoint_values=list(map(str,signs)),
                            squared_pole_interval=list(map(str,v)),
                            pole_interval=list(map(str,r)),
                            weight_interval=list(map(str,weight))))
    # Three disjoint sign-changing intervals exhaust the degree-three P.
    contribution = F(0),F(0)
    for r,w in intervals[1:]:
        r2 = imul(r,r)
        lower = idiv(iadd((F(1),F(1)),idiv((F(-1,6),F(-1,6)),r)),
                     iadd((F(4,3),F(4,3)),r2))
        contribution = iadd(contribution,imul(w,lower))
    target = idiv(iadd((F(9,10),F(9,10)),(-contribution[1],-contribution[0])),
                  intervals[0][1])
    assert F(98282,100000) < target[0] <= target[1] < F(98283,100000)
    return records,contribution,target


def mmul(a,b):
    return [[sum(a[j][k]*b[k][i] for k in range(3)) for i in range(3)] for j in range(3)]


def minverse(a):
    m = [list(row)+[F(i==j) for j in range(3)] for i,row in enumerate(a)]
    for j in range(3):
        k = next(k for k in range(j,3) if m[k][j])
        m[j],m[k] = m[k],m[j]
        pivot = m[j][j]
        m[j] = [x/pivot for x in m[j]]
        for k in range(3):
            if k!=j:
                c = m[k][j]
                m[k] = [x-c*y for x,y in zip(m[k],m[j])]
    return [row[3:] for row in m]


def matrix_polynomial(a,m):
    out = [[F(0)]*3 for _ in range(3)]
    power = [[F(i==j) for j in range(3)] for i in range(3)]
    for c in a:
        out = [[out[i][j]+c*power[i][j] for j in range(3)] for i in range(3)]
        power = mmul(power,m)
    return out


def partial_fractions():
    # Companion-matrix identity proves the complete rational decomposition,
    # using polynomial coefficients rather than sampled evaluation points.
    m = [[F(0),F(0),-D[0]/D[3]], [F(1),F(0),-D[1]/D[3]],
         [F(0),F(1),-D[2]/D[3]]]
    dp = matrix_polynomial([D[1],2*D[2],3*D[3]],m)
    w = mmul(matrix_polynomial(N,m),minverse(dp))
    a = [[[-m[i][j],F(i==j)] for j in range(3)] for i in range(3)]
    determinant = [F(0)]
    for perm in permutations(range(3)):
        sign = (-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
        term = [F(sign)]
        for i in range(3):
            term = mul(term,a[i][perm[i]])
        determinant = add(determinant,term)
    assert determinant == scale(D,1/D[3])
    adj = [[None]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            rows = [k for k in range(3) if k!=j]
            cols = [k for k in range(3) if k!=i]
            minor = add(mul(a[rows[0]][cols[0]],a[rows[1]][cols[1]]),
                        scale(mul(a[rows[0]][cols[1]],a[rows[1]][cols[0]]),-1))
            adj[i][j] = scale(minor,(-1)**(i+j))
    numerator = [F(0)]
    for i in range(3):
        for j in range(3):
            numerator = add(numerator,scale(adj[j][i],w[i][j]))
    assert numerator == scale(N,1/D[3])
    assert sum(w[i][i] for i in range(3)) == N[2]/D[3]
    return dict(companion_determinant=list(map(str,determinant)),
                trace_adjugate_numerator=list(map(str,numerator)),
                sum_weights=str(N[2]/D[3]))


def scalar_identities():
    square = mul(Q,Q)
    assert square == list(map(F,(6345361,-41472816,104885808,-131040168,
                                86095872,-28469952,3732624)))
    assert square[::2] == D and scale(square[1::2],-1) == N
    f_one = evaluate(square,F(1))/evaluate(D,F(1))
    assert f_one == F(76729,201059665)
    assert f_one-F(1,10) == F(-8011695,80423866)
    # Clearing the positive 2*r^3*(x^2+r^2) denominator gives
    # 2*r^3*x - 2*r*x*(x^2+r^2) + x^2*(x^2+r^2)
    # = x^2*(x-r)^2. Compare every x-polynomial coefficient in Q[r].
    # Independently expand (x^2+r^2)*(x^2-2*r*x)+2*r^3*x.
    expanded = [[],[],[],[],[]]
    for i,ai in enumerate(([F(0),F(0),F(1)],[F(0)],[F(1)])):
        for j,bj in enumerate(([F(0)],[F(0),F(-2)],[F(1)])):
            expanded[i+j] = add(expanded[i+j],mul(ai,bj))
    expanded[1] = add(expanded[1],[F(0),F(0),F(0),F(2)])
    # Independently square x-r, then multiply by x^2.
    rhs = [[F(0)],[F(0)],[F(0)],[F(0)],[F(0)]]
    for i,ai in enumerate(([F(0),F(-1)],[F(1)])):
        for j,bj in enumerate(([F(0),F(-1)],[F(1)])):
            rhs[i+j+2] = add(rhs[i+j+2],mul(ai,bj))
    def trim(p):
        p = list(p)
        while len(p)>1 and p[-1]==0:
            p.pop()
        return p
    assert list(map(trim,expanded)) == list(map(trim,rhs))
    block_cap = sum(weight*evaluate(square,x)/evaluate(D,x*x)
                    for x,weight in ((F(0),F(1,6)),(F(1),F(2,3)),(F(2),F(1,6))))
    assert block_cap >= F(1,6) > F(1,10)
    return dict(square_coefficients=list(map(str,square)),
                value_at_one=str(f_one),prime_removal_threshold=str(f_one-F(1,10)),
                low_moment_block_certificate_value=str(block_cap))


def two_moment_minorants():
    cases = []
    for r in (F(1,4),F(1,3),F(1),F(4,3),F(2),F(4),F(5)):
        s = (r*r+F(4,3))/(r+1)
        p = r*r-r*s
        denominator = 2*r*s*s
        minorant_numerator = [-(r-s)**2,2*s,F(-1)]
        remainder = add([F(0),denominator],
                        scale(mul(minorant_numerator,[r*r,F(0),F(1)]),-1))
        assert remainder==mul([p,-s,F(1)],[p,-s,F(1)])
        expected = (1-F(1,6)/r)/(r*r+F(4,3))
        assert (minorant_numerator[0]+minorant_numerator[1]
                +minorant_numerator[2]*F(4,3))/denominator==expected
        value_at_mean = evaluate([p,-s,F(1)],F(1))
        assert value_at_mean==F(-1,3)
        assert p-s+F(4,3)==0
        coarse = 1/r**2-F(2,3)/r**3
        assert expected-coarse==(r-F(4,3))**2/(2*r**3*(r*r+F(4,3)))
        cases.append(dict(pole=str(r),quadratic_s=str(s),quadratic_p=str(p),
                          exact_mean_lower_bound=str(expected)))
    return cases


def rational_one_point_certificate():
    """Prove the global majorant directly with a rational Sturm certificate."""
    d = [F(0)]*7
    for j,c in enumerate(D):
        d[2*j] = c
    a = [F(100912,100000),F(-83810,100000),F(22949,100000)]
    gden = [F(121,1600),F(0),F(1)]
    weight = F(41501,100000)
    residual = add(add(mul(mul(a,d),gden),
                       scale(mul([F(0),F(1)],d),-weight)),
                   scale(mul(mul(Q,Q),gden),F(-1)))
    integer_residual = list(map(F,(700223277072,16130581267790,38453491895885,
        -657884432134880,686791679385776,656143236526080,4216355565563136,
        -7275624514820640,3177114169954896,-500529947904000,137055981081600)))
    assert scale(residual,F(160000000))==integer_residual
    chain,signs,variations = sturm(integer_residual)
    assert [len(p)-1 for p in chain]==list(range(10,-1,-1))
    assert signs==[[1,-1,-1,1,-1,-1,1,-1,-1,-1,-1],
                  [1,1,-1,-1,-1,1,1,1,-1,1,-1]]
    assert variations==[5,5] and integer_residual[0]>0
    # The nonzero constant final remainder proves squarefreeness; Sturm's
    # theorem and positivity at zero then prove positivity everywhere.
    assert len(chain[-1])==1 and chain[-1][0]!=0
    # Independent coefficient comparison for the rational-pole SOS.
    sos = add(add(scale(mul([F(21),F(-128),F(-96)],
                           [F(21),F(-128),F(-96)]),F(55)),
                  scale(mul([F(0),F(55),F(-384)],
                            [F(0),F(55),F(-384)]),F(32))),
              [F(0),F(0),F(0),F(0),F(983808)])
    assert sos==add(scale(mul([F(1),F(0),F(16)],
                             [F(1),F(0),F(16)]),F(24255)),
                    [F(0),F(-295680)])
    expectation = a[0]+a[1]+F(4,3)*a[2]
    threshold = (expectation-F(1,10))/weight
    cap = expectation-weight*F(91,100)
    conversion = F(1)-2*cap
    assert expectation==F(71551,150000)
    assert threshold==F(113102,124503)<F(91,100)
    assert cap==F(2980427,30000000)<F(1,10)
    assert conversion==F(12019573,15000000)
    inverse_t = F(1600,1721),F(440,1721)
    inverse_t_squared = cmul(inverse_t,inverse_t)
    assert inverse_t_squared==(F(2366400,2961841),F(1408000,2961841))
    closing_target = F(91,100)-inverse_t[0]
    assert closing_target==F(-3389,172100)
    assert inverse_t_squared[1]*F(11,40)==F(387200,2961841)
    return dict(majorant_polynomial=list(map(str,a)),
                rational_pole='11/40',resolvent_weight=str(weight),
                cleared_residual_multiplier='160000000',
                cleared_residual_coefficients=list(map(str,integer_residual)),
                sturm_chain_coefficients=[list(map(str,p)) for p in chain],
                sturm_signs_at_infinities=signs,sturm_variations=variations,
                positive_value_at_zero=str(integer_residual[0]),
                sturm_self_check_cases=sturm_self_checks(),
                cleared_derivative_sos_coefficients=list(map(str,sos)),
                moment_expectation=str(expectation),
                sufficient_actual_resolvent_threshold=str(threshold),
                sufficient_clean_actual_resolvent_threshold='91/100',
                conditional_trace_cap_at_clean_threshold=str(cap),
                conditional_simple_zero_proportion=str(conversion),
                prime_removal_inverse_t=list(map(str,inverse_t)),
                two_prime_inverse_t_squared=list(map(str,inverse_t_squared)),
                sufficient_prime_removal_real_target=str(closing_target))


def cadd(a,b):
    return a[0]+b[0],a[1]+b[1]


def cmul(a,b):
    return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]


def cscale(a,c):
    return a[0]*c,a[1]*c


def quadratic_inner(a,h,b):
    out = F(0),F(0)
    for i in range(len(a)):
        for j in range(len(b)):
            out = cadd(out,cscale(cmul((a[i][0],-a[i][1]),b[j]),h[i][j]))
    return out


def moment_disk(moments,r,n,shifted=False):
    """Exact Bessel disk for a supplied (possibly shifted) moment matrix."""
    offset = int(shifted)
    h = [[moments[i+j+offset] for j in range(n+1)] for i in range(n+1)]
    minors = [det([row[:k] for row in h[:k]]) for k in range(1,n+2)]
    assert all(x>0 for x in minors)
    pivots = [x/y for x,y in zip(minors,[F(1)]+minors[:-1])]
    columns = [solve_linear(h,[F(i==j) for i in range(n+1)]) for j in range(n+1)]
    inverse = list(map(list,zip(*columns)))
    assert all(sum(h[i][k]*inverse[k][j] for k in range(n+1))==F(i==j)
               for i in range(n+1) for j in range(n+1))
    z = F(0),r
    powers = [(F(1),F(0))]
    for _ in range(n+offset):
        powers.append(cmul(powers[-1],z))
    q = powers[offset:n+offset+1]
    s = []
    for j in range(n+1):
        value = F(0),F(0)
        for k in range(j+offset):
            value = cadd(value,cscale(powers[k],moments[j+offset-1-k]))
        s.append(value)
    a = quadratic_inner(q,inverse,q)
    b = quadratic_inner(q,inverse,s)
    c = quadratic_inner(s,inverse,s)
    assert a[1]==c[1]==0 and a[0]>0
    if shifted:
        xnum,ynum,radius_num = F(1,2)-b[0],-b[1],F(1,2)
    else:
        xnum,ynum,radius_num = -b[0],1/(2*r)-b[1],1/(2*r)
    assert xnum*xnum+ynum*ynum-a[0]*c[0]==radius_num**2
    center = xnum/a[0],ynum/a[0]
    radius = radius_num/a[0]
    return dict(pole=str(r),degree=n,shifted=shifted,
                leading_principal_minors=list(map(str,minors)),
                ldl_pivots=list(map(str,pivots)),
                inverse_moment_matrix=[list(map(str,row)) for row in inverse],
                bessel_a=str(a[0]),bessel_b=list(map(str,b)),bessel_c=str(c[0]),
                center=list(map(str,center)),radius=str(radius),
                real_lower=str(center[0]-radius),real_upper=str(center[0]+radius))


def model_resolvent_compatibility():
    disks = [moment_disk(MODEL_MOMENTS,F(11,40),n) for n in range(1,5)]
    lower = F(disks[-1]['real_lower'])
    assert list(map(F,disks[-1]['ldl_pivots']))==[
        F(1),F(1,3),F(5,36),F(247,5040),F(2448223,104569920)]
    assert lower==F(23741200782961777600,26054764615856959631)>F(91,100)
    assert F(disks[0]['real_lower'])==(1-F(1,6)/F(11,40))/(F(11,40)**2+F(4,3))
    obstruction = moment_disk(MODEL_MOMENTS,F(1,4),3,shifted=True)
    upper = F(obstruction['real_upper'])
    assert upper==F(742724016,753835003)<1
    return dict(scope='Moment disks for the continuum model; higher actual arithmetic moments are not assumed.',
                supplied_model_moments={str(j):str(x) for j,x in MODEL_MOMENTS.items()},
                candidate_pole_disks=disks,
                model_margin_above_clean_target=str(lower-F(91,100)),
                positive_support_quarter_pole_disk=obstruction,
                quarter_pole_gap_to_one=str(1-upper))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    isolated,contribution,target = poles()
    record = dict(schema_version=2,generated_at_utc=datetime.now(timezone.utc).isoformat(),
                  checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  scalar=scalar_identities(),partial_fractions=partial_fractions(),
                  two_moment_quadratic_minorants=two_moment_minorants(),
                  rational_one_point_certificate=rational_one_point_certificate(),
                  model_resolvent_compatibility=model_resolvent_compatibility(),
                  dependency_sha256={name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
                                     for name in ('model_moments.py','christoffel_exact.py')},
                  isolated_poles=isolated,
                  large_pole_moment_contribution_interval=list(map(str,contribution)),
                  sufficient_small_pole_target_interval=list(map(str,target)),
                  arithmetic_resolvent_cap_proved=False,
                  scope='Exact polynomial and interval certificates; no prime-phase or covariance estimate.')
    if args.output:
        args.output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print('PASS: mirror square identities; three simple imaginary pole pairs; positive weights;')
    print('partial fractions; quadratic lower bounds; exact rational one-point Sturm majorant.')
    print('Model-moment disks support the 11i/40 target and exclude a stronger i/4 target.')
    print('The actual arithmetic resolvent inequality remains unproved.')


if __name__=='__main__':
    main()
