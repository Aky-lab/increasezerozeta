#!/usr/bin/env python3
"""Finite prime-matrix experiments with a specified C2 taper and aliasing bound.

Requires NumPy. These are compressed prime matrices, not zero-sum Weil
matrices or asymptotic arithmetic estimates. The analytical quadrature
bound excludes floating-point roundoff.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from math import ceil, floor, isfinite, log, pi, sqrt
from pathlib import Path

import numpy as np


def prime_powers(cutoff):
    """Return (n,p,exponent,log(p)/sqrt(n)) for every prime power <=cutoff."""
    sieve = bytearray(b'\x01')*(cutoff+1)
    sieve[:2] = b'\x00\x00'
    for p in range(2,int(sqrt(cutoff))+1):
        if sieve[p]:
            sieve[p*p::p] = b'\x00'*((cutoff-p*p)//p+1)
    out = []
    for p in range(2,cutoff+1):
        if not sieve[p]:
            continue
        n,exponent = p,1
        while n<=cutoff:
            out.append((n,p,exponent,log(p)/sqrt(n)))
            n *= p
            exponent += 1
    return sorted(out)


def ramp(x):
    v = np.clip(x,0.0,1.0)
    return v**3*(10.0+v*(-15.0+6.0*v))


def taper(x,length):
    return ramp(length/2+x)*ramp(length/2-x)


def overlap_coefficients(length,shift,modes,nodes):
    if length<=2 or nodes<2*modes or nodes%2:
        raise ValueError('Require L>2 and an even midpoint grid N>=2*d.')
    u = -length/2+(np.arange(nodes)+0.5)*(length/nodes)
    samples = taper(u+shift/2,length)*taper(u-shift/2,length)
    frequencies = np.arange(modes)
    phases = (1-2*(frequencies%2))*np.exp(-1j*pi*frequencies/nodes)
    transformed = np.fft.rfft(samples)[:modes]*phases/nodes
    # The centered overlap is real and even. Report the discarded imaginary
    # component separately instead of treating it as an error certificate.
    return transformed.real,float(np.max(np.abs(transformed.imag)))


def matrix(length,ell1,height,dimension,nodes,powers):
    indices = np.arange(dimension)
    differences = np.abs(indices[:,None]-indices[None,:])
    frequency_sums = height+pi*(indices[:,None]+indices[None,:])/length
    h = np.zeros((dimension,dimension))
    max_imaginary = 0.0
    for n,p,exponent,weight in powers:
        s = log(n)
        coefficients,imaginary = overlap_coefficients(length,s,dimension,nodes)
        max_imaginary = max(max_imaginary,imaginary)
        h += (2*weight/ell1)*coefficients[differences]*np.cos(frequency_sums*s)
    baseline,imaginary = overlap_coefficients(length,0.0,dimension,nodes)
    max_imaginary = max(max_imaginary,imaginary)
    a0 = baseline[differences]
    return h,a0,max_imaginary


def mirror_value(x):
    q = 2519+x*(-8232+x*(7368-1932*x))
    y = x*x
    d = 6345361+y*(104885808+y*(86095872+3732624*y))
    return q*q/d


def statistics(eigenvalues,pole):
    z = np.mean(1/(eigenvalues-1j*pole))
    first = float(np.mean(eigenvalues))
    second = float(np.mean(eigenvalues**2))
    return dict(first_moment=first,second_moment=second,
                resolvent_real=float(z.real),resolvent_imaginary=float(z.imag),
                mirror_trace=float(np.mean(mirror_value(eigenvalues))),
                rational_majorant_trace_bound=1.00912-.83810*first+.22949*second-.41501*float(z.real),
                nonpositive_fraction=float(np.mean(eigenvalues<=0)),
                minimum_eigenvalue=float(np.min(eigenvalues)),
                maximum_eigenvalue=float(np.max(eigenvalues)))


def direct_entry(length,ell1,height,a,b,nodes,powers):
    """Integrate the original, uncentered F_ab formula without FFT."""
    u = -length/2+(np.arange(nodes)+0.5)*(length/nodes)
    total = 0.0
    for n,p,exponent,weight in powers:
        s = log(n)
        gates = taper(u,length)*taper(s-u,length)
        phase = height*s+(2*pi/length)*(a*u+b*(s-u))
        total += (2*weight/ell1)*np.mean(gates*np.cos(phase))
    return float(total)


def closed_kernel_calibration():
    maximum = 0.0
    cases = 0
    length = 6.0
    for shift in (1.0,2.0,4.0):
        coefficients,_ = overlap_coefficients(length,shift,129,32768)
        for k in (0,1,3,9,128):
            if k==0:
                exact = 1-(shift+1)/length
            else:
                w = 2*pi*k/length
                psi = 120*((12-w*w)*np.sin(w/2)-6*w*np.cos(w/2))/w**5
                exact = psi*np.sin(pi*k*(1-(shift+1)/length))/(pi*k)
            maximum = max(maximum,abs(coefficients[k]-exact))
            cases += 1
    if maximum>1e-10:
        raise ArithmeticError('FFT coefficient calibration failed')
    return dict(cases=cases,maximum_absolute_error=float(maximum))


def probe(height,nodes_floor,pole,max_dimension):
    if not isfinite(height) or height<=0:
        raise ValueError('Heights must be finite and positive.')
    ell = log(height/(2*pi))
    if ell<=2:
        raise ValueError('The quintic flat-top taper requires T>2*pi*exp(2).')
    ell1 = ell+2*log(2)-1
    dimension = floor(ell*height/(2*pi))
    if dimension>max_dimension:
        raise ValueError(f'd={dimension} exceeds --max-dimension {max_dimension}')
    cutoff = floor(np.exp(ell))
    powers = prime_powers(cutoff)
    nodes = 2**ceil(log(max(nodes_floor,2*dimension),2))
    h,a0,imaginary = matrix(ell,ell1,height,dimension,nodes,powers)
    h2,a02,imaginary2 = matrix(ell,ell1,height,dimension,2*nodes,powers)
    eigenvalues = np.linalg.eigvalsh(np.eye(dimension)-h)
    eigenvalues2 = np.linalg.eigvalsh(np.eye(dimension)-h2)
    frozen = np.linalg.eigvalsh(a0-h)
    selected = sorted({(0,0),(0,dimension-1),(dimension//3,dimension//2),
                       (dimension-1,dimension-1)})
    entry_checks = [dict(a=a,b=b,fft_entry=float(h[a,b]),
                        direct_entry=direct_entry(ell,ell1,height,a,b,nodes,powers))
                    for a,b in selected]
    direct_max = max(abs(c['fft_entry']-c['direct_entry']) for c in entry_checks)
    if direct_max>1e-9:
        raise ArithmeticError('Original-entry integral and FFT matrix disagree')
    # For |k|<=N/2, the C2 Fourier envelope and a monotone integral
    # bound give |c_grid-c|<=15*L/(2*N^2). See the proof note.
    coefficient_error = 15*ell/(2*nodes**2)
    weight_sum = sum(row[3] for row in powers)
    entry_error = 2*weight_sum*coefficient_error/ell1
    normalized_hs_error = sqrt(dimension)*entry_error
    frozen_hs_error = sqrt(dimension)*(entry_error+coefficient_error)
    coarse = statistics(eigenvalues,pole)
    fine = statistics(eigenvalues2,pole)
    return dict(height=height,ell=ell,ell1=ell1,bandwidth=1,
                prime_cutoff=cutoff,prime_power_count=len(powers),
                ordinary_prime_count=sum(e==1 for n,p,e,w in powers),
                prime_power_weight_sum=weight_sum,dimension=dimension,
                midpoint_nodes=nodes,comparison_nodes=2*nodes,
                matrix_symmetry_error=float(np.max(np.abs(h-h.T))),
                centered_coefficient_imaginary_residual=max(imaginary,imaginary2),
                original_entry_integral_checks=entry_checks,
                original_entry_integral_maximum_error=direct_max,
                analytic_coefficient_aliasing_bound=coefficient_error,
                analytic_entry_aliasing_bound=entry_error,
                analytic_normalized_hs_aliasing_bound=normalized_hs_error,
                analytic_resolvent_aliasing_bound=normalized_hs_error/pole**2,
                analytic_mirror_trace_aliasing_bound=7*normalized_hs_error,
                analytic_frozen_normalized_hs_aliasing_bound=frozen_hs_error,
                analytic_frozen_resolvent_aliasing_bound=frozen_hs_error/pole**2,
                analytic_frozen_mirror_trace_aliasing_bound=7*frozen_hs_error,
                grid_difference_normalized_hs=float(np.linalg.norm(h-h2,'fro')/sqrt(dimension)),
                grid_difference_frozen_normalized_hs=float(np.linalg.norm((a0-h)-(a02-h2),'fro')/sqrt(dimension)),
                grid_difference_resolvent_real=abs(coarse['resolvent_real']-fine['resolvent_real']),
                grid_difference_mirror_trace=abs(coarse['mirror_trace']-fine['mirror_trace']),
                compressed_prime_statistics=coarse,
                refined_grid_prime_statistics=fine,
                frozen_archimedean_statistics=statistics(frozen,pole),
                frozen_taper_trace=float(np.trace(a0)/dimension),
                exact_frozen_taper_trace_formula=1-281/(231*ell),
                frozen_taper_normalized_hs_deficit=float(np.linalg.norm(np.eye(dimension)-a0,'fro')/sqrt(dimension)),
                compressed_prime_eigenvalues=eigenvalues.tolist(),
                frozen_archimedean_eigenvalues=frozen.tolist())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--heights',type=float,nargs='+',default=[128,256,512,1024])
    parser.add_argument('--nodes',type=int,default=32768)
    parser.add_argument('--max-dimension',type=int,default=2048)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    if args.nodes<2 or args.max_dimension<1:
        parser.error('--nodes must be at least two and --max-dimension positive')
    cases = []
    for height in args.heights:
        row = probe(height,args.nodes,11/40,args.max_dimension)
        cases.append(row)
        stats = row['compressed_prime_statistics']
        print(f"T={height:g} d={row['dimension']} primes={row['ordinary_prime_count']} "
              f"Re m={stats['resolvent_real']:.9f} mirror={stats['mirror_trace']:.9f}",flush=True)
    record = dict(schema_version=1,generated_at_utc=datetime.now(timezone.utc).isoformat(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  numpy_version=np.__version__,taper='fixed-width C2 quintic flat top',
                  matrix='I-H_T using exact finite entry formula and midpoint Fourier quadrature',
                  pole='11/40',arithmetic_cap_proved=False,numeric_interval_certified=False,
                  closed_kernel_calibration=closed_kernel_calibration(),
                  aliasing_bound_excludes_roundoff=True,cases=cases,
                  limitations=[
                      'Finite heights cannot establish a limiting arithmetic estimate.',
                      'The computed matrix is a compressed prime matrix, not the zero-sum Weil matrix.',
                      'The quintic C2 taper is explicit; no equality with an unspecified smooth taper is assumed.',
                      'Analytical bounds control midpoint aliasing only, not FFT, phase or eigensolver roundoff.',
                      'Frozen archimedean statistics omit variation of the exact gamma-factor density and the pole term.',
                      'No zeta-zero bound follows from these numerical experiments.',
                  ])
    if args.output:
        args.output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
