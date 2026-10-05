#!/usr/bin/env python3
"""Compile the finite Lean proofs and reject unapproved axiom dependencies."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]
VERSION='4.33.0'
SOURCE=ROOT/'formal/CountingBridge.lean'
ALLOWED_AXIOMS={'propext','Classical.choice','Quot.sound'}
THEOREMS=(
    'simple_accounting','simple_collar_accounting','distinct_collar_accounting',
    'eighty_percent_accounting','ninety_percent_distinct_accounting',
    'certificate_denominator_positive','bounded_certificate_nonnegative',
    'negative_certificate_numerator','bounded_certificate_negative_majorant',
    'count_le_certificate_sum','tenth_count_cap','threshold_count_cap',
    'finite_bounded_certificate_counting','cubic_model_cap','cubic_model_margin',
    'cubic_counting_conversion','cubic_model_trace_identity','rational_parameter_admissible',
    'mirror_denominator_identity','mirror_denominator_positive',
    'mirror_certificate_nonnegative','mirror_certificate_le_two',
    'mirror_certificate_negative_majorant','mirror_denominator_dominates',
    'mirror_threshold_count_cap','finite_mirror_certificate_counting',
    'resolvent_quadratic_identity','resolvent_quadratic_minorant',
)


def compiler_path(explicit=None):
    if explicit:
        return str(explicit)
    if os.environ.get('ZEROZETA_LEAN'):
        return os.environ['ZEROZETA_LEAN']
    # A direct installed compiler also works when the elan proxy environment
    # is unavailable. This lookup never installs or downloads a toolchain.
    home=Path(os.environ.get('ELAN_HOME',str(Path.home()/'.elan')))
    binary='lean.exe' if os.name=='nt' else 'lean'
    installed=home/'toolchains'/f'leanprover--lean4---v{VERSION}'/'bin'/binary
    if installed.is_file():
        return str(installed)
    found=shutil.which('lean')
    if not found:
        raise RuntimeError('Lean 4.33.0 is required; supply --lean-path or ZEROZETA_LEAN.')
    return found


def audit_output(returncode,stdout,stderr):
    if returncode!=0:
        raise RuntimeError(f'Lean compilation failed (exit {returncode}).\n{stdout}\n{stderr}')
    if re.search(r'\b(?:error|warning):',stdout+'\n'+stderr):
        raise RuntimeError('Lean emitted a diagnostic; no proof result is accepted.')
    if stderr.strip():
        raise RuntimeError('Unexpected Lean stderr; no proof result is accepted.')
    dependencies={}
    for line in stdout.splitlines():
        match=re.fullmatch(r"'ZeroZeta\.([A-Za-z0-9_]+)' depends on axioms: \[(.*)\]",line)
        empty=re.fullmatch(r"'ZeroZeta\.([A-Za-z0-9_]+)' does not depend on any axioms",line)
        if match:
            name=match.group(1)
            axioms=[a.strip() for a in match.group(2).split(',') if a.strip()]
        elif empty:
            name=empty.group(1)
            axioms=[]
        else:
            if line.strip():
                raise RuntimeError(f'Unexpected Lean output: {line}')
            continue
        if name in dependencies:
            raise RuntimeError(f'Duplicate axiom report: {name}')
        extra=set(axioms)-ALLOWED_AXIOMS
        if extra:
            raise RuntimeError(f'Unapproved axiom dependencies in {name}: {sorted(extra)}')
        dependencies[name]=axioms
    if set(dependencies)!=set(THEOREMS):
        raise RuntimeError('The complete expected theorem axiom audit was not returned.')
    return dependencies


def rejection_self_test(stdout,binary,source):
    fixtures={
        'compile_error':(1,stdout,''),
        'missing_theorem_report':(0,'\n'.join(stdout.splitlines()[:-1])+'\n',''),
        'admitted_dependency':(0,stdout.replace('[propext, Quot.sound]','[propext, sorryAx]',1),''),
        'native_dependency':(0,stdout.replace('[propext, Quot.sound]','[propext, ZeroZeta._native.fake]',1),''),
        'custom_axiom_dependency':(0,stdout.replace('[propext, Quot.sound]','[propext, ZeroZeta.unproved_cap]',1),''),
        'warning':(0,stdout+'warning: incomplete proof\n',''),
    }
    rejected=[]
    for name,args in fixtures.items():
        try:
            audit_output(*args)
        except RuntimeError:
            rejected.append(name)
        else:
            raise RuntimeError(f'The proof gate accepted rejection fixture {name}.')
    old='4*target ≤ 5*simple + 10*collar + 9*deficit'
    if source.count(old)!=1:
        raise RuntimeError('The deliberate false-bound test must target exactly one theorem.')
    wrong=source.replace(old,'5*target ≤ 5*simple + 10*collar + 9*deficit',1)
    owned=Path(tempfile.mkdtemp(prefix='zerozeta-lean-audit-',dir=ROOT))
    if not owned.resolve().is_relative_to(ROOT.resolve()):
        raise RuntimeError('Unexpected temporary proof directory.')
    mutation=owned/'IncorrectBound.lean'
    try:
        mutation.write_text(wrong,encoding='utf-8')
        run=subprocess.run([binary,'-DwarningAsError=true',str(mutation)],cwd=ROOT,
                           capture_output=True,text=True,encoding='utf-8',timeout=120)
        if run.returncode==0:
            raise RuntimeError('Lean unexpectedly accepted the false stronger counting bound.')
        try:
            audit_output(run.returncode,run.stdout,run.stderr)
        except RuntimeError:
            rejected.append('false_stronger_counting_bound')
        else:
            raise RuntimeError('The proof gate accepted a failed false-bound compilation.')
    finally:
        if mutation.exists():
            mutation.unlink()
        owned.rmdir()
    return rejected,run.returncode


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lean-path',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    source=SOURCE.read_text(encoding='utf-8')
    declared=re.findall(r'^theorem\s+([A-Za-z0-9_]+)',source,re.M)
    if len(declared)!=len(THEOREMS) or set(declared)!=set(THEOREMS):
        raise RuntimeError('The source theorem list does not match the complete audit list.')
    if re.findall(r'^import\s+(.+)$',source,re.M)!=['Std']:
        raise RuntimeError('Only the pinned standard-library import is expected.')
    if re.search(r'\b(?:sorry|admit|axiom|native_decide)\b',source):
        raise RuntimeError('An admitted-proof, axiom declaration, or native-evaluation token was found.')
    toolchain=(ROOT/'lean-toolchain').read_text().strip()
    if toolchain!=f'leanprover/lean4:v{VERSION}':
        raise RuntimeError('Unexpected Lean toolchain pin.')
    binary=compiler_path(args.lean_path)
    version=subprocess.run([binary,'--version'],capture_output=True,text=True,
                           encoding='utf-8',timeout=30)
    if version.returncode or not re.search(rf'Lean \(version {re.escape(VERSION)},',version.stdout):
        raise RuntimeError(f'Lean {VERSION} is required; compiler version check failed.')
    run=subprocess.run([binary,'-DwarningAsError=true','formal/CountingBridge.lean'],
                       cwd=ROOT,capture_output=True,text=True,encoding='utf-8',timeout=120)
    dependencies=audit_output(run.returncode,run.stdout,run.stderr)
    rejected,false_exit=rejection_self_test(run.stdout,binary,source)
    record={
        'generated_at_utc':datetime.now(timezone.utc).isoformat(),
        'status':'Lean kernel accepted the scalar and finite conditional counting proofs; arithmetic cap unproved',
        'lean_version':version.stdout.strip(),
        'source':'formal/CountingBridge.lean',
        'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'toolchain_sha256':hashlib.sha256((ROOT/'lean-toolchain').read_bytes()).hexdigest(),
        'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'theorem_count':len(THEOREMS),
        'axiom_dependencies':dependencies,
        'allowed_standard_axioms':sorted(ALLOWED_AXIOMS),
        'lean_exit_code':run.returncode,
        'lean_stdout':run.stdout,
        'rejected_gate_cases':rejected,
        'deliberate_false_bound_lean_exit_code':false_exit,
        'limitations':[
            'Matrix inertia and the identification of the list with actual matrix eigenvalues are explicit outside inputs.',
            'The zero-block decomposition, tail estimate and asymptotic transfer are not formalized.',
            'The actual certificate trace cap is a hypothesis, not an established arithmetic theorem.',
            'The generic ordered-field proof does not construct the real root parameter or a Mathlib real-number adapter.',
            'Boundedness, Lipschitz regularity and resolvent compression are outside this formal file.',
            'Axiom auditing does not replace review of the intended theorem statement.',
        ],
    }
    if args.output:
        args.output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(f'PASS: {len(THEOREMS)} Lean kernel-checked theorems, complete standard-axiom audit, {len(rejected)} gate rejection cases.')
    print('The actual arithmetic trace cap remains unproved.')


if __name__=='__main__':
    main()
