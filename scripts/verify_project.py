#!/usr/bin/env python3
"""Run bounded local checks and record their scope, hashes and outputs."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
CHECKS = (
    ("verify_lemmas.py", "finite local algebra and conditional core arithmetic"),
    ("verify_arithmetic_transport.py", "exact coefficient convolution, prime-power factors and overlap geometry"),
    ("verify_prime_pair_transport.py", "divisor-family aggregation and published-input prime-pair error bookkeeping"),
    ("verify_short_window_transport.py", "integer-window discretization, Abel weights and coupled-product marginal concentration"),
    ("verify_presieved_pairs.py", "finite presieved Fourier coefficients, CRT pair main terms, Rankin Euler factors and Parseval"),
    ("christoffel_exact.py", "exact moment consumption for supplied moments"),
    ("k8_target.py", "exact target geometry for supplied moments"),
    ("verify_spectator_reduction.py", "finite rational checks of actual-walk exact-v reduction"),
    ("verify_model_certificates.py", "recorded integer certificates and source provenance"),
    ("verify_pairing_model.py", "pair-cycle certificates, direct counts and exact cell integration"),
    ("verify_mixed_cycles.py", "mixed-cycle certificates and independent scalar frequency sums"),
    ("verify_reciprocity.py", "reduced reciprocity certificates and all earlier lattice counts"),
    ("verify_cue_model.py", "independent exact Weyl-integration checks of Haar-unitary Gram moments"),
    ("verify_cue_fluctuations.py", "connected multi-cycle identities, exact Weyl joint cumulants and logarithmic determinant identity"),
    ("verify_cue_bandwidth.py", "rectangular Gram normalizations, Weyl/connected identities and exact bandwidth variance"),
    ("verify_sine_bridge.py", "finite determinantal Campbell expansion, exact Fourier gauge and sine second-moment normalization"),
    ("bell7_gate.py", "seventh-order partition multiplicities"),
    ("bell8_orbits.py", "eighth-order dihedral orbits and multiplicities"),
    ("m7_ledger.py", "seventh-order continuum model ledger"),
    ("m8_ledger.py", "eighth-order continuum model ledger and targets"),
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    checks = []
    for name, scope in CHECKS:
        path = ROOT / "scripts" / name
        record = {
            "script": "scripts/" + name,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "scope": scope,
        }
        if name == "verify_spectator_reduction.py":
            dependency = ROOT / "scripts" / "spectator_reduction.py"
            record["dependency_sha256"] = {
                "scripts/spectator_reduction.py": hashlib.sha256(dependency.read_bytes()).hexdigest()
            }
        try:
            completed = subprocess.run(
                [sys.executable, "-B", str(path)], cwd=ROOT,
                capture_output=True, text=True, timeout=120,
            )
            record.update(exit_code=completed.returncode,
                          stdout=completed.stdout, stderr=completed.stderr,
                          passed=completed.returncode == 0)
        except subprocess.TimeoutExpired:
            record.update(passed=False, error="120 second timeout")
        checks.append(record)
    report = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "script_sources_sha256": {
            path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted((ROOT/"scripts").glob("*.py"))
        },
        "all_checks_passed": all(c["passed"] for c in checks),
        "limitations": [
            "Passing checks do not prove the surrounding analytic transport.",
            "Finite model moments do not establish arithmetic transport to zeta zeros.",
            "The eighth-order independent NumPy checks are verified from their stored record.",
            "Spectator checks are finite examples, not an exhaustive proof.",
            "Hashes identify the executed scripts, not a full repository revision.",
        ],
        "checks": checks,
    }
    payload = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
        print(f"Evidence written to {args.output}")
        for check in checks:
            print(f"{'PASS' if check['passed'] else 'FAIL'}: {check['script']}")
            if not check["passed"]:
                print(check.get("error", ""))
                print(check.get("stdout", ""))
                print(check.get("stderr", ""))
    else:
        for check in checks:
            print(f"{'PASS' if check['passed'] else 'FAIL'}: {check['script']}")
            if not check["passed"]:
                print(check.get("error", ""))
                print(check.get("stdout", ""))
                print(check.get("stderr", ""))
        print(f"{sum(c['passed'] for c in checks)}/{len(checks)} checks passed")
    return 0 if report["all_checks_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
