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
    ("christoffel_exact.py", "exact moment consumption for supplied moments"),
    ("k8_target.py", "exact target geometry for supplied moments"),
    ("verify_spectator_reduction.py", "finite rational checks of exact-v reduction"),
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
        "all_checks_passed": all(c["passed"] for c in checks),
        "limitations": [
            "Passing checks do not prove the surrounding analytic transport.",
            "Higher-moment fractions and the 80% target remain candidates.",
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
        print(payload, end="")
    return 0 if report["all_checks_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
