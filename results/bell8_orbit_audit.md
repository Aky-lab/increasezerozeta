# Bell(8) expected-orbit audit — 5 October 2026

The historical (4,2,2) expected size list had 22 entries summing to 218 for a class with 210 placements. The enumerator returned 22 orbits and 210 total placements but failed the inconsistent expected-list gate.

The corrected histogram is `{2:1, 4:4, 8:10, 16:7}`. It has 22 orbits and total size `2+4*4+8*10+16*7=210`. A separate check directly generated the placements by choosing the four-block and pairing its complement, then verified all dihedral orbits and the orbit-stabilizer identity. This agreed with the corrected histogram. The 22-orbit count and representatives are unchanged.

All four classes in `scripts/bell8_orbits.py` pass the enumeration checks. A preflight check verifies each expected list has the declared length and total before comparing it with the enumerated orbits.
