# Bell(8) dihedral orbit counts

The (4,2,2) class has 210 placements in 22 dihedral orbits. Its orbit-size histogram is `{2:1, 4:4, 8:10, 16:7}`, with total `2+4*4+8*10+16*7=210`.

A separate construction chooses the four-block and pairs its complement, then checks the dihedral orbits and the orbit-stabilizer identity. All four classes in `scripts/bell8_orbits.py` pass the enumeration checks. Each expected size list is checked for its declared length and total.
