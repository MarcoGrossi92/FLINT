#!/usr/bin/env python3
"""Generate test/chemistry/tables/WD-100K: database/WD restricted to 100..400 K, as a complete INPUT folder on
one temperature grid (thermo.dat, phase.txt, composition.txt, chemistry-info.txt, chemistry-Arrhenius.dat).

The rows are copied verbatim from the database/WD tables (value at T kelvin on row T), so the first
row of every table of the fixture is T = 100 K. test-tables loads it and checks that f_kf/f_kb return
the same rate as the full 1 K table at the same temperature; test-rhs-range uses it as a small INPUT
folder. Run from the repository root:
    python3 test/chemistry/tables/make_WD-100K.py
"""
import os, re
root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
src = os.path.join(root, 'database', 'WD'); dst = os.path.join(root, 'test', 'chemistry', 'tables', 'WD-100K')
T0, T1 = 100, 400

def zones(txt):
    idx = [m.start() for m in re.finditer(r'(?m)^ZONE T=', txt)]
    return txt[:idx[0]], [txt[idx[k]:(idx[k+1] if k+1 < len(idx) else len(txt))] for k in range(len(idx))]

def cut(z):
    L = z.split('\n'); rows = [r for r in L[2:] if r.strip()]
    keep = [r for r in rows if T0 <= float(r.split()[0]) <= T1]
    assert len(keep) == T1 - T0 + 1, len(keep)
    return L[0] + '\n' + 'I=%d, F=POINT\n' % len(keep) + '\n'.join(keep) + '\n'

os.makedirs(dst, exist_ok=True)
for name in ('chemistry-Arrhenius.dat', 'thermo.dat'):
    head, zs = zones(open(os.path.join(src, name)).read())
    open(os.path.join(dst, name), 'w').write(head + ''.join(cut(z) for z in zs))
for name in ('chemistry-info.txt', 'phase.txt', 'composition.txt'):
    open(os.path.join(dst, name), 'w').write(open(os.path.join(src, name)).read())
print('written', dst)
