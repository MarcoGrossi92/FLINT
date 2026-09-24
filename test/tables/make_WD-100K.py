#!/usr/bin/env python3
"""Generate test/tables/WD-100K: the WD chemistry tables of database/WD restricted to 100..400 K.

The rows are copied verbatim from database/WD/chemistry-Arrhenius.dat (rate at T kelvin on
row T), so the first row of the fixture is T = 100 K. test-tables loads it and checks that
f_kf/f_kb return the same rate as the
full 1 K table at the same temperature. Run from the repository root:
    python3 test/tables/make_WD-100K.py
"""
import os, re
root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src = os.path.join(root, 'database', 'WD'); dst = os.path.join(root, 'test', 'tables', 'WD-100K')
T0, T1 = 100, 400
txt = open(os.path.join(src, 'chemistry-Arrhenius.dat')).read()
idx = [m.start() for m in re.finditer(r'(?m)^ZONE T=', txt)]
head = txt[:idx[0]]
out = [head]
for k in range(len(idx)):
    z = txt[idx[k]:(idx[k+1] if k+1 < len(idx) else len(txt))].split('\n')
    rows = [r for r in z[2:] if r.strip()]
    keep = [r for r in rows if T0 <= float(r.split()[0]) <= T1]
    assert len(keep) == T1 - T0 + 1, len(keep)
    out.append(z[0] + '\n' + 'I=%d, F=POINT\n' % len(keep) + '\n'.join(keep) + '\n')
os.makedirs(dst, exist_ok=True)
open(os.path.join(dst, 'chemistry-Arrhenius.dat'), 'w').write(''.join(out))
open(os.path.join(dst, 'chemistry-info.txt'), 'w').write(open(os.path.join(src, 'chemistry-info.txt')).read())
print('written', dst)
