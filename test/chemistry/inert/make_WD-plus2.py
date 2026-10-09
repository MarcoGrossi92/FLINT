#!/usr/bin/env python3
"""Generate test/chemistry/inert/WD-plus2: database/WD with two species (N2, AR) appended after the five
slots of the WD routine, as a table writer does for inert or mixing species that the compiled
routine does not know. Thermo and rate tables are restricted to 1400..1600 K to keep the fixture
small (the thermo of the two appended species is a copy of the CO zone: only their inertness is
tested). Run from the repository root:  python3 test/chemistry/inert/make_WD-plus2.py
"""
import os, re
root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
src = os.path.join(root, 'database', 'WD'); dst = os.path.join(root, 'test', 'chemistry', 'inert', 'WD-plus2')
T0, T1 = 1400, 1600
extra = [('N2', '28.014000', '0 0 0 2 0'), ('AR', '39.948000', '0 0 0 0 1')]
os.makedirs(dst, exist_ok=True)

def zones(txt):
    idx = [m.start() for m in re.finditer(r'(?m)^ZONE T=', txt)]
    return txt[:idx[0]], [txt[idx[k]:(idx[k+1] if k+1 < len(idx) else len(txt))] for k in range(len(idx))]

def cut(z, title=None):
    L = z.split('\n'); rows = [r for r in L[2:] if r.strip()]
    keep = [r for r in rows if T0 <= float(r.split()[0]) <= T1]
    assert len(keep) == T1 - T0 + 1, len(keep)
    head = L[0] if title is None else 'ZONE T="%s"' % title
    return head + '\n' + 'I=%d, F=POINT\n' % len(keep) + '\n'.join(keep) + '\n'

# thermo.dat: 5 zones cut + 2 copies of the CO zone
head, zs = zones(open(os.path.join(src, 'thermo.dat')).read())
out = [head] + [cut(z) for z in zs] + [cut(zs[4], e[0]) for e in extra]
open(os.path.join(dst, 'thermo.dat'), 'w').write(''.join(out))
# chemistry-Arrhenius.dat: 3 zones cut
head, zs = zones(open(os.path.join(src, 'chemistry-Arrhenius.dat')).read())
open(os.path.join(dst, 'chemistry-Arrhenius.dat'), 'w').write(head + ''.join(cut(z) for z in zs))
# phase.txt: two rows appended
open(os.path.join(dst, 'phase.txt'), 'w').write(open(os.path.join(src, 'phase.txt')).read().rstrip('\n') + '\n' +
                                                ''.join('%s %s\n' % (e[0], e[1]) for e in extra))
# composition.txt: elements N and Ar added, rows widened
open(os.path.join(dst, 'composition.txt'), 'w').write('5\nC\nH\nO\nN\nAr\n\ncomposition\n' +
    ''.join(r.strip() + ' 0 0\n' for r in open(os.path.join(src, 'composition.txt')).read().split('composition\n')[1].split('\n') if r.strip()) +
    ''.join(e[2] + '\n' for e in extra))
# chemistry-info.txt: species count and two zero rows per reaction before the M row
info = open(os.path.join(src, 'chemistry-info.txt')).read().replace('N.ro species = 5', 'N.ro species = 7')
lines = []
for L in info.split('\n'):
    m = re.match(r'^(\d+) M ', L)
    if m:
        for e in extra: lines.append('%s %s 0.0 0.0 0.0' % (m.group(1), e[0]))
    lines.append(L)
open(os.path.join(dst, 'chemistry-info.txt'), 'w').write('\n'.join(lines))
print('written', dst)
