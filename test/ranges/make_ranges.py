#!/usr/bin/env python3
"""Generate the fixtures of test-ranges (the temperature-grid contract of the tables) from
database/WD, database/CORIA (one falloff-Troe zone) and database/Gerlinger (transport zones):
  thermo-1400          phase/composition of WD, thermo.dat 1400..1600 K
  rates-equal          WD rate tables 1400..1600 K                      -> accepted
  rates-narrow         WD rate tables 1450..1600 K                      -> refused (ios = 6)
  rates-troe-equal     reaction 3 declared Troe, Troe table 1400..1600  -> accepted
  rates-troe-mismatch  reaction 3 declared Troe, Troe table 1450..1600  -> refused (ios = 6)
  transport-equal      transport.dat 1400..1600 K (5 zones)             -> accepted
  transport-shifted    transport.dat 1420..1600 K                       -> refused (ios = 3)
  rates-short          WD rate tables 1400..1550 K (last row differs)   -> refused (ios = 6)
  rates-troe-short     Troe table 1400..1550 K (last row differs)       -> refused (ios = 6)
  rates-lind-equal     reaction 3 declared Lindemann (database/Pelucchi zone 1), 1400..1600 -> accepted
  rates-lind-mismatch  the same Lindemann table 1450..1600 K            -> refused (ios = 6)
  diffusion-equal      diffusion.dat 1400..1600 K (10 pairs, constant D) -> accepted
  diffusion-shifted    diffusion.dat 1420..1600 K                       -> refused (ios = 3)
The species of the transport zones are renamed to the WD species: only the grid is under test.
Run from the repository root:  python3 test/ranges/make_ranges.py
"""
import os, re
root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
db = os.path.join(root, 'database'); dst = os.path.join(root, 'test', 'ranges')
wd = os.path.join(db, 'WD')

def zones(txt):
    idx = [m.start() for m in re.finditer(r'(?m)^ZONE T=', txt)]
    return txt[:idx[0]], [txt[idx[k]:(idx[k+1] if k+1 < len(idx) else len(txt))] for k in range(len(idx))]

def cut(z, T0, T1, title=None):
    L = z.split('\n'); rows = [r for r in L[2:] if r.strip()]
    keep = [r for r in rows if T0 <= float(r.split()[0]) <= T1]
    assert len(keep) == T1 - T0 + 1, len(keep)
    head = L[0] if title is None else 'ZONE T="%s"' % title
    return head + '\n' + 'I=%d, F=POINT\n' % len(keep) + '\n'.join(keep) + '\n'

def w(folder, name, txt):
    os.makedirs(os.path.join(dst, folder), exist_ok=True)
    open(os.path.join(dst, folder, name), 'w').write(txt)

thead, tz = zones(open(os.path.join(wd, 'thermo.dat')).read())
w('thermo-1400', 'thermo.dat', thead + ''.join(cut(z, 1400, 1600) for z in tz))
for f in ('phase.txt', 'composition.txt'):
    w('thermo-1400', f, open(os.path.join(wd, f)).read())
info = open(os.path.join(wd, 'chemistry-info.txt')).read()
ahead, az = zones(open(os.path.join(wd, 'chemistry-Arrhenius.dat')).read())
for folder, T0 in (('rates-equal', 1400), ('rates-narrow', 1450)):
    w(folder, 'chemistry-info.txt', info)
    w(folder, 'chemistry-Arrhenius.dat', ahead + ''.join(cut(z, T0, 1600) for z in az))
fhead, fz = zones(open(os.path.join(db, 'CORIA', 'chemistry-Troe.dat')).read())
for folder, T0 in (('rates-troe-equal', 1400), ('rates-troe-mismatch', 1450)):
    w(folder, 'chemistry-info.txt', info.replace('3 Arrhenius', '3 Troe'))
    w(folder, 'chemistry-Arrhenius.dat', ahead + ''.join(cut(z, 1400, 1600) for z in az))
    w(folder, 'chemistry-Troe.dat', fhead + cut(fz[0], T0, 1600, 'Reaction 1'))
ghead, gz = zones(open(os.path.join(db, 'Gerlinger', 'transport.dat')).read())
names = [l.split()[0] for l in open(os.path.join(wd, 'phase.txt')).read().split('\n')[1:] if l.strip()]
for folder, T0 in (('transport-equal', 1400), ('transport-shifted', 1420)):
    w(folder, 'transport.dat', ghead + ''.join(cut(gz[i], T0, 1600, names[i]) for i in range(len(names))))
# last rows: the grid check compares the first AND the last row of every table
w('rates-short', 'chemistry-info.txt', info)
w('rates-short', 'chemistry-Arrhenius.dat', ahead + ''.join(cut(z, 1400, 1550) for z in az))
w('rates-troe-short', 'chemistry-info.txt', info.replace('3 Arrhenius', '3 Troe'))
w('rates-troe-short', 'chemistry-Arrhenius.dat', ahead + ''.join(cut(z, 1400, 1600) for z in az))
w('rates-troe-short', 'chemistry-Troe.dat', fhead + cut(fz[0], 1400, 1550, 'Reaction 1'))
lhead, lz = zones(open(os.path.join(db, 'Pelucchi', 'chemistry-Lindemann.dat')).read())
for folder, T0 in (('rates-lind-equal', 1400), ('rates-lind-mismatch', 1450)):
    w(folder, 'chemistry-info.txt', info.replace('3 Arrhenius', '3 Lindemann'))
    w(folder, 'chemistry-Arrhenius.dat', ahead + ''.join(cut(z, 1400, 1600) for z in az))
    w(folder, 'chemistry-Lindemann.dat', lhead + cut(lz[0], T0, 1600, 'Reaction 1'))
# binary diffusion coefficients: synthetic constant D_ij (only the grid is under test), 10 pairs of WD species
pairs = [(a, b) for i, a in enumerate(names) for b in names[i+1:]]
for folder, T0 in (('diffusion-equal', 1400), ('diffusion-shifted', 1420)):
    txt = 'TITLE = "Binary diffusion coefficients (Pref=101325 Pa)"\nVARIABLES = "Temperature", "Dij"\n'
    for a, b in pairs:
        txt += 'ZONE T="%s-%s"\nI=%d, F=POINT\n' % (a, b, 1600 - T0 + 1) + ''.join('%.1f 1.0e-05\n' % T for T in range(T0, 1601))
    w(folder, 'diffusion.dat', txt)
print('written', dst)
