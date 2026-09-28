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
  transport-lower      transport.dat 1380..1600 K (starts below)        -> accepted, rows of transport-equal
  rates-short          WD rate tables 1400..1550 K (last row differs)   -> refused (ios = 6)
  rates-troe-short     Troe table 1400..1550 K (last row differs)       -> refused (ios = 6)
  rates-lind-equal     reaction 3 declared Lindemann (database/Pelucchi zone 1), 1400..1600 -> accepted
  rates-lind-mismatch  the same Lindemann table 1450..1600 K            -> refused (ios = 6)
  diffusion-equal      diffusion.dat 1400..1600 K (10 pairs, constant D) -> accepted
  diffusion-shifted    diffusion.dat 1420..1600 K                       -> refused (ios = 3)
  diffusion-lower      diffusion.dat 1380..1600 K (D = 5e-05 below 1400 K) -> accepted, rows of diffusion-equal
  rates-missing-zone   3 Arrhenius reactions, 2 zones in the table      -> refused (ios = 4)
  rates-step2          201 rows on a 2 K step (1400, 1402, .. 1800 K)   -> refused (ios = 6)
  thermo-step2         thermo.dat on a 2 K step (1400..1800 K)          -> refused (ios = 4)
  rates-troe-negk      Troe table with k_inf < 0 at 1500 K              -> refused (ios = 5)
  rates-lind-negk      Lindemann table with k_0 < 0 at 1500 K           -> refused (ios = 5)
  rates-zone2-shift    zone 1 on 1400..1600 K, zones 2.. on 1401..1601 K -> refused (ios = 6)
  rates-gap            zone 1 without row 1500 K and row 1499 K twice   -> refused (ios = 6)
  thermo-gap           the same gap in zone 1 of thermo.dat             -> refused (ios = 4)
  rates-troe-nan       Troe table with F_cent = NaN at 1500 K           -> refused (ios = 5)
  transport-gap        transport.dat 1400..1600 K, last zone without row 1500 K and row 1499 K twice -> refused (ios = 3)
  diffusion-gap        diffusion.dat 1400..1600 K, the same gap in the last pair -> refused (ios = 3)
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
for folder, T0 in (('transport-equal', 1400), ('transport-shifted', 1420), ('transport-lower', 1380)):
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
# binary diffusion coefficients: synthetic D_ij (only the grid is under test), 10 pairs of WD species;
# constant from 1400 K, another value below (rows that only diffusion-lower has)
pairs = [(a, b) for i, a in enumerate(names) for b in names[i+1:]]
for folder, T0 in (('diffusion-equal', 1400), ('diffusion-shifted', 1420), ('diffusion-lower', 1380)):
    txt = 'TITLE = "Binary diffusion coefficients (Pref=101325 Pa)"\nVARIABLES = "Temperature", "Dij"\n'
    for a, b in pairs:
        txt += 'ZONE T="%s-%s"\nI=%d, F=POINT\n' % (a, b, 1600 - T0 + 1) + ''.join('%.1f %s\n' % (T, '1.0e-05' if T >= 1400 else '5.0e-05') for T in range(T0, 1601))
    w(folder, 'diffusion.dat', txt)
# zone count and 1 K step
w('rates-missing-zone', 'chemistry-info.txt', info)
w('rates-missing-zone', 'chemistry-Arrhenius.dat', ahead + ''.join(cut(z, 1400, 1600) for z in az[:2]))
def cut2(z, T0, T1):
    L = z.split('\n'); rows = [r for r in L[2:] if r.strip() and T0 <= float(r.split()[0]) <= T1 and int(float(r.split()[0])) % 2 == 0]
    return L[0] + '\n' + 'I=%d, F=POINT\n' % len(rows) + '\n'.join(rows) + '\n'
w('rates-step2', 'chemistry-info.txt', info)
w('rates-step2', 'chemistry-Arrhenius.dat', ahead + ''.join(cut2(z, 1400, 1800) for z in az))
w('thermo-step2', 'thermo.dat', thead + ''.join(cut2(z, 1400, 1800) for z in tz))
for f in ('phase.txt', 'composition.txt'):
    w('thermo-step2', f, open(os.path.join(wd, f)).read())
# a negative limiting rate coefficient in one row of a falloff table (column 1 = k_inf, 2 = k_0)
def negate(z, T, col):
    L = z.split('\n')
    for i in range(2, len(L)):
        t = L[i].split()
        if t and float(t[0]) == T:
            t[col] = '-' + t[col].lstrip('-'); L[i] = ' '.join(t)
    return '\n'.join(L)
w('rates-troe-negk', 'chemistry-info.txt', info.replace('3 Arrhenius', '3 Troe'))
w('rates-troe-negk', 'chemistry-Arrhenius.dat', ahead + ''.join(cut(z, 1400, 1600) for z in az))
w('rates-troe-negk', 'chemistry-Troe.dat', fhead + negate(cut(fz[0], 1400, 1600, 'Reaction 1'), 1500.0, 1))
w('rates-lind-negk', 'chemistry-info.txt', info.replace('3 Arrhenius', '3 Lindemann'))
w('rates-lind-negk', 'chemistry-Arrhenius.dat', ahead + ''.join(cut(z, 1400, 1600) for z in az))
w('rates-lind-negk', 'chemistry-Lindemann.dat', lhead + negate(cut(lz[0], 1400, 1600, 'Reaction 1'), 1500.0, 2))
# every zone and every row of the 1 K grid, not only the first and the last row of zone 1: a later
# zone on another grid, and an interior row missing with another one twice (the first row, the last
# row and the row count of zone 1 are those of the 1 K grid)
w('rates-zone2-shift', 'chemistry-info.txt', info)
w('rates-zone2-shift', 'chemistry-Arrhenius.dat', ahead + cut(az[0], 1400, 1600) + ''.join(cut(z, 1401, 1601) for z in az[1:]))
def gap(z, T):   # row T replaced by a copy of row T - 1
    L = z.split('\n')
    i = next(k for k in range(2, len(L)) if L[k].split() and float(L[k].split()[0]) == T)
    L[i] = L[i-1]
    return '\n'.join(L)
w('rates-gap', 'chemistry-info.txt', info)
w('rates-gap', 'chemistry-Arrhenius.dat', ahead + gap(cut(az[0], 1400, 1600), 1500.0) + ''.join(cut(z, 1400, 1600) for z in az[1:]))
w('thermo-gap', 'thermo.dat', thead + gap(cut(tz[0], 1400, 1600), 1500.0) + ''.join(cut(z, 1400, 1600) for z in tz[1:]))
for f in ('phase.txt', 'composition.txt'):
    w('thermo-gap', f, open(os.path.join(wd, f)).read())
# a NaN F_cent (column 4) in one row of a falloff-Troe table
def setcol(z, T, col, val):
    L = z.split('\n')
    for i in range(2, len(L)):
        t = L[i].split()
        if t and float(t[0]) == T:
            t[col] = val; L[i] = ' '.join(t)
    return '\n'.join(L)
w('rates-troe-nan', 'chemistry-info.txt', info.replace('3 Arrhenius', '3 Troe'))
w('rates-troe-nan', 'chemistry-Arrhenius.dat', ahead + ''.join(cut(z, 1400, 1600) for z in az))
w('rates-troe-nan', 'chemistry-Troe.dat', fhead + setcol(cut(fz[0], 1400, 1600, 'Reaction 1'), 1500.0, 4, 'NaN'))
# the same gap in the LAST zone of a transport and of a binary-diffusion table (every zone, every row)
w('transport-gap', 'transport.dat', ghead + ''.join(cut(gz[i], 1400, 1600, names[i]) for i in range(len(names) - 1))
  + gap(cut(gz[len(names) - 1], 1400, 1600, names[-1]), 1500.0))
txt = 'TITLE = "Binary diffusion coefficients (Pref=101325 Pa)"\nVARIABLES = "Temperature", "Dij"\n'
for k, (a, b) in enumerate(pairs):
    Ts = [T if (k < len(pairs) - 1 or T != 1500) else 1499 for T in range(1400, 1601)]
    txt += 'ZONE T="%s-%s"\nI=%d, F=POINT\n' % (a, b, len(Ts)) + ''.join('%.1f 1.0e-05\n' % T for T in Ts)
w('diffusion-gap', 'diffusion.dat', txt)
print('written', dst)
