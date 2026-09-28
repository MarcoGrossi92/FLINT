#!/usr/bin/env python3
"""Generate test/andersen/WD-Andersen: the INPUT folder of the WD-Andersen mechanism written by a table writer
(argument 1: that folder, tables from 1 K to 3200 K on a 1 K step;
argument 2: the yaml, WD-andersen.yaml with reaction 3 written as the explicit inverse of reaction 2:
orders {CO2: 1.0, H2O: 0.5, O2: -0.25}, negative-orders and nonreactant-orders true) restricted to
T0..T1 K on one grid (thermo.dat, chemistry-Arrhenius.dat; phase.txt, composition.txt and
chemistry-info.txt copied), and reference.txt: at every state (T integer = a table node, random
composition and pressure; plus states with a zero concentration of O2, H2O, CO2 and CH4, and one with
[O2] of order 1e-14 kmol/m3, where the power is evaluated: no threshold) the partial
densities [kg/m3] and the net production rates [kg/m3/s] of Cantera. The script asserts that the
table's kf equals Cantera's at every reference temperature (relative 5e-13: the table carries 13 significant digits) and prints kf2/kf3 against
the equilibrium constant of CO + 0.5 O2 <-> CO2 (the closure the Andersen orders provide).
Run from the repository root:
    python3 test/andersen/make_WD-Andersen.py <INPUT folder of the table writer> <WD-Andersen.yaml>
"""
import os, re, sys, numpy as np, cantera as ct
src, yaml = sys.argv[1], sys.argv[2]
root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
dst = os.path.join(root, 'test', 'andersen', 'WD-Andersen'); os.makedirs(dst, exist_ok=True)
T0, T1 = 1100, 1400
def cut(name):
    txt = open(os.path.join(src, name)).read()
    idx = [m.start() for m in re.finditer(r'(?m)^ZONE T=', txt)]
    out = [txt[:idx[0]]]
    for k in range(len(idx)):
        z = txt[idx[k]:(idx[k+1] if k+1 < len(idx) else len(txt))].split('\n')
        assert z[1].startswith('I='), z[1]
        keep = [r for r in z[2:] if r.strip() and T0 <= float(r.split()[0]) <= T1]
        assert len(keep) == T1 - T0 + 1, (name, len(keep))
        out.append(z[0] + '\n' + 'I=%d, F=POINT\n' % len(keep) + '\n'.join(keep) + '\n')
    open(os.path.join(dst, name), 'w').write(''.join(out))
    return out
tab = cut('chemistry-Arrhenius.dat'); cut('thermo.dat')
for name in ('phase.txt', 'composition.txt', 'chemistry-info.txt'):
    open(os.path.join(dst, name), 'w').write(open(os.path.join(src, name)).read())
g = ct.Solution(yaml, 'WD-Andersen')
names = [l.split()[0] for l in open(os.path.join(dst, 'phase.txt')).read().split('\n')[1:] if l.strip()]
assert names == list(g.species_names) == ['CH4', 'O2', 'CO2', 'H2O', 'CO'], (names, g.species_names)
r3 = g.reactions()[2]; assert r3.orders == {'CO2': 1.0, 'H2O': 0.5, 'O2': -0.25}, r3.orders
rng = np.random.RandomState(20260925)
states = []
for k in range(9):
    X = rng.uniform(0.02, 1.0, 5); X /= X.sum()
    states.append((int(rng.randint(T0 + 1, T1)), float(rng.choice([101325.0, 5.0e5, 2.0e6])), X, 'random'))
def special(T, p, zero, tag):
    X = rng.uniform(0.05, 1.0, 5); X[names.index(zero)] = 0.0; X /= X.sum(); states.append((T, p, X, tag))
special(1200, 101325.0, 'O2', 'O2zero'); special(1350, 5.0e5, 'O2', 'O2zero'); special(1250, 101325.0, 'H2O', 'H2Ozero')
special(1300, 2.0e6, 'CO2', 'CO2zero'); special(1150, 101325.0, 'CH4', 'CH4zero')
X = rng.uniform(0.05, 1.0, 5); X[names.index('O2')] = 1e-12; X /= X.sum(); states.append((1275, 101325.0, X, 'O2tiny'))
lines = ['# WD-Andersen (Andersen orders on reaction 3), Cantera %s; per state: "T tag", then roi(1:ns) [kg/m3], then the net production rates [kg/m3/s]' % ct.__version__,
         '%d %d' % (g.n_species, len(states))]
W = g.molecular_weights
for (T, p, X, tag) in states:
    g.TPX = float(T), p, X
    kf = g.forward_rate_constants
    for r in range(g.n_reactions):
        z = tab[r + 1].split('\n'); row = [x for x in z[2:] if x.strip() and float(x.split()[0]) == T][0].split()
        assert abs(float(row[1]) - kf[r]) <= 5e-13 * abs(kf[r]), (T, r + 1, row[1], kf[r])
    wd = W * g.net_production_rates
    assert np.isfinite(wd).all()   # an O2 = 0 state gives zero rates: every WD step needs O2
    lines.append('%d %s' % (T, tag)); lines.append(' '.join('%.16e' % v for v in g.density * g.Y)); lines.append(' '.join('%.16e' % v for v in wd))
open(os.path.join(dst, 'reference.txt'), 'w').write('\n'.join(lines) + '\n')
# closure check: the pair 2/3 at steady state gives [CO2]/([CO][O2]^0.5) = kf2/kf3, to be compared with Kc(CO + 0.5 O2 <-> CO2)
e = ct.Solution(thermo='ideal-gas', kinetics='gas', species=[g.species(s) for s in names],
                reactions=[ct.Reaction(equation='CO + 0.5 O2 <=> CO2', rate=ct.ArrheniusRate(1.0, 0.0, 0.0))])
for T in (T0, 1250, T1, 2000, 3000):
    g.TP = float(T), 101325.0; e.TP = float(T), 101325.0
    kf = g.forward_rate_constants; Kc = e.equilibrium_constants[0]
    print('T = %4d K: kf2/kf3 = %.4e, Kc(CO + 0.5 O2 <-> CO2) = %.4e, ratio %.4f' % (T, kf[1] / kf[2], Kc, kf[1] / kf[2] / Kc))
print('written', dst, 'states', len(states), 'sizes', {f: os.path.getsize(os.path.join(dst, f)) for f in os.listdir(dst)})
