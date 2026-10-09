#!/usr/bin/env python3
"""Generate test/chemistry/orders/JLR-frassoldati: the INPUT folder of the JLR-frassoldati mechanism written by a table writer
(argument 1: that folder, thermo from the NASA9 records, transport from Cantera, tables up to 3000 K; argument 2:
the yaml) restricted to T0..T1 K, with the
'Reaction orders' block (explicit reaction orders) appended to chemistry-info.txt for reaction 1
(orders: {CH4: 0.5, O2: 1.3}), and reference.txt: at three states (T integer, p = 1 atm) the partial
densities and the FORWARD production rates [kg/m3/s] of Cantera with the yaml orders (what the block must
reproduce), of the integer-rounded stoichiometric law (what the general procedure computed without the
block before test-stoich) and of Cantera with the yaml orders removed (the mass-action law of the
stoichiometric coefficients: what the general procedure computes without the block). Forward only: the reverse rate constants of the table (kb = kf/Kc) come from the thermo
database chosen by the table writer (NASA9 here), not from the thermo of the yaml, so the reverse rates
of the two sides differ by construction; the driver zeroes kb_tab and compares the forward part, which
is where the reaction orders act. The script checks that the table's kf equals Cantera's at the three
temperatures (relative 1e-10) and prints the kb ratio. Run from the repository root:
    python3 test/chemistry/orders/make_JLR-frassoldati.py <INPUT folder of the table writer> <JLR-frassoldati.yaml>
"""
import os, re, sys, numpy as np, cantera as ct
src, yaml = sys.argv[1], sys.argv[2]
root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
dst = os.path.join(root, 'test', 'chemistry', 'orders', 'JLR-frassoldati'); os.makedirs(dst, exist_ok=True)
T0, T1 = 1200, 1800
txt = open(os.path.join(src, 'chemistry-Arrhenius.dat')).read()
idx = [m.start() for m in re.finditer(r'(?m)^ZONE T=', txt)]
out = [txt[:idx[0]]]
for k in range(len(idx)):
    z = txt[idx[k]:(idx[k+1] if k+1 < len(idx) else len(txt))].split('\n')
    keep = [r for r in z[2:] if r.strip() and T0 <= float(r.split()[0]) <= T1]
    assert len(keep) == T1 - T0 + 1
    out.append(z[0] + '\n' + 'I=%d, F=POINT\n' % len(keep) + '\n'.join(keep) + '\n')
open(os.path.join(dst, 'chemistry-Arrhenius.dat'), 'w').write(''.join(out))
# a table writer that writes the 'Reaction orders' block has already appended it: keep the part before it
info = open(os.path.join(src, 'chemistry-info.txt')).read().split('\nReaction orders\n')[0].rstrip('\n') + '\n'
open(os.path.join(dst, 'chemistry-info.txt'), 'w').write(info + '\nReaction orders\n2\n1 CH4 0.5\n1 O2 1.3\n')
open(os.path.join(dst, 'chemistry-info-noblock.txt'), 'w').write(info)
open(os.path.join(dst, 'phase.txt'), 'w').write(open(os.path.join(src, 'phase.txt')).read())
g = ct.Solution(yaml, 'JLR-Frassoldati')
names = [l.split()[0] for l in open(os.path.join(dst, 'phase.txt')).read().split('\n')[1:] if l.strip()]
assert names == list(g.species_names), (names, g.species_names)
X = {'O2': 0.25, 'CH4': 0.10, 'H2O': 0.15, 'CO': 0.08, 'CO2': 0.10, 'H2': 0.12, 'H': 0.02, 'O': 0.03, 'OH': 0.05}
nu1, nu2 = g.reactant_stoich_coeffs, g.product_stoich_coeffs
lines = ['# JLR-frassoldati, p = 101325 Pa, Cantera %s; per state: T, then roi(1:ns) [kg/m3], then forward wdot with the yaml orders, then forward wdot with nint(stoichiometric) orders, then forward wdot without the yaml orders [kg/m3/s]' % ct.__version__, '%d %d' % (g.n_species, 3)]
rs = []   # independent copies of the reactions without their orders (setting .orders on g.reactions() edits g's objects)
for r in g.reactions():
    d = dict(r.input_data)
    for key in ('orders', 'negative-orders', 'nonreactant-orders'): d.pop(key, None)
    rs.append(ct.Reaction.from_dict(d, g))
g0 = ct.Solution(thermo='ideal-gas', kinetics='gas', species=g.species(), reactions=rs)   # the yaml without its orders
for T in (1300.0, 1500.0, 1750.0):
    g.TPX = T, 101325.0, X
    c, W = g.concentrations, g.molecular_weights
    kf, kb = g.forward_rate_constants, g.reverse_rate_constants
    wd_orders = W * ((nu2 - nu1) @ g.forward_rates_of_progress)     # yaml orders, forward only
    wd_nint = np.zeros(g.n_species)
    for r in range(g.n_reactions):
        e1 = np.floor(nu1[:, r] + 0.5)      # Fortran NINT (half away from zero, coefficients >= 0)
        wd_nint += W * (nu2[:, r] - nu1[:, r]) * kf[r] * np.prod(c ** e1)
    # the table (row T, T integer) must carry Cantera's kf; kb differs (thermo source), print the ratio
    for r in range(g.n_reactions):
        z = out[r + 1].split('\n'); row = [x for x in z[2:] if x.strip() and float(x.split()[0]) == T][0].split()
        kft, kbt = float(row[1]), float(row[2])
        assert abs(kft - kf[r]) <= 1e-10 * abs(kf[r]), (T, r + 1, kft, kf[r])
        print('T = %.0f K reaction %d: kf_tab/kf_cantera - 1 = %.1e, kb_tab/kb_cantera = %s' % (T, r + 1, kft / kf[r] - 1, ('%.4f' % (kbt / kb[r])) if kb[r] > 0 else 'irreversible'))
    g0.TPX = T, 101325.0, X
    wd_stoich = W * ((nu2 - nu1) @ g0.forward_rates_of_progress)   # stoichiometric coefficients as orders, forward only
    assert np.allclose(wd_stoich, W * ((nu2 - nu1) @ (kf * np.prod(c[:, None] ** nu1, axis=0))), rtol=1e-12, atol=0)
    for arr in ((T,), g.density * g.Y, wd_orders, wd_nint, wd_stoich):
        lines.append(' '.join('%.16e' % v for v in arr))
open(os.path.join(dst, 'reference.txt'), 'w').write('\n'.join(lines) + '\n')
print('written', dst, 'reactions', g.n_reactions, 'irreversible', sum(1 for r in g.reactions() if not r.reversible), 'orders r1', g.reactions()[0].orders)
