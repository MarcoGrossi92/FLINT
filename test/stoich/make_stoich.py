#!/usr/bin/env python3
"""Generate test/stoich: mechanisms with FRACTIONAL stoichiometric coefficients and no explicit orders, for test-stoich.
(Cantera flags proportional equations as duplicates: they carry `duplicate: true`.) Cantera's law is mass action with the real stoichiometric coefficients: forward [A]^nu'_A..., reverse [P]^nu''_P... with
kb = kf/Kc (Kc of the same, possibly fractional, reaction). Five mechanisms (SI units; thermo and transport of Cantera's
h2o2.yaml):
  stoich-frac   H2 + 0.5 O2 <=> H2O;  0.5 H2 + 0.25 O2 <=> 0.5 H2O      (fractional reactants; nint(0.25) = 0)
  stoich-prod   H2O <=> H2 + 0.5 O2;  1.5 H2O <=> 1.5 H2 + 0.75 O2       (fractional products; nint(1.5) = 2)
  stoich-3b     H2 + 0.5 O2 + M <=> H2O + M                              (three-body)
  stoich-troe   H + O2 <=> O + OH;  1.5 H2 + 0.75 O2 (+M) <=> 1.5 H2O (+M)   (Troe, both sides fractional; states 1000..2100 K)
  stoich-lind   H + O2 <=> O + OH;  1.5 H2O (+M) <=> 1.5 H2 + 0.75 O2 (+M)   (Lindemann, both sides fractional; states 1000..2100 K)
  stoich-int    2 H2 + O2 <=> 2 H2O  (control: integer coefficients, the Arrhenius numbers of reaction 1 of stoich-frac)
(every fixture directory stays below 200 KB: at most two 1 K rate tables on 1000..2501 K, or 1000..2101 K with a falloff one)
Steps (run from the repository root, Cantera 3.0.1):
  1) python3 test/stoich/make_stoich.py yaml             writes test/stoich/yaml/<name>.yaml
  2) write the INPUT tables of each <name>.yaml with a table writer, on a 1 K step from Tmin = 1000 K to
     Tmax = 2501 K (2101 K for stoich-troe and stoich-lind: one row above the highest reference state, since the
     rate at T interpolates the rows T and T + 1), with the thermo of the yaml: the tables' kb = kf/Kc then come
     from the thermo that Cantera uses below
  3) python3 test/stoich/make_stoich.py fixture <name>=<its INPUT folder> ...
     copies phase.txt, chemistry-info.txt and chemistry-{Arrhenius,Troe,Lindemann}.dat into test/stoich/<name>/ and writes
     reference.txt: species names, molecular weights, then per state T [K], p [Pa], flag (1 = Cantera equilibrium composition),
     the partial densities [kg/m3], Cantera's net production rates [kg/m3/s] and the gross rates W*(creation + destruction)
     [kg/m3/s] (the scale of the comparison). States: T integer (a table node) in 1000..2500 K, 0.5..10 atm, lean and rich
     mixtures, compositions within 1e-3 of equilibrium, and equilibrium compositions (equilibrate TP). The script asserts that
     the tables carry Cantera's kf and kb (Arrhenius) and k0, kinf, Kc (falloff) at every reference temperature (relative 5e-12).
"""
import os, re, sys, numpy as np, cantera as ct
here = os.path.dirname(os.path.abspath(__file__))
T0, T1 = 1000, 2500
SP4 = '[H2, O2, H2O, N2]'; SP7 = '[H2, O2, H2O, N2, H, O, OH]'
R1 = '{A: 1.0e+10, b: 0.0, Ea: 1.2e+08}'
HO2 = '- {equation: H + O2 <=> O + OH, rate-constant: {A: 3.52e+13, b: -0.7, Ea: 7.1e+07}}\n'
MECH = {
 'stoich-frac': (SP4, '- {equation: H2 + 0.5 O2 <=> H2O, rate-constant: %s, duplicate: true}\n' % R1 +
                 '- {equation: 0.5 H2 + 0.25 O2 <=> 0.5 H2O, rate-constant: {A: 1.0e+08, b: 0.3, Ea: 1.1e+08}, duplicate: true}\n'),
 'stoich-prod': (SP4, '- {equation: H2O <=> H2 + 0.5 O2, rate-constant: {A: 1.0e+12, b: 0.5, Ea: 4.0e+08}, duplicate: true}\n'
                 '- {equation: 1.5 H2O <=> 1.5 H2 + 0.75 O2, rate-constant: {A: 1.0e+11, b: 0.2, Ea: 3.8e+08}, duplicate: true}\n'),
 'stoich-3b':   (SP4, '- {equation: H2 + 0.5 O2 + M <=> H2O + M, type: three-body, rate-constant: {A: 1.0e+12, b: -0.5, Ea: 1.0e+08},'
                 ' efficiencies: {H2O: 6.0, N2: 0.8}}\n'),
 'stoich-troe': (SP7, HO2 + '- {equation: 1.5 H2 + 0.75 O2 (+M) <=> 1.5 H2O (+M), type: falloff, low-P-rate-constant: {A: 1.0e+13, b: 0.0,'
                 ' Ea: 1.0e+08}, high-P-rate-constant: %s, Troe: {A: 0.5, T3: 100.0, T1: 2000.0, T2: 5000.0},'
                 ' efficiencies: {H2O: 6.0}}\n' % R1),
 'stoich-lind': (SP7, HO2 + '- {equation: 1.5 H2O (+M) <=> 1.5 H2 + 0.75 O2 (+M), type: falloff, low-P-rate-constant: {A: 1.0e+15, b: 0.0,'
                 ' Ea: 3.5e+08}, high-P-rate-constant: {A: 1.0e+12, b: 0.5, Ea: 4.0e+08}, efficiencies: {H2O: 6.0}}\n'),
 'stoich-int':  (SP4, '- {equation: 2 H2 + O2 <=> 2 H2O, rate-constant: %s}\n' % R1),
}
def build(name):
    sp, rx = MECH[name]
    return ct.Solution(yaml='units: {length: m, quantity: kmol, activation-energy: J/kmol}\nphases:\n- name: %s\n  thermo: ideal-gas\n'
        '  elements: [H, O, N]\n  species: [{h2o2.yaml/species: %s}]\n  kinetics: gas\n  reactions: [R]\n  transport: mixture-averaged\n'
        '  state: {T: 300.0, P: 1 atm}\nR:\n%s' % (name, sp, rx))
def zones(path):
    txt = open(path).read(); idx = [m.start() for m in re.finditer(r'(?m)^ZONE T=', txt)]
    return [[[float(v) for v in r.split()] for r in txt[idx[k]:(idx[k+1] if k+1 < len(idx) else len(txt))].split('\n')[2:] if r.strip()]
            for k in range(len(idx))]
if sys.argv[1] == 'yaml':
    os.makedirs(os.path.join(here, 'yaml'), exist_ok=True)
    for name in MECH:
        build(name).write_yaml(os.path.join(here, 'yaml', name + '.yaml'))
    sys.exit(0)
for arg in sys.argv[2:]:
    name, src = arg.split('=', 1)
    g = ct.Solution(os.path.join(here, 'yaml', name + '.yaml'))
    ref = build(name)
    assert [str(r) for r in g.reactions()] == [str(r) for r in ref.reactions()]
    dst = os.path.join(here, name); os.makedirs(dst, exist_ok=True)
    for f in ('phase.txt', 'chemistry-info.txt', 'chemistry-Arrhenius.dat', 'chemistry-Troe.dat', 'chemistry-Lindemann.dat'):
        if os.path.exists(os.path.join(src, f)):
            open(os.path.join(dst, f), 'w').write(open(os.path.join(src, f)).read())
    names = [l.split()[0] for l in open(os.path.join(dst, 'phase.txt')).read().split('\n')[1:] if l.strip()]
    assert names == list(g.species_names), (names, g.species_names)
    ZA = zones(os.path.join(dst, 'chemistry-Arrhenius.dat'))
    ZT = zones(os.path.join(dst, 'chemistry-Troe.dat')) if os.path.exists(os.path.join(dst, 'chemistry-Troe.dat')) else []
    ZL = zones(os.path.join(dst, 'chemistry-Lindemann.dat')) if os.path.exists(os.path.join(dst, 'chemistry-Lindemann.dat')) else []
    rad = {'H': 2e-4, 'O': 1e-4, 'OH': 1e-3} if 'OH' in names else {}
    lean = dict({'H2': 0.10, 'O2': 0.20, 'H2O': 0.15, 'N2': 0.55}, **rad); rich = dict({'H2': 0.30, 'O2': 0.05, 'H2O': 0.10, 'N2': 0.55}, **rad)
    states = []
    Tmax = 2100 if name in ('stoich-troe', 'stoich-lind') else 2500
    for k, T in enumerate(range(1000, Tmax + 1, 150)):
        for m, X in enumerate((lean, rich)):
            states.append((float(T), [0.5, 1.0, 2.0, 5.0, 10.0][(2*k + m) % 5] * ct.one_atm, X, 0))
    for T, p in ((1200., ct.one_atm), (1700., 2*ct.one_atm), (Tmax - 300., 0.5*ct.one_atm), (float(Tmax), 10*ct.one_atm)):
        g.TPX = T, p, lean; g.equilibrate('TP'); Xe = g.X * (1 + 1e-3*np.where(np.arange(g.n_species) % 2 == 0, 1.0, -1.0))
        states.append((T, p, Xe / Xe.sum(), 0))
    for T, p, X in ((1100., ct.one_atm, lean), (1600., 5*ct.one_atm, rich), (Tmax - 400., ct.one_atm, lean), (float(Tmax), 0.5*ct.one_atm, rich)):
        g.TPX = T, p, X; g.equilibrate('TP'); states.append((T, p, g.X.copy(), 1))
    W = g.molecular_weights
    lines = ['# %s: Cantera %s, tables of a table writer; names, W [kg/kmol]; per state: T p flag / roi [kg/m3] / net wdot [kg/m3/s] / gross W*(creation+destruction) [kg/m3/s]'
             % (name, ct.__version__), '%d %d' % (g.n_species, len(states)), ' '.join(names),
             ' '.join('%.16e' % v for v in W)]
    worst_eq = 0.0
    for T, p, X, flag in states:
        g.TPX = T, p, X
        # the tables must carry Cantera's constants at this node (Arrhenius kf, kb; falloff k0, kinf, Kc)
        ia = it = il = 0
        for i, r in enumerate(g.reactions()):
            if r.reaction_type == 'falloff-Troe':
                row = [z for z in ZT[it] if z[0] == T][0]; it += 1
                for a, b in ((row[1], r.rate.high_rate(T)), (row[2], r.rate.low_rate(T)), (row[3], g.equilibrium_constants[i])):
                    assert abs(a - b) <= 5e-12 * abs(b), (name, T, i, a, b)
            elif r.reaction_type == 'falloff-Lindemann':
                row = [z for z in ZL[il] if z[0] == T][0]; il += 1
                for a, b in ((row[1], r.rate.high_rate(T)), (row[2], r.rate.low_rate(T)), (row[3], g.equilibrium_constants[i])):
                    assert abs(a - b) <= 5e-12 * abs(b), (name, T, i, a, b)
            else:
                row = [z for z in ZA[ia] if z[0] == T][0]; ia += 1
                for a, b in ((row[1], g.forward_rate_constants[i]), (row[2], g.reverse_rate_constants[i])):
                    assert abs(a - b) <= 5e-12 * abs(b), (name, T, i, a, b)
        wd = W * g.net_production_rates; gr = W * (g.creation_rates + g.destruction_rates)
        if flag: worst_eq = max(worst_eq, np.max(np.abs(wd)) / np.max(gr))
        lines += ['%.1f %.10e %d' % (T, p, flag)] + [' '.join('%.16e' % v for v in a) for a in (g.density * g.Y, wd, gr)]
    open(os.path.join(dst, 'reference.txt'), 'w').write('\n'.join(lines) + '\n')
    print('%s: %d states, %d reactions (%s); Cantera residual at its equilibrium states max|wdot|/max(gross) = %.1e'
          % (name, len(states), g.n_reactions, ', '.join(r.reaction_type for r in g.reactions()), worst_eq))
