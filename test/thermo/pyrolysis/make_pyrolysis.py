#!/usr/bin/env python3
"""Generate the fixtures of test-pyrolysis-gas (the pyrolysis-gas tables of FLINT_Lib_PyrolysisGas):
  good.dat            three tables, read and mapped onto the solver list of the driver  -> accepted (ios = 0)
                       (a full-line comment on line 1, a trailing one on the line 'condensed 0' of bilin)
                       bilin   gas-only, 3 pressures x 4 temperatures, species CO H2 C2H2,acetylene H2O: every
                               column is a + b P + c Q + d P Q, P = (p - 1e5)/2e5, Q = (T - 1000)/500, with the
                               dyadic coefficients COEF below, so that the bilinear interpolation reproduces it
                               (exactly at the points of the driver where P and Q are multiples of 1/4)
                       single  strip, one pressure (2e5 Pa) x 3 temperatures, species CO H2 CH4, CH4 = 0 and
                               not a species of the solver (accepted: its column is identically 0)
                       cond    condensed, condensed 1, 2 pressures x 2 temperatures, species CO H2 + X_C
  bad-format.dat      'format 2'                                                       -> refused (ios = 3)
  bad-species.dat     bilin with species XYZ (not in the solver list) in place of H2O  -> refused (ios = 5)
  bad-grid.dat        bilin with the temperatures 1000 1500 1400 3000                  -> refused (ios = 6)
  bad-sum.dat         bilin with row 5 summing to 1 + 1e-6                             -> refused (ios = 7)
  bad-short.dat       bilin without its last data row                                  -> refused (ios = 8)
  bad-extra.dat       bilin with a fifth value (0) on every row, with condensed 0      -> refused (ios = 8)
  bad-noend.dat       bilin without 'end table'                                        -> refused (ios = 9)
  bad-condflag.dat    bilin with 'condensed 1' and treatment gas-only                  -> refused (ios = 4)
  bad-makeup.dat      bilin with elemental mass fractions summing to 0.99              -> refused (ios = 4)
  bad-negative.dat    bilin with 0.25 moved from H2O to CO in row 5: H2O -0.04296875,
                      the row still summing to 1 (dyadic: exactly)                     -> refused (ios = 10)
  bad-dupname.dat     bilin, then the table single named bilin                         -> refused (ios = 11)
The values are written as ATLAS GPB writes them (ES24.16, Fortran's E format: 1.6 digits after the point).
Run from the repository root:  python3 test/thermo/pyrolysis/make_pyrolysis.py
"""
import os
import sys

root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
dst = os.path.join(root, 'test', 'thermo', 'pyrolysis')
if len(sys.argv) > 1:          # another folder (a scratch check of the generator)
    dst = sys.argv[1]

# bilin: column k = a + b P + c Q + d P Q (keep in step with COEF of src/test/thermo/test-pyrolysis-gas.f90)
SPECIES = ['CO', 'H2', 'C2H2,acetylene', 'H2O']
COEF = [(0.40625, 1/64, -1/128, 1/512),       # CO
        (0.125, -1/128, 1/256, -1/1024),       # H2
        (0.25, 1/256, 1/64, -1/512)]           # C2H2,acetylene; H2O = 1 - the three
COEF.append(tuple((1.0 if i == 0 else 0.0) - sum(c[i] for c in COEF) for i in range(4)))
P_GRID = [1.0e5, 3.0e5, 7.0e5]
T_GRID = [1000.0, 1500.0, 2500.0, 3000.0]


def e(x):
    """ES24.16 as Fortran writes it."""
    s = f'{x:.16E}'
    mant, exp = s.split('E')
    return f'{mant}E{int(exp):+03d}'.rjust(24)


def bilin(p, T):
    P, Q = (p - 1.0e5) / 2.0e5, (T - 1000.0) / 500.0
    return [a + b * P + c * Q + d * P * Q for a, b, c, d in COEF]


def table(name, treatment, species, pgrid, tgrid, rows, condensed=0, makeup=(('C', 0.5), ('H', 0.25), ('O', 0.25)),
          end=True):
    out = [f'table {name}', f'treatment {treatment}',
           f'elements {len(makeup)} ' + ' '.join(f'{el} {e(m).strip()}' for el, m in makeup),
           f'condensed {condensed}', f'species {len(species)}'] + species
    out += [f'pressures {len(pgrid)}', ' '.join(e(p) for p in pgrid),
            f'temperatures {len(tgrid)}', ' '.join(e(T) for T in tgrid), 'data']
    out += [' '.join(e(v) for v in r) for r in rows]
    if end:
        out.append('end table')
    return out


def bilin_rows(pgrid=P_GRID, tgrid=T_GRID):
    return [bilin(p, T) for p in pgrid for T in tgrid]      # pressure-major


def write(fname, blocks, fmt='1'):
    lines = ['# hydra pyrolysis-gas composition table', f'format {fmt}', f'tables {len(blocks)}']
    for b in blocks:
        lines += b
    with open(os.path.join(dst, fname), 'w') as f:
        f.write('\n'.join(lines) + '\n')


def main():
    os.makedirs(dst, exist_ok=True)
    for r in bilin_rows():
        assert all(v > 0.0 for v in r) and abs(sum(r) - 1.0) <= 1e-15, r
    bil = table('bilin', 'gas-only', SPECIES, P_GRID, T_GRID, bilin_rows())
    # single: y_CO = 0.75 - Q/16, Q = (T - 1000)/1000; y_H2 = 1 - y_CO; CH4 = 0 (keep in step with the driver)
    srows = [[0.75 - q / 16.0, 0.25 + q / 16.0, 0.0] for q in (0.0, 1.0, 2.0)]
    sng = table('single', 'strip', ['CO', 'H2', 'CH4'], [2.0e5], [1000.0, 2000.0, 3000.0], srows)
    # cond: X_C = 0.25 + R/8, y_CO = 0.5 - R/16, y_H2 = 0.25 - R/16, R = 2 (row of the 2nd temperature) + (p index - 1)
    crows = []
    for ip in range(2):
        for iT in range(2):
            r = iT * 2 + ip
            crows.append([0.5 - r / 16.0, 0.25 - r / 16.0, 0.25 + r / 8.0])
    cnd = table('cond', 'condensed', ['CO', 'H2'], [1.0e5, 2.0e5], [2000.0, 3000.0], crows, condensed=1)
    # a comment after a value (a '#' after a blank) on the condensed line of bilin
    write('good.dat', [[l + '   # no X_C column' if l == 'condensed 0' else l for l in bil], sng, cnd])

    write('bad-format.dat', [bil], fmt='2')
    write('bad-species.dat', [table('bilin', 'gas-only', SPECIES[:3] + ['XYZ'], P_GRID, T_GRID, bilin_rows())])
    write('bad-grid.dat', [table('bilin', 'gas-only', SPECIES, P_GRID, [1000.0, 1500.0, 1400.0, 3000.0],
                                 bilin_rows())])
    rows = bilin_rows()
    rows[4] = rows[4][:3] + [rows[4][3] + 1.0e-6]
    write('bad-sum.dat', [table('bilin', 'gas-only', SPECIES, P_GRID, T_GRID, rows)])
    write('bad-short.dat', [table('bilin', 'gas-only', SPECIES, P_GRID, T_GRID, bilin_rows()[:-1])])
    write('bad-extra.dat', [table('bilin', 'gas-only', SPECIES, P_GRID, T_GRID, [r + [0.0] for r in bilin_rows()])])
    write('bad-noend.dat', [table('bilin', 'gas-only', SPECIES, P_GRID, T_GRID, bilin_rows(), end=False)])
    write('bad-condflag.dat', [table('bilin', 'gas-only', SPECIES, P_GRID, T_GRID, bilin_rows(), condensed=1)])
    write('bad-makeup.dat', [table('bilin', 'gas-only', SPECIES, P_GRID, T_GRID, bilin_rows(),
                                   makeup=(('C', 0.5), ('H', 0.24), ('O', 0.25)))])
    rows = bilin_rows()
    rows[4] = [rows[4][0] + 0.25] + rows[4][1:3] + [rows[4][3] - 0.25]
    write('bad-negative.dat', [table('bilin', 'gas-only', SPECIES, P_GRID, T_GRID, rows)])
    write('bad-dupname.dat', [bil, table('bilin', 'strip', ['CO', 'H2', 'CH4'], [2.0e5], [1000.0, 2000.0, 3000.0], srows)])


if __name__ == '__main__':
    main()
