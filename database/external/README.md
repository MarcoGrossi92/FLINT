# External mechanisms

This folder holds chemical mechanisms that are kept out of FLINT's repository, for example because they are
sensitive material. Git ignores everything in it except this README and the submodule `TSR-Family`. FLINT compiles
the mechanisms it finds here and hooks them in `Assign_Mechanism` like its own, with no change to its sources.

## Layout

One folder per family:

```
database/external/<family>/
  <mechanism>/                     tables (phase.txt, thermo, chemistry-*, chemistry-info.txt) and Cantera file
  FLINT/src/*.f90                  routines, as utils/YTF.py writes them (module <X>_mod, subroutine <X>)
  FLINT/mechanism_contract.json    their records of the mechanism contract (format of
                                   src/lib/Lib_ChemMech/mechanism_contract.json)
```

The configure compiles `FLINT/src/*.f90` of every family and writes the module `FLINT_Lib_Chemistry_external`
into the build folder from the records, with `python3 utils/mechanism_contract.py`. That module maps each mechanism
name to its routine and gives `check_mechanism_contract` its record. Without any family, FLINT compiles the stub
`src/lib/external/Lib_Chemistry_external.f90`, which hooks nothing.

A record names the routine, its species slots and its reaction counts, and a mechanism name must not be one that
FLINT hooks already. To write the fingerprints of the routines into the records, and to check them:

```bash
python3 utils/mechanism_contract.py --fingerprints
python3 utils/mechanism_contract.py --check
```

`-DFLINT_EXTERNAL_DIR=<path>` takes the families from another folder.

## TSR mechanisms

`TSR-Family` is the private submodule of the TSR mechanisms (`TSR-CDF-13`, `TSR-GP-24`, `TSR-Rich-31`), available
on request. With access, check it out from the root of FLINT and configure:

```bash
git submodule update --init --checkout database/external/TSR-Family
git config submodule.TSR-Family.update checkout   # optional: git submodule update keeps it in step from now on
```
