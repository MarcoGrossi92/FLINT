# FLINT Input File Format

FLINT supports **two distinct input paths** for defining chemical mechanisms and thermodynamic data, depending on your use case and whether Cantera is available.

## Quick Navigation

### 🚀 **Choose Your Path:**

| Path | Format | Best For | See |
|------|--------|----------|-----|
| **YAML** | Human-readable | Development, testing, flexibility | [INPUT_YAML_FORMAT.md](INPUT_YAML_FORMAT.md) |
| **Native** | Pre-tabulated | Production, high-performance, no Cantera | [INPUT_NATIVE_FORMAT.md](INPUT_NATIVE_FORMAT.md) |

---

## Path 1: YAML Format (Cantera-Enabled)

**Use this if:**
- ✅ You have Cantera installed
- ✅ You're developing or testing mechanisms
- ✅ You need to modify reaction rates or species properties
- ✅ You want human-readable input files
- ✅ You prioritize flexibility over raw speed

**Features:**
- Single file (`.yaml`) format
- Complete specification in one place
- Easy to edit and version control
- Full Cantera interoperability

**→ [Read Full YAML Specification](INPUT_YAML_FORMAT.md)**

---

## Path 2: Native Files (Non-Cantera)

**Use this if:**
- ✅ You have a production-validated mechanism
- ✅ You need maximum performance (no parsing overhead)
- ✅ Cantera is not available or not desired
- ✅ You want pre-tabulated data for speed

**Features:**
- Four to seven pre-computed files (`.txt`, `.ini`, `.dat`)
- Fast I/O with pre-tabulated properties
- Production-ready

**→ [Read Full Native Format Specification](INPUT_NATIVE_FORMAT.md)**

---

## Comparison Table

| Aspect | YAML | Native |
|--------|------|--------|
| **Cantera Required** | Yes | No |
| **Human-Editable** | Yes | No |
| **Number of Files** | 1 | 4/7 |
| **File Size** | Small | Large |
| **Flexibility** | High (modify any parameter) | Low (regenerate from YAML) |
| **Use Case** | Development | Production |

---

## Documentation Structure

```
INPUT_FORMATS.md (this file)
├── INPUT_YAML_FORMAT.md
│   ├── Units specification
│   ├── Phases configuration
│   ├── Species definitions (thermo & transport)
│   ├── Reactions
│   └── Complete YAML example
├── INPUT_NATIVE_FORMAT.md
│   ├── File structure overview
│   ├── Detailed format for each of 6 files
│   ├── Examples (CORIA mechanism)
│   └── Usage workflow
├── NATIVE_FILES_SUMMARY.md
│   └── Quick reference tables and verification
```

---

## File Structure Overview

### YAML Format
```
project/
└── INPUT/
    └── mechanism.yaml              # Complete definition
```

### Native Format
```
project/
└── INPUT/
    ├── phase.txt               # Species list
    ├── chemistry-info.txt      # Chemistry - Metadata & stoichiometry
    ├── chemistry-Arrhenius.dat # Chemistry - Reaction rates
    ├── chemistry-Troe.dat      # Chemistry - Fall-off parameters (optional)
    ├── chemistry-Lindemann.dat # Chemistry - Fall-off parameters (optional)
    ├── transport.dat           # Transport properties (optional)
    └── thermo.dat              # Thermodynamic properties
```

---

## Quick Examples

### YAML Example (snippet)
```yaml
units:
  energy: cal
  activation-energy: cal/mol

species:
  - name: H2
    composition: {H: 2}
    thermo:
      model: NASA7
      data: [...]
    transport:
      diameter: 2.92
      well-depth: 38.0

reactions:
  - reactants: "H2 + O2"
    products: "2 OH"
    rate-coefficients:
      A: 1.7e13
      b: 0.0
      Ea: 48000.0
```

### Native Example (snippet)
```
# phase.txt
ideal-gas phase
H2      2.016000
O2     31.998000
OH     17.007000

# chemistry-Arrhenius.dat (Tecplot format)
TITLE = "Arrhenius data"
VARIABLES = "Temperature", "kf", "kb"
ZONE T=Reaction1
I=15000, F=POINT
1.0  0.0e+00  0.0e+00
2.0  1.5e-30  2.3e-40
...
```

---
