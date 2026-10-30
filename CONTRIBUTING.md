# Contributing to FLINT

This document provides guidelines and instructions for contributing to the project.

## Ways to Contribute

There are many ways to contribute to FLINT:

- 🐛 **Report bugs** - Help us identify and fix issues
- 💡 **Suggest features** - Share ideas for improvements
- 📝 **Improve documentation** - Fix typos, clarify explanations, add examples
- 🔧 **Submit code** - Fix bugs or implement new features
- 🧪 **Add tests** - Improve test coverage
- 🎨 **Share examples** - Contribute usage examples

## Getting Help

Questions about installing or using FLINT go to the [issue tracker](https://github.com/MarcoGrossi92/FLINT/issues),
as bug reports and feature requests do: open an issue with the template **Question**, **Bug report** or **Feature
request**. Look through the existing issues and the [documentation](https://marcogrossi92.github.io/FLINT/) first.

## Getting Started

### 1. Fork the Repository

Click the "Fork" button on [GitHub](https://github.com/MarcoGrossi92/FLINT) to create your own copy.

### 2. Clone Your Fork

```bash
git clone https://github.com/YOUR_USERNAME/FLINT.git
cd FLINT
git remote add upstream https://github.com/MarcoGrossi92/FLINT.git
```

### 3. Build and Test

```bash
./install.sh build --compilers=gnu
cd build && ctest -L quick
```

The first configure fetches ORION and OSLO from GitHub, so it needs the network; the next ones use the releases
recorded in `deps.lock`. See [Installation](https://marcogrossi92.github.io/FLINT/getting-started/installation/)
for the options (TecIO, SUNDIALS, Cantera) and [Testing](https://marcogrossi92.github.io/FLINT/development/testing/)
for the test suite.

### 4. Create a Branch

```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/your-bug-fix
```

Use descriptive branch names:

- `feature/add-hdf5-support`
- `fix/tecplot-binary-crash`
- `docs/improve-installation`

## Development Workflow

### Making Changes

1. **Write clear, focused commits.** One change per commit, with a subject that names the part of FLINT it touches
   and a body that says why the change is needed and what it changes:
   ```
   Loaders: refuse a table whose temperature step is not 1 K

   Why: a table on a 2 K step was accepted and row T held the value at 2T - 1 K.
   Before and after: such a table is refused with an [ERROR] line and ios = 6;
   tables on a 1 K step are read as before.
   ```

2. **Follow coding standards:**
   - Fortran: free-form format, meaningful variable names, comments for complex logic
   - Use 2 spaces for indentation in Fortran
   - Write new code like the code around it
   - FLINT is built with gfortran and with ifx: a change to the sources or to the build must compile with both

3. **Add tests for new features and fixes.** A test is a driver in `src/test/<area>/`, registered in
   `src/test/CMakeLists.txt` with `flint_add_test`, with its fixtures in `test/<area>/`. Run `ctest -L quick` before
   every pull request, and the whole suite (`ctest`) after a change to the loaders, the thermodynamics or the
   kinetics.

4. **Update documentation:**
   - Add docstrings to new functions
   - Update relevant .md files in `docs/`
   - Add examples if appropriate

5. **Add an entry to the changelog.** Every change that a user of FLINT sees goes under `[Unreleased]` in
   [CHANGELOG.md](https://github.com/MarcoGrossi92/FLINT/blob/main/CHANGELOG.md), in its section (Added, Changed,
   Deprecated, Removed, Fixed). Say what changes for the user and, for a change of results, which results and by
   how much. Begin the entry with **Breaking.** when a project that uses FLINT has to change with it;
   [Semantic Versioning](https://marcogrossi92.github.io/FLINT/development/versioning/) says which changes are.

### Working on ORION or OSLO with FLINT

FLINT builds releases of ORION and OSLO. To build a local checkout instead, for a change that touches both:

```bash
./install.sh build --include-orion=../ORION     # or --include-oslo=../OSLO
# with CMake: -DFETCHCONTENT_SOURCE_DIR_ORION=/path/to/ORION
```

The checkout is built as it is, uncommitted changes included. The change then reaches FLINT through a release:

1. Merge the change into ORION and publish a release (a tag `vX.Y.Z`).
2. If FLINT needs it, raise the requirement in FLINT's `CMakeLists.txt` (`VERSION 2.2`), with a comment that says
   why, in the pull request of the FLINT change. The GitHub workflow builds releases, not your checkout, so that
   pull request passes only once the release exists.
3. Go back to the releases: `cmake -UFETCHCONTENT_SOURCE_DIR_ORION <build>`, or a build without `--include-orion`.

### Keeping Your Fork Updated

```bash
# Fetch upstream changes
git fetch upstream

# Merge into your branch
git checkout main
git merge upstream/main

# Update your fork
git push origin main
```

## Submitting Changes

**1. Push Your Branch**

```bash
git push origin feature/your-feature-name
```

**2. Create a Pull Request**

1. Go to your fork on GitHub
2. Click "Pull Request"
3. Select your branch
4. Fill out the PR template

**Pull Request Guidelines**

A good pull request includes:

- **Clear title:** "Add HDF5 format support" not "Update files"
- **Description:** What changes were made and why
- **References:** Link to related issues (`Fixes #123`)
- **Tests:** Proof that changes work
- **Documentation:** Updates to relevant docs
- **Changelog:** An entry under `[Unreleased]`

The GitHub workflow runs `ctest -L quick` on every pull request.

**PR Template:**

```markdown
## Description
Brief description of changes

## Motivation
Why is this change needed?

## Changes Made
- Added X feature
- Fixed Y bug
- Updated Z documentation

## Testing
How were these changes tested?

## Checklist
- [ ] `ctest -L quick` passes (and `ctest`, for changes to loaders, thermodynamics or kinetics)
- [ ] Documentation updated
- [ ] CHANGELOG.md entry under [Unreleased]
- [ ] Code follows style guidelines
- [ ] Commit messages are clear
```

**Code Review Process**

1. Maintainers will review your PR
2. Address any requested changes
3. Push updates to the same branch
4. Once approved, your PR will be merged!

Maintainers publish the releases, as described in
[Semantic Versioning](https://marcogrossi92.github.io/FLINT/development/versioning/).

## Reporting Issues

The [issue tracker](https://github.com/MarcoGrossi92/FLINT/issues/new/choose) has a template for each kind of issue,
which asks for what follows.

### Bug Reports

When reporting bugs, include:

1. **System information:**
   - OS and version
   - Compiler and version
   - Build configuration
   - The releases of ORION and OSLO (in `deps.lock`)

2. **Steps to reproduce:**

3. **Expected vs actual behavior:**
   - What should happen
   - What actually happens
   - Error messages or output

4. **Sample data** (if possible):
   - Minimal example that reproduces the issue

### Feature Requests

For feature requests, describe:

1. **Use case:** Why is this feature needed?
2. **Proposed solution:** How should it work?
3. **Alternatives:** What alternatives exist?
4. **Impact:** Who would benefit?

## Recognition

Significant contributions may earn you:

- Commit access
- Listed as project contributor in
  [AUTHORS.md](https://github.com/MarcoGrossi92/FLINT/blob/main/AUTHORS.md)

We appreciate your contributions! 🎉
