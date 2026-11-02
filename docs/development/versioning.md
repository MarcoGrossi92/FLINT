# Semantic Versioning

FLINT, ORION and OSLO number their releases `MAJOR.MINOR.PATCH`, following
[Semantic Versioning 2.0.0](https://semver.org/spec/v2.0.0.html). A release is a tag `vX.Y.Z` on `main`. The number
tells a project that uses FLINT what it can take without a change of its own: a minor or a patch release reaches it
with no change, a major release may need one.

This page says what each kind of release of FLINT may change, how FLINT requires its own libraries, and how a release
is made. The changes of each release are listed in
[CHANGELOG.md](https://github.com/MarcoGrossi92/FLINT/blob/main/CHANGELOG.md).

## What a release may change

A project that uses FLINT relies on its public interface:

- **The Fortran interface** that the [User Guide](../user/index.md) describes: the modules, the procedures, their
  arguments and their return codes, and the module variables that a solver reads (the species, the tables).
- **The input files**: the files of an INPUT folder in the [native](../user/input/native.md) and
  [YAML](../user/input/yaml.md) formats, what their content means, and the mechanism names that select a compiled
  routine.
- **The build interface**: the target `FLINT` (`FLINT::FLINT`), the list `FLINT_SRC` and `link_FLINT_dependencies`,
  the options `USE_CANTERA`, `USE_OPENMP` and `USE_MPI`, and the commands and options of `install.sh`.
- **The requirements**: the releases of ORION and OSLO that FLINT accepts.

The test drivers, their labels and their fixtures, and what the documentation does not describe, are not part of it.

| Release | When | For example |
|---|---|---|
| Major, 2.1.0 → 3.0.0 | A project that uses FLINT may have to change: something it relies on is removed, renamed or means something else. | Remove or rename a procedure, an argument, a module variable, a target or an option. Change the meaning of a return code or of an input field. Refuse an input that the previous release read. Require a new major release of ORION or OSLO. |
| Minor, 2.1.0 → 2.2.0 | Something new, and every project keeps working. | Add a procedure, an optional argument, a mechanism, an optional block of an input file, an option. Require a later minor release of ORION or OSLO. Deprecate something, which the next major release may remove. |
| Patch, 2.1.0 → 2.1.1 | A fix, with no new interface. | Correct a wrong result, a crash, the refusal of a valid input. Documentation and tests. |

The changelog marks the entries of a major release with **Breaking.**, and says what a project has to change.

### Changes of results

A fix that changes numerical results, such as the exact gas constant of 1.1.0, is a patch or a minor release, not a
major one: a wrong result is not part of the interface. Its changelog entry says which results change and by how
much, so that a user can tell which computations to run again.

### The releases of ORION and OSLO

A build tree holds one copy of ORION and one of OSLO, whatever the number of projects that use them. A solver that
requires ORION 2 and uses FLINT shares ORION with it, so:

- **A new major release of ORION or OSLO in FLINT's requirement is a major release of FLINT.** If a minor release of
  FLINT required ORION 3, the next configure of that solver would stop, with no change of its own.
- **A later minor release is a minor release of FLINT.** Raising the requirement from ORION 2.1 to 2.2 asks the
  build tree for 2.2.0 or a later 2.x, which satisfies the solver's requirement of ORION 2 too.

FLINT 2.0.0 is a major release for this reason among others: FLINT 1 used ORION 1, FLINT 2 requires ORION 2.

## How FLINT requires its libraries

`CMakeLists.txt` requires ORION and OSLO by a range of releases, with `semver_dependency` of
`cmake/SemVerDeps.cmake`:

```cmake
semver_dependency(ORION GIT_REPOSITORY https://github.com/MarcoGrossi92/ORION.git VERSION 2.1)
```

A requirement accepts the releases from its own up to the next change of its leftmost nonzero number, the one that a
release of the library may not change without breaking its users:

| Requirement | Accepts |
|---|---|
| `2.1` | 2.1.0 and later, up to 3.0.0 excluded |
| `2` | 2.0.0 and later, up to 3.0.0 excluded |
| `2.1.3` | 2.1.3 and later, up to 3.0.0 excluded |
| `0.3` | 0.3.0 and later, up to 0.4.0 excluded (before 1.0.0, a minor release may break) |
| `0.0.3` | 0.0.3 only |
| `1.7 2` | 1.7.0 and later 1.x, or 2.x: the values are alternatives |

Only the tags `vX.Y.Z` (or `X.Y.Z`) are releases. A tag with a suffix, such as `v2.2.0-rc1`, is a pre-release: it
is never selected, and a project tries it by building a local checkout of it (see
[Installation](../getting-started/installation.md#orion-and-oslo)). `GIT_TAG <tag>` in place of `VERSION` pins one
tag, for a library that is not released in this way.

FLINT requires the oldest minor release that has everything it uses, and raises the requirement when it starts to use
something newer, with a comment that says what: `2.1` for ORION because `install.sh` and the GitHub workflow pass
`ORION_TECIO_DIR`, which ORION 2.1.0 added.

### Which release a build gets

The first `semver_dependency` of a library in a build tree selects its release; a later one, from a project added
inside the tree, only checks it, and the configure stops, naming both requirements, when it is not satisfied. The
release is, in order:

1. the local checkout given with `-DFETCHCONTENT_SOURCE_DIR_<NAME>=<folder>`, as it is;
2. the release in `deps.lock`, in the top-level source folder, if it satisfies the requirement;
3. the newest release on GitHub that satisfies the requirement.

The top-level project of a tree requires all its libraries before it makes any of them available, so that it is the
one that selects. [Installation](../getting-started/installation.md#orion-and-oslo) gives the commands: the lock file,
`-DSEMVER_UPDATE` for newer releases, a local checkout.

An application at the top of its build tree may commit its `deps.lock`, so that every build of a commit gets the same
releases. A library such as FLINT does not: in a larger tree the top-level project chooses.

## Requiring FLINT

A project requires FLINT in the same way, with its own copy of `cmake/SemVerDeps.cmake`:

```cmake
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SemVerDeps.cmake)
semver_dependency(FLINT GIT_REPOSITORY https://github.com/MarcoGrossi92/FLINT.git VERSION 2.0
                  GIT_SUBMODULES "")      # lib/cantera is cloned by FLINT, if USE_CANTERA is on
FetchContent_MakeAvailable(FLINT)
```

With `VERSION 2.0` it gets every 2.x release of FLINT from 2.0.0 on, and never a 3.0.0 that could break it.
[Using FLINT in a CMake Project](../getting-started/installation.md#using-flint-in-a-cmake-project) describes the
targets and the flags.

## Making a release

1. Choose the number from the entries under `[Unreleased]` in `CHANGELOG.md`: major if one is **Breaking**,
   otherwise minor if one adds something, otherwise patch.
2. In `CHANGELOG.md`, rename `[Unreleased]` to `[X.Y.Z] - YYYY-MM-DD`, and open a new, empty `[Unreleased]` above it.
3. In `CITATION.cff`, set `version` and `date-released`.
4. Merge into `main` through a pull request.
5. Tag the merge commit and push the tag:
   ```bash
   git fetch origin
   git tag vX.Y.Z origin/main
   git push origin vX.Y.Z
   ```
6. Optionally, publish a GitHub release from the tag, with the section of the changelog as its notes.

A project that uses FLINT gets the release at its next selection: a new build folder, or
`cmake -DSEMVER_UPDATE=FLINT <build>`.

!!! warning "A published tag stays"
    Once pushed, a tag is never moved or deleted: the `deps.lock` files of other projects name it, and a build of
    theirs would get other sources under the same number. A mistake in a release is corrected by the next release.
