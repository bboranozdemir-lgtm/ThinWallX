# Sectalix release guide

This document summarizes the release workflow for the Sectalix Python package.

## Local verification

Before preparing a release, run the complete test suite:

```bash
python -m pytest -W error
```

The continuous-integration workflow runs the same test suite on Ubuntu, Windows, and macOS with Python 3.10 through 3.13.

A successful local run is useful but does not replace the cross-platform CI matrix.

## Build distributions

Sectalix uses the setuptools PEP 517 backend. Build tooling should include a compatible setuptools version and `wheel`.

For a manual local build, use a clean environment and create both a wheel and a source distribution. The resulting archives should be checked for:

- the expected package version
- the `sectalix` console entry point
- packaged theory, CLI, verification, and schema files
- absence of obsolete package paths
- successful installation into a clean temporary environment

The repository packaging tests exercise these checks automatically.

## GitHub Actions publication

The repository contains a tag-triggered publication workflow. The intended release sequence is:

1. Ensure the release commit is on `main` and CI is green.
2. Confirm that the package version in `pyproject.toml` matches the intended release tag.
3. Create and push the release tag.
4. Allow the workflow to rebuild the distributions from the tagged source.
5. Review the package metadata and workflow result before publication.

The publication workflow uses GitHub's OpenID Connect / Trusted Publishing path rather than storing a long-lived PyPI API token in the repository.

## Versioning and compatibility

Sectalix v1.0.1 followed the earlier ThinWallX naming used during development. The package, imports, CLI, and public repository use the Sectalix name.

The v0.8 JSON wire-format identifier `"format": "thinwallx"` is retained intentionally so that previously created files remain readable. Changing that identifier would be a separate interchange-format compatibility decision rather than a cosmetic rename.

## Release checks

Before publishing a new version, verify at minimum:

- `python -m pytest -W error` passes locally
- the GitHub Actions test matrix passes
- package metadata matches the intended version
- wheel and source distribution build successfully
- a clean environment can install the built wheel and invoke `sectalix --version`
- documentation links resolve
- examples used in the README still match the current public API
- no generated local files, credentials, or development-only notes are included in the release

## Notes on numerical changes

A packaging or documentation release should not silently change numerical mechanics. If a release modifies section-analysis formulas, topology logic, numerical tolerances, or stress recovery, the change should include:

- a clear technical explanation
- a regression test for the previous failure or limitation
- an independent benchmark or reference case where practical
- an update to the theory or verification documentation when the public model changes

This separation helps distinguish release engineering changes from changes to the underlying mechanics.
