# Sectalix v1.0.1 release and Trusted Publishing guide

The historical ThinWallX v1.0.0 tag remains unchanged. Sectalix v1.0.1 is the
rename/release-preparation version; create its tag only after review. Protect main
and release tags on GitHub.

## Local verification and CI

Run `python -m pytest -W error`: the frozen v1.0 core plus the v1.0.1
compatibility regression contains 560 tests.
The CI workflow targets Ubuntu, Windows and macOS with Python 3.10–3.13 (12 jobs).
These hosted jobs can only be confirmed after pushing to GitHub; a local Windows
run does not establish that all matrix combinations pass.

The tag-triggered publish.yml first calls the same full matrix, then builds fresh
wheel/sdist into release-dist (not the tracked historical dist directory), checks
tag/version equality and runs strict Twine metadata validation. Verified archives
are retained as a GitHub Actions artifact for 30 days. There is intentionally no
automatic upload is performed locally; the tag workflow uses OIDC only after the
protected `pypi` environment is approved.

The packaging tests invoke setuptools directly, so CI installs that backend and
wheel explicitly in addition to the requested editable test/plot extras. This is
build/test tooling, not a new application dependency.

The minimum backend is setuptools 77.0.3, which supports the SPDX license-string
format. Checkout and setup-python use their Node 24-based v6 actions; reverting
to checkout v4/setup-python v5 is not a migration away from Node 20.

Project-wide deprecation suppression is disabled. The packaging subprocess exempts
only the known distutils import deprecation message on Windows Python 3.10/3.11.
All other messages remain errors. Passing means no unsuppressed warnings; it does
not imply that an exempt dependency warning was never emitted.

The residual-torque regression uses 32 machine epsilons times the sum of absolute
torque terms, without an absolute floor, and exercises three length and three load
scales. These narrow test maintenance changes were explicitly authorized without
changing the numerical core. Platform rounding differences do not by themselves
prove that FMA caused the original discrepancy.

## Connect a new GitHub repository

1. Create an empty repository in your GitHub account. Do not initialize it with a
   README, license or gitignore, because local history already contains these project files.
2. In this project directory replace YOUR_ACCOUNT with the actual account:

```sh
git remote add origin https://github.com/YOUR_ACCOUNT/Sectalix.git
git push -u origin main --tags
```

3. Check Actions: CI and Release artifacts. Do not describe the matrix as green
   until all jobs finish successfully. Download the release artifact before it expires.
4. Configure branch/ruleset protection, restrict updates/deletions of v* tags and
   require successful CI for changes. Do not force-move v1.0.0 after distribution.

## PyPI: separate publication decision

Sectalix is officially distributed under the MIT License (see `LICENSE` in the repository
root and `license = "MIT"` in `pyproject.toml`). Ensure the package name ownership is confirmed
on PyPI before first upload. Any post-tag metadata or release updates follow semantic versioning.

Trusted Publishing is the required release path. Configure the publisher below,
approve the protected `pypi` environment, and let the tag workflow publish the
verified artifact. Manual upload is only a fallback after an explicit release
review; never use arbitrary files from a working directory. For a manual dry run:

```sh
python -m pip install twine
python -m twine check --strict release-dist/*
python -m twine upload --repository testpypi release-dist/*
# Only after explicit publication review:
python -m twine upload release-dist/*
```

Supply credentials through Twine's secure prompting/keyring mechanism, never as
committed text. TestPyPI and PyPI are separate accounts/credentials. Uploaded PyPI
version filenames cannot simply be replaced; review before uploading.

The publishing job uses a protected `pypi` GitHub Environment and
`id-token: write`. Configure required reviewers before creating the tag. On PyPI,
register a Trusted Publisher for owner `berkeboranozdemir`, repository `Sectalix`,
workflow `.github/workflows/publish.yml`, and environment `pypi`.

Official references:

- [GitHub Python CI](https://docs.github.com/en/actions/tutorials/build-and-test-code/python)
- [Python Packaging publication guide](https://packaging.python.org/en/latest/guides/publishing-package-distribution-releases-using-github-actions-ci-cd-workflows/)
- [PyPI Trusted Publishing](https://docs.pypi.org/trusted-publishers/using-a-publisher/)
