# v1.0.0 release seal and distribution guide

The user authorized formal freezing. ACTIVE_PHASE.md records FROZEN status; the
computational source and tests are immutable for this release. The local annotated
tag v1.0.0 identifies the release commit. A Git tag is a reference, not a cryptographic
signature or a write-protection mechanism. Protect main and release tags on GitHub.

## Local verification and CI

Run `python -m pytest -W error`: the frozen suite contains 559 tests.
The CI workflow targets Ubuntu, Windows and macOS with Python 3.10–3.13 (12 jobs).
These hosted jobs can only be confirmed after pushing to GitHub; a local Windows
run does not establish that all matrix combinations pass.

The tag-triggered publish.yml first calls the same full matrix, then builds fresh
wheel/sdist into release-dist (not the tracked historical dist directory), checks
tag/version equality and runs strict Twine metadata validation. Verified archives
are retained as a GitHub Actions artifact for 30 days. There is intentionally no
automatic PyPI upload, credential or id-token permission.

The packaging tests invoke setuptools directly, so CI installs that backend and
wheel explicitly in addition to the requested editable test/plot extras. This is
build/test tooling, not a new application dependency.

## Connect a new GitHub repository

1. Create an empty repository in your GitHub account. Do not initialize it with a
   README, license or gitignore, because local history already contains these project files.
2. In this project directory replace YOUR_ACCOUNT with the actual account:

```sh
git remote add origin https://github.com/YOUR_ACCOUNT/ThinWallX.git
git push -u origin main --tags
```

3. Check Actions: CI and Release artifacts. Do not describe the matrix as green
   until all jobs finish successfully. Download the release artifact before it expires.
4. Configure branch/ruleset protection, restrict updates/deletions of v* tags and
   require successful CI for changes. Do not force-move v1.0.0 after distribution.

## PyPI: separate publication decision

ThinWallX is officially distributed under the MIT License (see `LICENSE` in the repository
root and `license = "MIT"` in `pyproject.toml`). Ensure the package name ownership is confirmed
on PyPI before first upload. Any post-tag metadata or release updates follow semantic versioning.

For manual publication, use the verified GitHub release artifact, not arbitrary
files from a working directory. Set up a PyPI account with 2FA and a scoped upload
credential; do not commit it. First validate and optionally exercise TestPyPI:

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

An alternative future configuration is PyPI Trusted Publishing with a protected
GitHub environment and a separate minimal upload job. It is not enabled here.

Official references:

- [GitHub Python CI](https://docs.github.com/en/actions/tutorials/build-and-test-code/python)
- [Python Packaging publication guide](https://packaging.python.org/en/latest/guides/publishing-package-distribution-releases-using-github-actions-ci-cd-workflows/)
- [PyPI Trusted Publishing](https://docs.pypi.org/trusted-publishers/using-a-publisher/)
