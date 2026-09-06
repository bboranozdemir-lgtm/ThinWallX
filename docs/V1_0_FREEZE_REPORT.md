# ThinWallX v1.0.0 — formal freeze

Authorized by the project owner on 2026-09-06.

- STATUS: FROZEN — v1.0 Production-Quality Technical Release Complete.
- Final local command: `python -m pytest -W error`.
- Result: **559 passed in 44.67 seconds**, zero failures and warnings.
- Environment: Windows, Python 3.14.6, pytest 8.4.2.
- Pre/post content hashes confirm all protected source and test files are unchanged.
- ACTIVE_PHASE.md and ROADMAP.md record the release seal; v0.9 is described
  explicitly as preparation incorporated into v1.0, not an invented separate release.
- .gitignore excludes caches, build intermediates and credentials; dist and docs remain tracked.
- CI: 3 operating systems x 4 Python versions. Hosted execution is pending GitHub push.
- Tag workflow: full matrix, fresh build, strict metadata check, downloadable artifacts.
- No remote push, PyPI upload, new mechanics or next phase was performed.

The initial commit and annotated v1.0.0 tag provide local release identity; neither
is cryptographic signing or server-side protection. Configure remote protections
when connecting GitHub. Existing dist archives are retained; the release workflow
rebuilds from the tagged checkout into a separate directory.

See [release instructions](RELEASE_GUIDE.md) for GitHub/PyPI steps and the unresolved
license selection required before open-source distribution. No license was chosen
on the owner's behalf.
