Changelog
=========

pysolnp 2022.3.08
-------------------
Fixed bug that would give incorrect results for completely unconstrained problems.
Also removed support for Python 2.7 and 3.5, and added support for Python 3.10.

pysolnp 2021.4.30
-------------------
Serious issue found in releases 2021.3.8, 2021.4.25 and 2021.4.26 that caused incorrect output.
This has been fixed in this release and previous releases have been deprecated.

Older Releases
-------------------

**pysolnp 2021.4.26** [Deprecated due to bug in output]

Fixed bug where the converged flag would only be set correctly when the debug was set to true.

**pysolnp 2021.4.25** [Deprecated due to bug in output]

No changes, re-release due to issue with source code build in previous version.

**pysolnp 2021.3.8** [Deprecated due to bug in output]

- Fixed issues where build would fail on Windows with newer versions of pip
- Added outputs:

  1.  converged : Boolean that indicates if the algorithm converged or not
  2.  hessian_matrix : A nested list of the last Hessian Matrix used by pysolnp

**pysolnp 2021.1.27** [Deprecated due to build issues]

Add wheel to Python 3.9 for all platforms and fix below issues:

- A change to pip means that meta-data version must match file-name version. This broke the source code build of pysolnp with recent versions of pip.

**pygosolnp 2021.1.24**

Initial release of the PYGOSOLNP library to PyPi.
Pure Python 3.6+ library so no precompiled binaries released.

**pysolnp 2020.4.11**

Initial release of the PYSOLNP library to PyPi.
Release includes precompiled wheels for Python 2.7, 3.5-3.8 (excluding Python 3.5 for Windows due to compilation issues).

**cppsolnp 2020.4.11**

Initial release of the C++ SOLNP library.
