# Security policy

## Supported versions

Fixes go into the next release from `main`. Older tags are not patched.

## Reporting a vulnerability

Please report privately through GitHub's
[security advisories](https://github.com/vinecopulib/vinecopulib/security/advisories/new)
rather than a public issue.

vinecopulib is a numerical library: it parses model files (JSON and CBOR) and
otherwise operates on in-memory numeric data. The most plausible reports are
crashes or memory errors from malformed model files, or from data whose shape
violates a documented precondition. Both are worth reporting.

Please include the vinecopulib version, compiler and platform, and the smallest
input that reproduces the problem.
