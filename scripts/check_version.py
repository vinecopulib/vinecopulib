#!/usr/bin/env python3
"""Check that the version is stated consistently everywhere it appears.

The version lives in several hand-maintained places: ``CMakeLists.txt``,
``version.hpp``, the Doxyfiles, ``NEWS.md``, ``CITATION.cff`` and
``.zenodo.json``. This checks that they agree.

Run with no arguments to check internal consistency (the ``version`` CI job on
every pull request)::

    python3 scripts/check_version.py

Add ``--released`` to additionally require that the top ``NEWS.md`` heading
carries a real date rather than ``(unreleased)``, and ``--tag vX.Y.Z`` to
require that a git tag matches the declared version (the release workflow)::

    python3 scripts/check_version.py --released --tag "$GITHUB_REF_NAME"

Exits 0 if every check passes, 1 otherwise, printing one line per check.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# vinecopulib encodes the version as major * 100000 + minor * 100 + patch; see
# the comment block in include/vinecopulib/version.hpp.
MAJOR_SCALE = 100000
MINOR_SCALE = 100


class Failures(list):
    """Collected problems, so one run reports every mismatch at once."""

    def check(self, condition: bool, message: str) -> bool:
        status = "ok  " if condition else "FAIL"
        print(f"  [{status}] {message}")
        if not condition:
            self.append(message)
        return condition


def read(relative_path: str) -> str | None:
    path = REPO_ROOT / relative_path
    if not path.is_file():
        return None
    return path.read_text(encoding="utf-8")


def encode(major: int, minor: int, patch: int) -> int:
    return major * MAJOR_SCALE + minor * MINOR_SCALE + patch


def lib_version_string(major: int, minor: int, patch: int) -> str:
    """"x_y" when the patch level is 0, "x_y_z" otherwise."""
    if patch == 0:
        return f"{major}_{minor}"
    return f"{major}_{minor}_{patch}"


def parse_project_version(failures: Failures) -> tuple[int, int, int] | None:
    """The single source of truth: project(vinecopulib VERSION X.Y.Z)."""
    text = read("CMakeLists.txt")
    if text is None:
        failures.check(False, "CMakeLists.txt is missing")
        return None
    match = re.search(
        r"project\s*\(\s*vinecopulib\s+VERSION\s+(\d+)\.(\d+)\.(\d+)",
        text,
        re.IGNORECASE,
    )
    if match is None:
        failures.check(False, "CMakeLists.txt: no project(vinecopulib VERSION X.Y.Z)")
        return None
    version = tuple(int(g) for g in match.groups())
    print(f"  project version: {version[0]}.{version[1]}.{version[2]}")
    return version  # type: ignore[return-value]


def check_version_header(failures: Failures, major: int, minor: int, patch: int) -> None:
    text = read("include/vinecopulib/version.hpp")
    if text is None:
        failures.check(False, "include/vinecopulib/version.hpp is missing")
        return

    match = re.search(r"#define\s+VINECOPULIB_VERSION\s+(\S+)", text)
    if not failures.check(match is not None, "version.hpp defines VINECOPULIB_VERSION"):
        return
    literal = match.group(1)  # type: ignore[union-attr]

    # A leading zero would make this octal, breaking the arithmetic that
    # version.hpp documents.
    is_decimal = re.fullmatch(r"0|[1-9]\d*", literal) is not None
    failures.check(
        is_decimal,
        f"VINECOPULIB_VERSION is a plain decimal literal (found {literal!r}"
        + ("" if is_decimal else "; a leading zero makes it octal")
        + ")",
    )

    if is_decimal:
        expected = encode(major, minor, patch)
        failures.check(
            int(literal) == expected,
            f"VINECOPULIB_VERSION == {expected} (found {literal})",
        )

    match = re.search(r'#define\s+VINECOPULIB_LIB_VERSION\s+"([^"]*)"', text)
    if failures.check(
        match is not None, "version.hpp defines VINECOPULIB_LIB_VERSION"
    ):
        expected_lib = lib_version_string(major, minor, patch)
        found_lib = match.group(1)  # type: ignore[union-attr]
        failures.check(
            found_lib == expected_lib,
            f'VINECOPULIB_LIB_VERSION == "{expected_lib}" (found "{found_lib}")',
        )


def check_doxyfiles(failures: Failures, version: str) -> None:
    for relative_path in ("docs/Doxyfile", "docs/Doxyfile.in"):
        text = read(relative_path)
        if text is None:
            continue
        match = re.search(r"^\s*PROJECT_NUMBER\s*=\s*(.*)$", text, re.MULTILINE)
        if match is None:
            failures.check(False, f"{relative_path}: no PROJECT_NUMBER")
            continue
        found = match.group(1).strip()
        if found.startswith("@") and found.endswith("@"):
            # Substituted by CMake from PROJECT_VERSION.
            print(f"  [ok  ] {relative_path}: PROJECT_NUMBER templated ({found})")
            continue
        failures.check(
            found == version,
            f"{relative_path}: PROJECT_NUMBER == {version} (found {found or '<empty>'})",
        )


def parse_news_heading(failures: Failures) -> tuple[str, str] | None:
    """Return (version, date-or-unreleased) from the top NEWS.md heading."""
    text = read("NEWS.md")
    if text is None:
        failures.check(False, "NEWS.md is missing")
        return None
    for line in text.splitlines():
        if not line.startswith("## "):
            continue
        match = re.match(r"##\s+vinecopulib\s+(\S+)\s+\((.+)\)\s*$", line)
        if match is None:
            failures.check(
                False,
                "NEWS.md: first heading is not '## vinecopulib X.Y.Z (date)'"
                f" (found {line.strip()!r})",
            )
            return None
        return match.group(1), match.group(2).strip()
    failures.check(False, "NEWS.md: no '## vinecopulib ...' heading found")
    return None


def check_citation(failures: Failures, version: str) -> None:
    text = read("CITATION.cff")
    if text is None:
        return
    match = re.search(r"^version:\s*['\"]?([^'\"\s]+)", text, re.MULTILINE)
    if match is None:
        failures.check(False, "CITATION.cff: no version field")
        return
    failures.check(
        match.group(1) == version,
        f"CITATION.cff: version == {version} (found {match.group(1)})",
    )


def check_zenodo(failures: Failures, version: str) -> None:
    text = read(".zenodo.json")
    if text is None:
        return
    match = re.search(r'"version"\s*:\s*"([^"]+)"', text)
    if match is None:
        # Optional: Zenodo falls back to the release tag.
        return
    failures.check(
        match.group(1) == version,
        f".zenodo.json: version == {version} (found {match.group(1)})",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--released",
        action="store_true",
        help="also require the top NEWS.md heading to carry a real date",
    )
    parser.add_argument(
        "--tag",
        default=None,
        metavar="vX.Y.Z",
        help="also require this git tag to match the declared version",
    )
    args = parser.parse_args()

    failures = Failures()
    print("Checking version consistency:")

    parsed = parse_project_version(failures)
    if parsed is None:
        print("\nFAILED: cannot determine the project version.")
        return 1
    major, minor, patch = parsed
    version = f"{major}.{minor}.{patch}"

    check_version_header(failures, major, minor, patch)
    check_doxyfiles(failures, version)
    check_citation(failures, version)
    check_zenodo(failures, version)

    news = parse_news_heading(failures)
    if news is not None:
        news_version, news_date = news
        failures.check(
            news_version == version,
            f"NEWS.md: top heading version == {version} (found {news_version})",
        )
        if args.released:
            failures.check(
                news_date.lower() != "unreleased",
                f"NEWS.md: top heading is dated, not '(unreleased)' (found"
                f" '({news_date})')",
            )
        else:
            print(f"  news date: ({news_date})")

    if args.tag is not None:
        expected_tag = f"v{version}"
        failures.check(
            args.tag == expected_tag,
            f"git tag == {expected_tag} (found {args.tag})",
        )

    if failures:
        print(f"\nFAILED: {len(failures)} check(s) did not pass.")
        return 1
    print(f"\nOK: version {version} is stated consistently.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
