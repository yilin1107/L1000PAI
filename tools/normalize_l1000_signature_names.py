#!/usr/bin/env python3
"""Normalize GMT first-column naming for L1000PAI resources.

This script rewrites GMT term names into one of two canonical forms:

1. ``perturbagen up`` or ``perturbagen dn``
2. ``perturbagen(cellline) up`` or ``perturbagen(cellline) dn``

Structured L1000 collections contain a fixed metadata prefix before the
perturbagen name. For those files, the cell line is extracted from the second
underscore-delimited token and the perturbagen is extracted from the remaining
tokens after stripping trailing dose annotations such as ``10uM`` or
``100ng_per_ml``.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


STRUCTURED_GMT_FILES = {
    "l1000_aby.gmt",
    "l1000_cp.gmt",
    "l1000_lig.gmt",
    "l1000_oe.gmt",
    "l1000_shRNA.gmt",
    "l1000_siRNA.gmt",
    "l1000_xpr.gmt",
}

SPACE_DIRECTION_RE = re.compile(r"^(?P<core>.+?)\s+(?P<direction>up|down|dn)$", re.IGNORECASE)
HYPHEN_DIRECTION_RE = re.compile(r"^(?P<core>.+?)-(?P<direction>up|down|dn)$", re.IGNORECASE)
DOSE_TOKEN_RE = re.compile(
    r"^\d+(?:\.\d+)?(?:[numk]?M|mg|ug|ng|pg|uL|mL|L)$",
    re.IGNORECASE,
)


def normalize_direction(direction: str) -> str:
    direction = direction.lower()
    if direction == "down":
        return "dn"
    return direction


def split_direction(term: str) -> tuple[str, str]:
    """Split a term into the base name and a normalized direction token."""
    match = SPACE_DIRECTION_RE.match(term)
    if match:
        return match.group("core"), normalize_direction(match.group("direction"))

    match = HYPHEN_DIRECTION_RE.match(term)
    if match:
        return match.group("core"), normalize_direction(match.group("direction"))

    raise ValueError(f"Could not detect direction suffix in term: {term}")


def strip_trailing_dose(parts: list[str]) -> list[str]:
    """Drop trailing dose annotations while preserving perturbagen names."""
    if len(parts) >= 3 and DOSE_TOKEN_RE.match(parts[-3]) and [p.lower() for p in parts[-2:]] == ["per", "ml"]:
        return parts[:-3]
    if parts and DOSE_TOKEN_RE.match(parts[-1]):
        return parts[:-1]
    return parts


def normalize_structured_term(term: str) -> str:
    core, direction = split_direction(term)
    parts = core.split("_")
    if len(parts) < 5:
        raise ValueError(f"Unexpected structured term format: {term}")

    cellline = parts[1]
    perturbagen_parts = strip_trailing_dose(parts[4:])
    perturbagen = "_".join(perturbagen_parts).strip()
    if not perturbagen:
        raise ValueError(f"Could not extract perturbagen from term: {term}")

    return f"{perturbagen}({cellline}) {direction}"


def normalize_simple_term(term: str) -> str:
    core, direction = split_direction(term)
    return f"{core} {direction}"


def normalize_term(term: str, file_name: str) -> str:
    if file_name in STRUCTURED_GMT_FILES:
        return normalize_structured_term(term)
    return normalize_simple_term(term)


def normalize_gmt_file(path: Path, write: bool) -> tuple[int, int]:
    changed = 0
    total = 0
    output_lines: list[str] = []

    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line:
                output_lines.append(raw_line)
                continue

            fields = line.split("\t")
            original_term = fields[0]
            normalized_term = normalize_term(original_term, path.name)
            if normalized_term != original_term:
                changed += 1
            fields[0] = normalized_term
            output_lines.append("\t".join(fields) + "\n")
            total += 1

    if write and changed:
        with path.open("w", encoding="utf-8", newline="") as handle:
            handle.writelines(output_lines)

    return changed, total


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Directory containing the L1000 resource files")
    parser.add_argument(
        "--write",
        action="store_true",
        help="Rewrite GMT files in place. Without this flag the script only reports planned changes.",
    )
    args = parser.parse_args()

    directory = args.directory
    if not directory.exists():
        raise SystemExit(f"Directory not found: {directory}")

    gmt_files = sorted(p for p in directory.rglob("*.gmt") if p.is_file())
    if not gmt_files:
        raise SystemExit(f"No GMT files found under: {directory}")

    for path in gmt_files:
        changed, total = normalize_gmt_file(path, write=args.write)
        action = "rewrote" if args.write else "would rewrite"
        print(f"{action}\t{path.name}\tchanged={changed}\ttotal={total}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
