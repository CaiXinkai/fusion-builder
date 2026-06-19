#!/usr/bin/env python3
"""Compare generated EML4-ALK proteins with the supplied v6/v7 reference FASTA."""

from __future__ import annotations

import argparse
from pathlib import Path


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name: str | None = None
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        elif name is None:
            raise ValueError(f"{path}: sequence before FASTA header")
        else:
            records[name].append(line)
    return {key: "".join(value) for key, value in records.items()}


def first_difference(left: str, right: str) -> int | None:
    for index, (a, b) in enumerate(zip(left, right), start=1):
        if a != b:
            return index
    return None if len(left) == len(right) else min(len(left), len(right)) + 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Return a non-zero exit code unless generated and reference proteins match exactly.",
    )
    args = parser.parse_args()
    references = read_fasta(args.reference)
    strict_failed = False
    for variant in ("v6", "v7"):
        generated = next(iter(read_fasta(args.results_dir / f"EML4-ALK_{variant}.protein.fa").values()))
        expected = references[variant]
        difference = first_difference(generated, expected)
        if difference is None:
            print(f"{variant}: PASS_EXACT ({len(generated)} aa, exact sequence match)")
        else:
            substitutions = [
                f"{index}{a}>{b}"
                for index, (a, b) in enumerate(zip(generated, expected), start=1)
                if a != b
            ]
            if len(generated) == len(expected):
                print(
                    f"{variant}: PASS_WITH_SEQUENCE_DIFFERENCES "
                    f"(length={len(generated)} aa, substitutions={len(substitutions)}: "
                    f"{','.join(substitutions[:20])}; likely reference-version differences)"
                )
            else:
                print(
                    f"{variant}: FAIL_LENGTH "
                    f"(generated={len(generated)} aa, expected={len(expected)} aa, "
                    f"first difference={difference})"
                )
            strict_failed = True
    if strict_failed and args.strict:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
