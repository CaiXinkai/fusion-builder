#!/usr/bin/env python3
"""Reproduce the supplied EML4-ALK v6 and v7 fusion sequences."""

from __future__ import annotations

import argparse
from pathlib import Path

from fusion_builder import build_fusion, write_fasta


VARIANTS = {
    "v6": {"five_exon": 13, "three_offset": -69, "reported_offset": -69},
    # ALK exon 20 has GTF CDS phase 2. The reported +12 coding-base offset
    # therefore corresponds to +14 from the complete spliced exon boundary.
    "v7": {"five_exon": 14, "three_offset": 14, "reported_offset": 12},
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for name, rule in VARIANTS.items():
        result = build_fusion(
            gtf=args.input_dir / "EML4_and_ALK.gtf",
            five_gene_fasta=args.input_dir / "EML4_gene.fna",
            three_gene_fasta=args.input_dir / "ALK_gene.fna",
            five_transcript_id="ENST00000318522.10",
            three_transcript_id="ENST00000389048.8",
            five_exon=rule["five_exon"],
            three_exon=20,
            three_offset=rule["three_offset"],
        )
        annotation = (
            f"EML4-ALK_{name} EML4:{result['five_transcript']}:exon{rule['five_exon']} "
            f"ALK:{result['three_transcript']}:exon20{rule['three_offset']:+d} "
            f"reported_coding_offset={rule['reported_offset']:+d}"
        )
        write_fasta(args.output_dir / f"EML4-ALK_{name}.transcript.fa", annotation, str(result["transcript"]))
        write_fasta(args.output_dir / f"EML4-ALK_{name}.cds.fa", annotation, str(result["cds"]))
        write_fasta(args.output_dir / f"EML4-ALK_{name}.protein.fa", annotation, str(result["protein"]))
        print(
            f"{name}: transcript={len(str(result['transcript']))} nt, "
            f"CDS={len(str(result['cds']))} nt, protein={len(str(result['protein']))} aa, "
            f"stop={'yes' if result['found_stop'] else 'no'}"
        )


if __name__ == "__main__":
    main()
