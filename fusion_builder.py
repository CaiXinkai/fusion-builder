#!/usr/bin/env python3
"""Build a fusion transcript, CDS and protein from a GTF and gene-region FASTA files."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass(frozen=True)
class Feature:
    start: int
    end: int
    exon_number: int


@dataclass
class Transcript:
    transcript_id: str
    gene_name: str
    chrom: str
    strand: str
    exons: list[Feature]
    cds: list[Feature]


@dataclass(frozen=True)
class GeneRegion:
    sequence: str
    genomic_low: int
    genomic_high: int
    strand: str

    def extract(self, start: int, end: int) -> str:
        if not (self.genomic_low <= start <= end <= self.genomic_high):
            raise ValueError(
                f"Interval {start}-{end} lies outside FASTA region "
                f"{self.genomic_low}-{self.genomic_high}"
            )
        if self.strand == "+":
            left = start - self.genomic_low
            right = end - self.genomic_low + 1
        else:
            left = self.genomic_high - end
            right = self.genomic_high - start + 1
        return self.sequence[left:right]


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        key, value = item.split(maxsplit=1)
        attributes[key] = value.strip('"')
    return attributes


def load_transcript(gtf_path: Path, transcript_id: str) -> Transcript:
    wanted = transcript_id.split(".")[0]
    transcript: Transcript | None = None
    with gtf_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] not in {"transcript", "exon", "CDS"}:
                continue
            attrs = parse_attributes(fields[8])
            current = attrs.get("transcript_id", "")
            if current.split(".")[0] != wanted:
                continue
            if transcript is None:
                transcript = Transcript(
                    transcript_id=current,
                    gene_name=attrs.get("gene_name", ""),
                    chrom=fields[0],
                    strand=fields[6],
                    exons=[],
                    cds=[],
                )
            if fields[2] in {"exon", "CDS"}:
                feature = Feature(
                    start=int(fields[3]),
                    end=int(fields[4]),
                    exon_number=int(attrs["exon_number"]),
                )
                (transcript.exons if fields[2] == "exon" else transcript.cds).append(feature)
    if transcript is None or not transcript.exons or not transcript.cds:
        raise ValueError(f"Transcript {transcript_id!r} with exon/CDS features not found")
    transcript.exons.sort(key=lambda x: x.exon_number)
    transcript.cds.sort(key=lambda x: x.exon_number)
    return transcript


def load_gene_region(path: Path) -> GeneRegion:
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().strip()
        sequence_lines: list[str] = []
        for line in handle:
            if line.startswith(">"):
                break
            sequence_lines.append(line.strip())
        sequence = "".join(sequence_lines).upper()
    match = re.match(r">[^:]+:(c?)(\d+)-(\d+)", header)
    if not match:
        raise ValueError(f"Cannot parse genomic interval from FASTA header: {header}")
    complement, first, second = match.groups()
    first_i, second_i = int(first), int(second)
    low, high = sorted((first_i, second_i))
    strand = "-" if complement == "c" or first_i > second_i else "+"
    expected = high - low + 1
    if len(sequence) != expected:
        raise ValueError(f"{path}: sequence length {len(sequence)} != interval length {expected}")
    return GeneRegion(sequence, low, high, strand)


def reconstruct_transcript(tx: Transcript, region: GeneRegion) -> tuple[str, dict[int, tuple[int, int]]]:
    if tx.strand != region.strand:
        raise ValueError(f"Strand mismatch for {tx.transcript_id}: GTF={tx.strand}, FASTA={region.strand}")
    pieces: list[str] = []
    bounds: dict[int, tuple[int, int]] = {}
    cursor = 0
    for exon in tx.exons:
        piece = region.extract(exon.start, exon.end)
        bounds[exon.exon_number] = (cursor, cursor + len(piece))
        pieces.append(piece)
        cursor += len(piece)
    return "".join(pieces), bounds


def genomic_base_to_transcript_offset(
    tx: Transcript, bounds: dict[int, tuple[int, int]], genomic_position: int
) -> int:
    for exon in tx.exons:
        if exon.start <= genomic_position <= exon.end:
            exon_start = bounds[exon.exon_number][0]
            within = (
                genomic_position - exon.start
                if tx.strand == "+"
                else exon.end - genomic_position
            )
            return exon_start + within
    raise ValueError(f"Genomic position {genomic_position} is not in an exon of {tx.transcript_id}")


def cds_start_offset(tx: Transcript, bounds: dict[int, tuple[int, int]]) -> int:
    first_cds = min(tx.cds, key=lambda x: x.exon_number)
    genomic_start = first_cds.start if tx.strand == "+" else first_cds.end
    return genomic_base_to_transcript_offset(tx, bounds, genomic_start)


def translate_to_stop(sequence: str) -> tuple[str, str, bool]:
    amino_acids: list[str] = []
    cds_end = len(sequence) - (len(sequence) % 3)
    found_stop = False
    for index in range(0, cds_end, 3):
        codon = sequence[index:index + 3]
        aa = CODON_TABLE.get(codon, "X")
        if aa == "*":
            cds_end = index + 3
            found_stop = True
            break
        amino_acids.append(aa)
    return sequence[:cds_end], "".join(amino_acids), found_stop


def wrap_fasta(sequence: str, width: int = 70) -> str:
    return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


def write_fasta(path: Path, header: str, sequence: str) -> None:
    path.write_text(f">{header}\n{wrap_fasta(sequence)}\n", encoding="ascii")


def build_fusion(
    gtf: Path,
    five_gene_fasta: Path,
    three_gene_fasta: Path,
    five_transcript_id: str,
    three_transcript_id: str,
    five_exon: int,
    three_exon: int,
    three_offset: int,
) -> dict[str, object]:
    five_tx = load_transcript(gtf, five_transcript_id)
    three_tx = load_transcript(gtf, three_transcript_id)
    five_seq, five_bounds = reconstruct_transcript(five_tx, load_gene_region(five_gene_fasta))
    three_seq, three_bounds = reconstruct_transcript(three_tx, load_gene_region(three_gene_fasta))

    if five_exon not in five_bounds or three_exon not in three_bounds:
        raise ValueError("Requested exon number is absent from the selected transcript")
    donor_end = five_bounds[five_exon][1]
    acceptor_exon_start = three_bounds[three_exon][0]
    acceptor_start = acceptor_exon_start + three_offset
    if acceptor_start < 0 or acceptor_start > len(three_seq):
        raise ValueError("Three-prime acceptor offset is outside the transcript sequence")

    if three_offset < 0:
        # Intronic sequence is absent from spliced RNA, so recover it directly from
        # the transcript-oriented gene-region FASTA.
        exon = next(x for x in three_tx.exons if x.exon_number == three_exon)
        if three_tx.strand == "+":
            intron = load_gene_region(three_gene_fasta).extract(exon.start + three_offset, exon.start - 1)
        else:
            intron = load_gene_region(three_gene_fasta).extract(exon.end + 1, exon.end - three_offset)
        three_suffix = intron + three_seq[acceptor_exon_start:]
    else:
        three_suffix = three_seq[acceptor_start:]

    fusion_transcript = five_seq[:donor_end] + three_suffix
    start = cds_start_offset(five_tx, five_bounds)
    fusion_cds, protein, found_stop = translate_to_stop(fusion_transcript[start:])
    return {
        "transcript": fusion_transcript,
        "cds": fusion_cds,
        "protein": protein,
        "found_stop": found_stop,
        "cds_start_1based": start + 1,
        "junction_1based": donor_end,
        "five_transcript": five_tx.transcript_id,
        "three_transcript": three_tx.transcript_id,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build fusion transcript/CDS/protein sequences from exon-relative breakpoints."
    )
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--five-gene-fasta", type=Path, required=True)
    parser.add_argument("--three-gene-fasta", type=Path, required=True)
    parser.add_argument("--five-transcript", required=True)
    parser.add_argument("--three-transcript", required=True)
    parser.add_argument("--five-exon", type=int, required=True)
    parser.add_argument("--three-exon", type=int, required=True)
    parser.add_argument(
        "--three-offset",
        type=int,
        default=0,
        help="0-based offset from the 3' exon start in transcript direction; negative retains upstream intron",
    )
    parser.add_argument("--name", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    result = build_fusion(
        args.gtf,
        args.five_gene_fasta,
        args.three_gene_fasta,
        args.five_transcript,
        args.three_transcript,
        args.five_exon,
        args.three_exon,
        args.three_offset,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    common = (
        f"{args.name} five={result['five_transcript']}:exon{args.five_exon} "
        f"three={result['three_transcript']}:exon{args.three_exon}{args.three_offset:+d}"
    )
    write_fasta(args.output_dir / f"{args.name}.transcript.fa", common, str(result["transcript"]))
    write_fasta(args.output_dir / f"{args.name}.cds.fa", common, str(result["cds"]))
    write_fasta(args.output_dir / f"{args.name}.protein.fa", common, str(result["protein"]))
    print(
        f"{args.name}: transcript={len(str(result['transcript']))} nt, "
        f"CDS={len(str(result['cds']))} nt, protein={len(str(result['protein']))} aa, "
        f"stop={'yes' if result['found_stop'] else 'no'}"
    )


if __name__ == "__main__":
    main()
