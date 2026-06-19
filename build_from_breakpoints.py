#!/usr/bin/env python3
"""Build candidate fusion sequences from arbitrary hg38 genomic breakpoints."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path

from fusion_builder import Feature, Transcript, cds_start_offset, translate_to_stop, write_fasta


@dataclass(frozen=True)
class Breakpoint:
    chrom: str
    position: int
    strand: str


@dataclass
class AnnotatedTranscript(Transcript):
    gene_id: str = ""
    transcript_name: str = ""
    transcript_type: str = ""
    tags: set[str] = field(default_factory=set)

    @property
    def low(self) -> int:
        return min(exon.start for exon in self.exons)

    @property
    def high(self) -> int:
        return max(exon.end for exon in self.exons)


class IndexedFasta:
    """Minimal samtools-compatible .fai reader using 1-based closed intervals."""

    def __init__(self, fasta: Path):
        self.fasta = fasta
        self.index_path = Path(str(fasta) + ".fai")
        if not self.index_path.exists():
            raise FileNotFoundError(
                f"Missing FASTA index: {self.index_path}. Create it with: samtools faidx {fasta}"
            )
        self.index: dict[str, tuple[int, int, int, int]] = {}
        for line in self.index_path.read_text(encoding="utf-8").splitlines():
            name, length, offset, bases_per_line, bytes_per_line, *_ = line.split("\t")
            self.index[name] = (
                int(length), int(offset), int(bases_per_line), int(bytes_per_line)
            )

    def resolve_chrom(self, chrom: str) -> str:
        alternatives = [chrom]
        alternatives.append(chrom[3:] if chrom.startswith("chr") else f"chr{chrom}")
        for candidate in alternatives:
            if candidate in self.index:
                return candidate
        raise KeyError(f"Chromosome {chrom!r} is absent from {self.fasta}")

    def fetch(self, chrom: str, start: int, end: int) -> str:
        chrom = self.resolve_chrom(chrom)
        length, offset, line_bases, line_bytes = self.index[chrom]
        if not (1 <= start <= end <= length):
            raise ValueError(f"Invalid interval {chrom}:{start}-{end}; chromosome length={length}")
        zero_start = start - 1
        zero_end = end
        byte_start = offset + (zero_start // line_bases) * line_bytes + zero_start % line_bases
        byte_end = offset + (zero_end // line_bases) * line_bytes + zero_end % line_bases
        with self.fasta.open("rb") as handle:
            handle.seek(byte_start)
            raw = handle.read(byte_end - byte_start + line_bytes)
        sequence = re.sub(rb"\s+", b"", raw).decode("ascii").upper()
        return sequence[: end - start + 1]


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def oriented_fetch(genome: IndexedFasta, tx: Transcript, start: int, end: int) -> str:
    sequence = genome.fetch(tx.chrom, start, end)
    return sequence if tx.strand == "+" else reverse_complement(sequence)


def parse_breakpoint(value: str) -> Breakpoint:
    match = re.fullmatch(r"(.+):(\d+):([+-])", value.strip())
    if not match:
        raise argparse.ArgumentTypeError(
            f"Invalid breakpoint {value!r}; expected CHROM:1-based-position:+|-"
        )
    chrom, position, strand = match.groups()
    return Breakpoint(chrom, int(position), strand)


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def parse_attrs(text: str) -> tuple[dict[str, str], set[str]]:
    attrs: dict[str, str] = {}
    tags: set[str] = set()
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        key, value = item.split(maxsplit=1)
        value = value.strip('"')
        if key == "tag":
            tags.add(value)
        else:
            attrs[key] = value
    return attrs, tags


def chrom_equivalent(left: str, right: str) -> bool:
    return left.removeprefix("chr") == right.removeprefix("chr")


def load_candidate_transcripts(
    gtf: Path,
    breakpoint: Breakpoint,
    gene: str | None,
    transcript_id: str | None,
    ignore_strand: bool,
) -> list[AnnotatedTranscript]:
    records: dict[str, AnnotatedTranscript] = {}
    wanted_id = transcript_id.split(".")[0] if transcript_id else None
    with open_text(gtf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] not in {"transcript", "exon", "CDS"}:
                continue
            if not chrom_equivalent(fields[0], breakpoint.chrom):
                continue
            if not ignore_strand and fields[6] != breakpoint.strand:
                continue
            attrs, tags = parse_attrs(fields[8])
            current = attrs.get("transcript_id", "")
            if not current:
                continue
            if wanted_id and current.split(".")[0] != wanted_id:
                continue
            gene_name = attrs.get("gene_name", attrs.get("gene_id", ""))
            gene_id = attrs.get("gene_id", "")
            if gene and gene not in {gene_name, gene_id, gene_id.split(".")[0]}:
                continue
            tx = records.setdefault(
                current,
                AnnotatedTranscript(
                    transcript_id=current,
                    gene_name=gene_name,
                    chrom=fields[0],
                    strand=fields[6],
                    exons=[],
                    cds=[],
                    gene_id=gene_id,
                    transcript_name=attrs.get("transcript_name", ""),
                    transcript_type=attrs.get("transcript_type", attrs.get("transcript_biotype", "")),
                    tags=set(),
                ),
            )
            tx.tags.update(tags)
            if fields[2] in {"exon", "CDS"}:
                exon_number = attrs.get("exon_number")
                if exon_number is None:
                    continue
                feature = Feature(int(fields[3]), int(fields[4]), int(exon_number))
                (tx.exons if fields[2] == "exon" else tx.cds).append(feature)
    candidates = []
    for tx in records.values():
        if not tx.exons or not tx.cds:
            continue
        tx.exons.sort(key=lambda item: item.exon_number)
        tx.cds.sort(key=lambda item: item.exon_number)
        if tx.low <= breakpoint.position <= tx.high:
            candidates.append(tx)
    return candidates


def transcript_rank(tx: AnnotatedTranscript) -> tuple[int, int, int, str]:
    preferred = bool({"MANE_Select", "Ensembl_canonical"} & tx.tags)
    basic = "basic" in tx.tags
    coding = tx.transcript_type == "protein_coding"
    return (not preferred, not basic, not coding, tx.transcript_id)


def breakpoint_context(tx: Transcript, position: int) -> str:
    for exon in tx.exons:
        if exon.start <= position <= exon.end:
            offset = position - exon.start + 1 if tx.strand == "+" else exon.end - position + 1
            return f"exon{exon.exon_number}:{offset}/{exon.end - exon.start + 1}"
    ordered = tx.exons
    for left, right in zip(ordered, ordered[1:]):
        low = min(left.start, right.start)
        high = max(left.end, right.end)
        if low < position < high:
            return f"intron{left.exon_number}-{right.exon_number}"
    return "transcript_boundary"


def prefix_exact(tx: Transcript, genome: IndexedFasta, position: int) -> tuple[str, int]:
    pieces: list[str] = []
    junction = 0
    for index, exon in enumerate(tx.exons):
        if exon.start <= position <= exon.end:
            start, end = (exon.start, position) if tx.strand == "+" else (position, exon.end)
            piece = oriented_fetch(genome, tx, start, end)
            pieces.append(piece)
            junction += len(piece)
            return "".join(pieces), junction
        passed = exon.end < position if tx.strand == "+" else exon.start > position
        if passed:
            piece = oriented_fetch(genome, tx, exon.start, exon.end)
            pieces.append(piece)
            junction += len(piece)
            if index + 1 < len(tx.exons):
                next_exon = tx.exons[index + 1]
                in_intron = (
                    exon.end < position < next_exon.start
                    if tx.strand == "+"
                    else next_exon.end < position < exon.start
                )
                if in_intron:
                    start, end = (
                        (exon.end + 1, position)
                        if tx.strand == "+"
                        else (position, exon.start - 1)
                    )
                    piece = oriented_fetch(genome, tx, start, end)
                    pieces.append(piece)
                    junction += len(piece)
                    return "".join(pieces), junction
    raise ValueError(f"Cannot place 5' breakpoint {position} in {tx.transcript_id}")


def prefix_spliced(tx: Transcript, genome: IndexedFasta, position: int) -> tuple[str, int]:
    pieces: list[str] = []
    for exon in tx.exons:
        if exon.start <= position <= exon.end:
            start, end = (exon.start, position) if tx.strand == "+" else (position, exon.end)
            pieces.append(oriented_fetch(genome, tx, start, end))
            return "".join(pieces), sum(map(len, pieces))
        passed = exon.end < position if tx.strand == "+" else exon.start > position
        if passed:
            pieces.append(oriented_fetch(genome, tx, exon.start, exon.end))
        else:
            break
    if not pieces:
        raise ValueError(f"No annotated 5' exon before breakpoint in {tx.transcript_id}")
    return "".join(pieces), sum(map(len, pieces))


def suffix_exact(tx: Transcript, genome: IndexedFasta, position: int) -> str:
    pieces: list[str] = []
    for index, exon in enumerate(tx.exons):
        if exon.start <= position <= exon.end:
            start, end = (position, exon.end) if tx.strand == "+" else (exon.start, position)
            pieces.append(oriented_fetch(genome, tx, start, end))
            pieces.extend(
                oriented_fetch(genome, tx, later.start, later.end)
                for later in tx.exons[index + 1:]
            )
            return "".join(pieces)
        before = position < exon.start if tx.strand == "+" else position > exon.end
        if before:
            start, end = (
                (position, exon.start - 1)
                if tx.strand == "+"
                else (exon.end + 1, position)
            )
            pieces.append(oriented_fetch(genome, tx, start, end))
            pieces.extend(
                oriented_fetch(genome, tx, later.start, later.end)
                for later in tx.exons[index:]
            )
            return "".join(pieces)
    raise ValueError(f"Cannot place 3' breakpoint {position} in {tx.transcript_id}")


def suffix_spliced(tx: Transcript, genome: IndexedFasta, position: int) -> str:
    pieces: list[str] = []
    for index, exon in enumerate(tx.exons):
        if exon.start <= position <= exon.end:
            start, end = (position, exon.end) if tx.strand == "+" else (exon.start, position)
            pieces.append(oriented_fetch(genome, tx, start, end))
            pieces.extend(
                oriented_fetch(genome, tx, later.start, later.end)
                for later in tx.exons[index + 1:]
            )
            return "".join(pieces)
        before = position < exon.start if tx.strand == "+" else position > exon.end
        if before:
            return "".join(
                oriented_fetch(genome, tx, later.start, later.end)
                for later in tx.exons[index:]
            )
    raise ValueError(f"No annotated 3' exon after breakpoint in {tx.transcript_id}")


def full_spliced(tx: Transcript, genome: IndexedFasta) -> tuple[str, dict[int, tuple[int, int]]]:
    pieces: list[str] = []
    bounds: dict[int, tuple[int, int]] = {}
    cursor = 0
    for exon in tx.exons:
        piece = oriented_fetch(genome, tx, exon.start, exon.end)
        bounds[exon.exon_number] = (cursor, cursor + len(piece))
        cursor += len(piece)
        pieces.append(piece)
    return "".join(pieces), bounds


def build_candidate(
    five_tx: AnnotatedTranscript,
    three_tx: AnnotatedTranscript,
    left: Breakpoint,
    right: Breakpoint,
    genome: IndexedFasta,
    mode: str,
) -> dict[str, object]:
    if mode == "exact":
        prefix, junction = prefix_exact(five_tx, genome, left.position)
        suffix = suffix_exact(three_tx, genome, right.position)
    else:
        prefix, junction = prefix_spliced(five_tx, genome, left.position)
        suffix = suffix_spliced(three_tx, genome, right.position)
    transcript = prefix + suffix
    _, five_bounds = full_spliced(five_tx, genome)
    native_cds_start = cds_start_offset(five_tx, five_bounds)
    # Prefix sequence before the CDS start is identical in exact and splice-aware
    # modes unless the breakpoint lies upstream of the start codon.
    if len(prefix) <= native_cds_start:
        return {
            "transcript": transcript, "cds": "", "protein": "", "found_stop": False,
            "junction": junction, "status": "breakpoint_before_5prime_start_codon",
        }
    cds, protein, found_stop = translate_to_stop(transcript[native_cds_start:])
    return {
        "transcript": transcript, "cds": cds, "protein": protein,
        "found_stop": found_stop, "junction": junction, "status": "translated",
    }


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build fusion transcript/CDS/protein candidates from arbitrary genomic breakpoints."
    )
    parser.add_argument("--genome", type=Path, required=True, help="Reference FASTA with a .fai index")
    parser.add_argument("--gtf", type=Path, required=True, help="Matching GTF, optionally gzip-compressed")
    parser.add_argument("--left-breakpoint", type=parse_breakpoint, required=True)
    parser.add_argument("--right-breakpoint", type=parse_breakpoint, required=True)
    parser.add_argument("--left-gene")
    parser.add_argument("--right-gene")
    parser.add_argument("--left-transcript")
    parser.add_argument("--right-transcript")
    parser.add_argument("--mode", choices=["exact", "spliced", "both"], default="both")
    parser.add_argument("--max-transcripts-per-side", type=int, default=5)
    parser.add_argument("--ignore-breakpoint-strand", action="store_true")
    parser.add_argument("--name", default="fusion")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    genome = IndexedFasta(args.genome)
    left_txs = load_candidate_transcripts(
        args.gtf, args.left_breakpoint, args.left_gene,
        args.left_transcript, args.ignore_breakpoint_strand,
    )
    right_txs = load_candidate_transcripts(
        args.gtf, args.right_breakpoint, args.right_gene,
        args.right_transcript, args.ignore_breakpoint_strand,
    )
    left_txs = sorted(left_txs, key=transcript_rank)[:args.max_transcripts_per_side]
    right_txs = sorted(right_txs, key=transcript_rank)[:args.max_transcripts_per_side]
    if not left_txs:
        raise SystemExit("No compatible 5' transcripts found at the left breakpoint")
    if not right_txs:
        raise SystemExit("No compatible 3' transcripts found at the right breakpoint")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    modes = ["exact", "spliced"] if args.mode == "both" else [args.mode]
    rows: list[dict[str, object]] = []
    for five_tx in left_txs:
        for three_tx in right_txs:
            for mode in modes:
                candidate_id = safe_name(
                    f"{args.name}.{five_tx.transcript_id}.{three_tx.transcript_id}.{mode}"
                )
                try:
                    result = build_candidate(
                        five_tx, three_tx, args.left_breakpoint,
                        args.right_breakpoint, genome, mode,
                    )
                    header = (
                        f"{candidate_id} five={five_tx.gene_name}:{five_tx.transcript_id} "
                        f"three={three_tx.gene_name}:{three_tx.transcript_id} mode={mode} "
                        f"breakpoints={args.left_breakpoint.chrom}:{args.left_breakpoint.position}:"
                        f"{args.left_breakpoint.strand}|{args.right_breakpoint.chrom}:"
                        f"{args.right_breakpoint.position}:{args.right_breakpoint.strand}"
                    )
                    write_fasta(args.output_dir / f"{candidate_id}.transcript.fa", header, str(result["transcript"]))
                    if result["cds"]:
                        write_fasta(args.output_dir / f"{candidate_id}.cds.fa", header, str(result["cds"]))
                        write_fasta(args.output_dir / f"{candidate_id}.protein.fa", header, str(result["protein"]))
                    error = ""
                except ValueError as exc:
                    result = {"transcript": "", "cds": "", "protein": "", "junction": "", "found_stop": False, "status": "error"}
                    error = str(exc)
                rows.append({
                    "candidate_id": candidate_id,
                    "mode": mode,
                    "five_gene": five_tx.gene_name,
                    "five_transcript": five_tx.transcript_id,
                    "five_context": breakpoint_context(five_tx, args.left_breakpoint.position),
                    "three_gene": three_tx.gene_name,
                    "three_transcript": three_tx.transcript_id,
                    "three_context": breakpoint_context(three_tx, args.right_breakpoint.position),
                    "transcript_length": len(str(result["transcript"])),
                    "junction_after_base": result["junction"],
                    "cds_length": len(str(result["cds"])),
                    "protein_length": len(str(result["protein"])),
                    "stop_found": result["found_stop"],
                    "status": result["status"],
                    "error": error,
                })

    report = args.output_dir / f"{safe_name(args.name)}.candidates.tsv"
    with report.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(
        f"Wrote {len(rows)} candidates from {len(left_txs)} x {len(right_txs)} transcripts "
        f"to {args.output_dir}"
    )
    print(f"Candidate report: {report}")


if __name__ == "__main__":
    main()
