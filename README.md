# Fusion builder prototype

This prototype reconstructs spliced transcripts from a GTF plus transcript-oriented
gene-region FASTA files, joins them at exon-relative breakpoints, and emits the
fusion transcript, CDS, and translated protein.

The 3-prime offset is measured from the first base of the acceptor exon in
transcript direction:

- `0`: begin at the first exon base.
- Positive: skip bases at the start of the exon.
- Negative: retain that many bases immediately upstream in the intron.

For the supplied EML4-ALK examples:

- v6 is `EML4 exon 13 -> ALK exon 20 offset -69` (`E13;ins69;A20`).
- v7 is reported as 12 bases after the coding portion of ALK exon 20. Because
  this GTF gives exon 20 a CDS phase of 2, the equivalent offset from the full
  spliced exon boundary is `+14`.

Run:

```powershell
python .\reproduce_eml4_alk.py `
  --input-dir E:\fusion_builder\EML4-ALK `
  --output-dir .\results

python .\validate_reference.py `
  --reference E:\fusion_builder\EML4-ALK\EML4-ALK_v6v7.fa `
  --results-dir .\results
```

The validator reports `PASS_WITH_SEQUENCE_DIFFERENCES` when protein lengths
match but a small number of amino acids differ because the reference protein
uses another sequence version. Add `--strict` only when an exact character-for-
character protein match is required.

The reproduction script fixes the MANE/Ensembl canonical transcripts:

- EML4: `ENST00000318522.10`
- ALK: `ENST00000389048.8`

`EML4-ALK_v6v7.fa` contains protein sequences despite its generic `.fa`
extension. Generated transcript, CDS, and protein files are written separately.

The supplied reference proteins appear to use a sequence version different from
the GRCh38/GTF gene sequences. The validator reports exact amino-acid differences
instead of silently changing genome-derived bases to force a match.

## Arbitrary genomic breakpoints

`build_from_breakpoints.py` accepts hg38-style 1-based genomic breakpoints such
as those reported by CTAT-LR-fusion:

```powershell
python .\build_from_breakpoints.py `
  --genome D:\reference\GRCh38.primary_assembly.genome.fa `
  --gtf D:\reference\gencode.v49.annotation.gtf `
  --left-breakpoint "chr17:81512136:-" `
  --right-breakpoint "chr7:5528545:-" `
  --left-gene LEFT_GENE `
  --right-gene RIGHT_GENE `
  --name LEFT_GENE--RIGHT_GENE `
  --mode both `
  --output-dir .\breakpoint_results
```

The FASTA must have a samtools-compatible index:

```powershell
samtools faidx D:\reference\GRCh38.primary_assembly.genome.fa
```

The left breakpoint is treated as the 5-prime fusion partner and the right
breakpoint as the 3-prime partner. The final `+` or `-` is used to select GTF
transcripts on that genomic strand.

Modes:

- `exact`: preserves partial exons and the sequence between an intronic
  breakpoint and the nearest transcript exon. This represents the literal
  genomic junction before assuming an RNA splice.
- `spliced`: preserves partial exons, but for an intronic breakpoint uses the
  nearest compatible annotated donor/acceptor exon boundary.
- `both`: emits both interpretations. This is the recommended default when only
  genomic breakpoints are known.

If transcript IDs are not supplied, up to five transcripts per partner are
ranked with MANE Select / Ensembl canonical transcripts first. Use
`--left-transcript` and `--right-transcript` to force a specific isoform.

The `*.candidates.tsv` report records breakpoint context, transcript choices,
sequence lengths, translation status, and whether a stop codon was found.
