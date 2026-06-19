#!/usr/bin/env python3

import unittest

from build_from_breakpoints import (
    prefix_exact,
    prefix_spliced,
    reverse_complement,
    suffix_exact,
    suffix_spliced,
)
from fusion_builder import Feature, Transcript


class FakeGenome:
    def __init__(self, sequence: str):
        self.sequence = sequence

    def fetch(self, chrom: str, start: int, end: int) -> str:
        return self.sequence[start - 1:end]


class BreakpointSlicingTests(unittest.TestCase):
    def setUp(self):
        self.sequence = "ACGT" * 20
        self.genome = FakeGenome(self.sequence)

    def genomic(self, start, end):
        return self.sequence[start - 1:end]

    def test_plus_strand_intronic_breakpoint(self):
        tx = Transcript(
            "plus", "PLUS", "chr1", "+",
            [Feature(1, 5, 1), Feature(11, 15, 2)],
            [Feature(1, 5, 1), Feature(11, 15, 2)],
        )
        self.assertEqual(prefix_exact(tx, self.genome, 8)[0], self.genomic(1, 8))
        self.assertEqual(prefix_spliced(tx, self.genome, 8)[0], self.genomic(1, 5))
        self.assertEqual(
            suffix_exact(tx, self.genome, 8),
            self.genomic(8, 10) + self.genomic(11, 15),
        )
        self.assertEqual(suffix_spliced(tx, self.genome, 8), self.genomic(11, 15))

    def test_minus_strand_intronic_breakpoint(self):
        tx = Transcript(
            "minus", "MINUS", "chr1", "-",
            [Feature(21, 25, 1), Feature(11, 15, 2)],
            [Feature(21, 25, 1), Feature(11, 15, 2)],
        )
        expected_prefix = reverse_complement(self.genomic(21, 25)) + reverse_complement(self.genomic(18, 20))
        expected_suffix = reverse_complement(self.genomic(16, 18)) + reverse_complement(self.genomic(11, 15))
        self.assertEqual(prefix_exact(tx, self.genome, 18)[0], expected_prefix)
        self.assertEqual(
            prefix_spliced(tx, self.genome, 18)[0],
            reverse_complement(self.genomic(21, 25)),
        )
        self.assertEqual(suffix_exact(tx, self.genome, 18), expected_suffix)
        self.assertEqual(
            suffix_spliced(tx, self.genome, 18),
            reverse_complement(self.genomic(11, 15)),
        )

    def test_partial_exon_is_preserved(self):
        plus = Transcript(
            "plus", "PLUS", "chr1", "+",
            [Feature(1, 9, 1)], [Feature(1, 9, 1)],
        )
        minus = Transcript(
            "minus", "MINUS", "chr1", "-",
            [Feature(11, 19, 1)], [Feature(11, 19, 1)],
        )
        self.assertEqual(prefix_exact(plus, self.genome, 5)[0], self.genomic(1, 5))
        self.assertEqual(suffix_exact(plus, self.genome, 5), self.genomic(5, 9))
        self.assertEqual(
            prefix_exact(minus, self.genome, 15)[0],
            reverse_complement(self.genomic(15, 19)),
        )
        self.assertEqual(
            suffix_exact(minus, self.genome, 15),
            reverse_complement(self.genomic(11, 15)),
        )


if __name__ == "__main__":
    unittest.main()
