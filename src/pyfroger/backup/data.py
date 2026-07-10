"""Python module for the transformation of input data
into an 'apparat en liste'.
"""
from dataclasses import dataclass
import csv
from collections import defaultdict


@dataclass
class VariantGroup:
    """Bears the list of witnesses that contains the given reading.
    See 'La critique des textes et son automatisation', p. 82, for the
    rational.
    """
    locus: str  # Variant location
    reading: str  # Reading
    reading_type: str  # Categorization of reading
    genealogical: bool  # Whether or not the variant reading can be considered
    # meaningful for genealogy
    witnesses: frozenset[str]  # Witnesses bearing this lesson

    def __repr__(self) -> str:
        wits = ", ".join(sorted(self.witnesses))
        return f"VariantGroup('{self.locus}', '{self.reading}', {{{wits}}})"


class ListApparatus:
    """Implements the apparat 'en liste' (p. 82)"""

    def __init__(self,
                 groups: list[VariantGroup]) -> None:
        self.groups = groups
        self.witnesses: set[str] = {w for g in groups for w in g.witnesses}

    def __repr__(self) -> str:
        lines = [f"Apparat en liste of ({len(self.witnesses)} witnesses, "
                 f"{len(self.groups)} groups)"]
        loci: dict[str, list[VariantGroup]] = {}
        for g in self.groups:
            loci.setdefault(g.locus, []).append(g)
        for locus, groups in loci.items():
            lines.append(f"  {locus}:")
            for g in groups:
                wits = ", ".join(sorted(g.witnesses))
                lines.append(f"    \"{g.reading}\": {wits}")
        return "\n".join(lines)

    def group_by_frequency(self) -> dict[frozenset, int]:
        """Returns the frequency of the apparition of each group of
        manuscripts and the corresponding readings."""
        frequency = {}
        for group in self.groups:
            print("===")
            print(group)
            if group.witnesses not in frequency:
                frequency[group.witnesses] = {"frequency": 1,
                                              "readings": [group.reading]}
            else:
                frequency[group.witnesses]["frequency"] += 1
                frequency[group.witnesses]["readings"].append(group.reading)
        return frequency

    @classmethod
    def from_csv(cls,
                 path: str) -> "ListApparatus":
        """Creates an 'apparat en liste' from a CSV file presenting
        alignments.

        Two types of CSV are accepted:
            - A raw alignment matrix, with witness names as
            columns and the value of each reading for every locus displayed
            in column. This is ideal if you work from a collation file using
            only the same language.

            For example:
                | MS_A | MS_B | MS_C |
            locus_1 | reading11 | reading12 | reading13
            locus_2 | reading21 | reading22 | reading23

            - An already encoded matrix marking the value of a witness
            for a given locus, with witness names as
            columns and the value of each reading for every locus displayed
            in column. This is required for working with multi-lingual
            data.

            For example:
                        | MS_A | MS_B | MS_C |
                locus_1 | a | a | b |
                locus_2 | a | b | a |
        """
        with open(path, encoding="utf-8") as f:
            reader = csv.reader(f)
            witnesses = next(reader)[1:]

            groups: list[VariantGroup] = []
            for row in reader:
                locus = row[0]
                # Group witnesses by shared reading
                reading_to_witnesses: dict[str, set[str]] = defaultdict(set)
                for witness, reading in zip(witnesses, row[1:]):
                    reading = reading.strip()
                    if reading:  # empty cells are lacunae
                        reading_to_witnesses[reading].add(witness)

                for reading, wits in reading_to_witnesses.items():
                    groups.append(VariantGroup(
                        locus=locus,
                        reading=reading,
                        reading_type="",
                        genealogical=True,
                        witnesses=frozenset(wits),
                    ))

        return cls(groups)

    @classmethod
    def from_tei(cls, path: str) -> "ListApparatus":
        """Creates an 'apparat en liste' from a TEI file presenting
        alignments."""
