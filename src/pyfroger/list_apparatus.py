"""Python module for the transformation of input data
into an 'apparat en liste'.
"""
from dataclasses import dataclass
import csv
from collections import defaultdict


@dataclass
class VariantGroup:
    """Bears the list of witnesses that share a given erroneous reading
    (in the sense of different from the selected reference).

    See 'La critique des textes et son automatisation', p. 82, for the
    rationale.
    """
    locus: str
    reading: str
    witnesses: frozenset[str]
    reference: str = ""  # Reading from the reference manuscript

    @property
    def is_error(self) -> bool:
        """True if this reading differs from the reference."""
        return self.reading != self.reference

    @property
    def is_singular(self) -> bool:
        """True if only one witness bears this reading (individual error)."""
        return len(self.witnesses) == 1

    def __repr__(self) -> str:
        wits = ", ".join(sorted(self.witnesses))
        return f"VariantGroup('{self.locus}', reading: '{self.reading}', "\
               f"witnesses: {{{wits}}}, reference reading: '{self.reference}'"


class ListApparatus:
    """Implements the apparat 'en liste' (p. 82) by comparing every manuscript
    with the reference one and creating iteratively constrasting groups.
    """

    def __init__(self,
                 groups: list[VariantGroup],
                 ref_witness: str = "ref") -> None:
        """Create a list apparatus from a list of variant group."""
        # The different variant groups
        self.groups = groups
        # The available witnesses
        self.witnesses: set[str] = {w for g in groups for w in g.witnesses}
        # The reference witness
        self.reference = ref_witness

    def error_groups(self,
                     min_size: int = 1,
                     ) -> list[VariantGroup]:
        """Returns the list of variant groups that differs
        from the reference manuscript.

        Args:
            min_size: Minimum witness count (1 = include singular readings,
                2 = shared errors only).
        """
        return [
            g for g in self.groups
            if g.is_error
            and len(g.witnesses) >= min_size
        ]

    def group_by_frequency(self,
                           min_group_size: int = 1,
                           min_frequency: int = 1
                           ) -> dict[frozenset[str], dict]:
        """Returns the frequency of each group group across loci.

        This is the input for the Froger method: each entry maps a
        frozenset of witnesses to the number of loci where they share
        the same (erroneous, in the sense of against the reference) reading.

        Args:
            min_group_size: Minimum number of witnesses in a group.
            min_frequency: The minimum frequency of apparition in order to
                consider the group variant.

        Returns:
            Dict mapping frozenset[str] → {
                "frequency": int,
                "readings": list[str],
                "loci": list[str],
            }
        """
        frequency: dict[frozenset[str], dict] = {}
        for group in self.groups:
            if group.is_error:
                # Filter: minimum size
                if len(group.witnesses) < min_group_size:
                    continue

                key = group.witnesses
                if key not in frequency:
                    frequency[key] = {
                        "frequency": 1,
                        "readings": [group.reading],
                        "reference": [group.reference],
                        "loci": [group.locus],
                    }
                else:
                    frequency[key]["frequency"] += 1
                    frequency[key]["readings"].append(group.reading)
                    frequency[key]["loci"].append(group.locus)
                    frequency[key]["reference"].append(group.reference)

        # Filter out groups
        filtered_frequency = {}
        for group in frequency.keys():
            if frequency[group]["frequency"] >= min_frequency:
                filtered_frequency[group] = frequency[group]
        return filtered_frequency

    @classmethod
    def from_csv(cls,
                 path: str,
                 ref_witness: str = "ref",
                 delimiter: str = ",",
                 exclude: list[str] = []) -> "ListApparatus":
        """Creates an apparat en liste from a CSV collation file.
        Optionally include some locus by ID.

        The format should look like this:
                locus,    MS_A,      MS_B,      MS_C
                locus_1,  reading_a, reading_b, reading_c
            
        The manuscript taken as the reference can be specified in columns
        'ref_witness', which defaults to 'ref'.
        """
        with open(path, encoding="utf-8") as f:
            reader = csv.reader(f, delimiter=delimiter)
            header = [h.strip() for h in next(reader)]

            # ── Identify column roles ──
            ref_idx: int | None = None
            witness_cols: list[tuple[int, str]] = []

            for i, col in enumerate(header):
                low = col.lower()
                if i == 0:
                    continue
                elif low == ref_witness.lower():
                    ref_idx = i
                else:
                    witness_cols.append((i, col))

            # If the reference column was not found, break
            if not ref_idx:
                raise ValueError("Reference column not found, abort.")

            # ── Parse rows ──
            groups: list[VariantGroup] = []
            for row in reader:
                locus = row[0].strip()
                ref = row[ref_idx].strip() if ref_idx is not None else ""
                if locus not in exclude:
                    # Group witnesses by shared reading
                    reading_to_wits: dict[str, set[str]] = defaultdict(set)
                    for idx, wit_name in witness_cols:
                        reading = row[idx].strip().lower() if idx < len(row) else ""
                        if reading:
                            reading_to_wits[reading].add(wit_name)

                    for reading, wits in reading_to_wits.items():
                        groups.append(VariantGroup(
                            locus=locus,
                            reading=reading,
                            witnesses=frozenset(wits),
                            reference=ref,
                        ))

        return cls(groups=groups,
                   ref_witness=ref_witness)

    @classmethod
    def from_tei(cls, path: str) -> "ListApparatus":
        """Creates an apparat en liste from a TEI file."""
        raise NotImplementedError
