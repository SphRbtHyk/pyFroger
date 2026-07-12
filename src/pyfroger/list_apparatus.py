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
    rationale of building the actual apparatus.
    """
    locus: str  # The position of the variant in the collation
    reading: str  # The value of the reading
    witnesses: frozenset[str]  # The set of witnesses bearing this reading
    reference: str = ""  # Reading from the reference manuscript

    @property
    def is_error(self) -> bool:
        """Check if the reading differs from the reference."""
        return self.reading != self.reference

    @property
    def is_singular(self) -> bool:
        """Check if only one witness shared this reading (individual error,
        specific to only one witness)."""
        return len(self.witnesses) == 1

    def __repr__(self) -> str:
        wits = ", ".join(sorted(self.witnesses))
        return f"VariantGroup('{self.locus}', reading: '{self.reading}', "\
            f"witnesses: {{{wits}}}, reference reading: '{self.reference}'"


class ListApparatus:
    """Implements the apparat 'en liste' (p. 82) by comparing the readings of
    every manuscript with the reference one and iteratively creating
    constrasting groups.
    """

    def __init__(self,
                 groups: list[VariantGroup],
                 ref_witness: str = "ref") -> None:
        """Create a list apparatus from a list of variant group.

        Args:
            groups (list[VariantGroup]): A list of the different variant
                groups computed on the witnesses.
            ref_witness (str): The ID of the witness used for reference.
        """
        # The different variant groups
        self.groups = groups
        # The available witnesses
        self.witnesses: set[str] = {w for g in self.groups for w in g.witnesses}
        # The reference witness
        self.reference = ref_witness

    def error_groups(self,
                     min_size: int = 1,
                     ) -> list[VariantGroup]:
        """Returns the list of variant groups that differs
        from the reference manuscript.

        Args:
            min_size: Minimum witness count (1 = include singular readings,
                >= 2 = shared errors only).
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
        """Returns the frequency of each group group across the
        different loci.

        This will be the input of the Froger's method, as a list of witnesses
        and the number of readings that support them.

        Args:
            min_group_size (int): Minimum number of witnesses in a group.
            min_frequency (int): The minimum frequency of apparition of a 
            given group in order to keep the variant group.

        Returns:
            Dict mapping frozenset[str] → {
                "frequency": int,
                "readings": list[str],
                "loci": list[str],
            }: A dictionary containing for every group its frequency,
            the list of its readings and their location in the collation.
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

        # Filter out groups that have a low frequency
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
                 exclude: list[str] = [],
                 lacuna_sign: str = "",
                 drop_lacunae: bool = False) -> "ListApparatus":
        """Creates an apparat en liste from a CSV collation file.
        Unwanted locus can be removed through their ID.

        The format should look like this:
                locus,    MS_A,      MS_B,      MS_C
                locus_1,  reading_a, reading_b, reading_c

        The manuscript taken as the reference for computing erroneous readings
        can be specified in columns 'ref_witness', which defaults to 'ref'.

        Args:
            path (str): The path to the CSV file to load.
            ref_witness (str): The ID of the witness to take as a reference.
            delimiter (str):  The delimited of the CV.
            exlude (list[str]): An optional list of loci (by ID) to exlude
                from the collation process.
            lacuna_sign (str): The sign indicating that there is a lacuna in
                the collation. Defaults to an empty string.
            drop_lacunae (bool): Whether or not a reading should be dropped if
                any of the variants present a lacuna.

        Returns:
            ListApparatus: A ListApparatus object.
        """
        with open(path, encoding="utf-8") as f:
            reader = csv.reader(f, delimiter=delimiter)
            header = [h.strip() for h in next(reader)]

            # Establish columns
            ref_idx = None
            witness_cols = []

            for i, col in enumerate(header):
                low = col.lower()
                if i == 0:
                    continue
                elif low == ref_witness.lower():
                    ref_idx = i
                else:
                    witness_cols.append((i, col))

            # If the reference column was not found, break
            if ref_idx is None:
                raise ValueError(
                    f"Reference MS {ref_witness} not found, abort.")

            def cell(row: list[str], idx: int) -> str:
                """Read a (possibly missing) cell, stripped and lowercased."""
                return row[idx].strip().lower() if idx < len(row) else ""

            # Recursively parse the rows
            groups = []
            ix = 0
            for row in reader:
                if not row:
                    continue
                # Check if there is any lacuna in the row
                is_lacuna = any(cell(row, idx) == lacuna_sign
                                for idx, _ in witness_cols)
                locus = row[0].strip()
                ref = cell(row, ref_idx)

                # If the reference is missing for this locus,
                # there can be no meaningful comparison and the
                # locus should be dropped altogether
                if ref == lacuna_sign:
                    continue
                if drop_lacunae:
                    filter_mask = locus not in exclude and not is_lacuna
                else:
                    filter_mask = locus not in exclude
                if filter_mask:
                    # Group witnesses by shared reading
                    reading_to_wits = defaultdict(set)
                    for idx, wit_name in witness_cols:
                        reading = cell(row, idx)
                        # If there is a reading and it is not missing
                        if reading and not reading == lacuna_sign:
                            reading_to_wits[reading].add(wit_name)

                    for reading, wits in reading_to_wits.items():
                        groups.append(VariantGroup(
                            locus=locus,
                            reading=reading,
                            witnesses=frozenset(wits),
                            reference=ref,
                        ))
                ix += 1

        return cls(groups=groups,
                   ref_witness=ref_witness)

