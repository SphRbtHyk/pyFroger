"""Apply the 'Méthode par les groupes' on the different available datasets.
"""
from pprint import pprint

from pyfroger.list_apparatus import ListApparatus
from pyfroger.group_method import UnorientedStemma

import pandas as pd

# Add locus ID to froger data
tmp_data_file = "tmp/notre_besoin.csv"
new_lines = ""
with open("../../datasets/artificial_scribes/notre_besoin_test.tsv") as f:
    for ix, line in enumerate(f.readlines()):
        if ix == 0:
            new_lines += "locus\t" + line
        else:
            new_lines += f"locus_{ix}\t" + line
with open(tmp_data_file, "w") as f:
    f.write("\t".join(new_lines.split("\t"))) #[:40*14-2]))

# considered_MS = ["p1", "p4", "p9", "p7", "p3", "p11", "p10", "p5", "p8", "p14", "p6", "p15", "p2", "p13", "p12", "p16"]

tmp_csv = pd.read_csv(tmp_data_file, sep="\t") #[["locus"] + considered_MS]
tmp_csv.to_csv(tmp_data_file, index=False, sep=";")

apparatus = ListApparatus.from_csv(
    "tmp/notre_besoin.csv",
    ref_witness="T1",
    delimiter=";")

stemma = UnorientedStemma(list_apparatus=apparatus,
                          min_frequency=1)

if __name__ == "__main__":
    print("------")
    # print(apparatus.error_groups())
    # print(stemma.children_of)
    stemma.build_tree(resolve_anomaly_threshold=1)
    # pprint(stemma.children_of)

    print(stemma.display_tree())
    # print(stemma.display_groups())
    # print(stemma.contaminated_witnesses)
    # print(stemma.display_variant_groups())

    # Reconstruct without offending
    