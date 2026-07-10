
from pyfroger.list_apparatus import ListApparatus
from pyfroger.group_method import UnorientedStemma


if __name__ == "__main__":
    print("------")
    apparatus = ListApparatus.from_csv(
        "demo.csv",
        ref_witness="MS_C",
        delimiter=",")

    stemma = UnorientedStemma(list_apparatus=apparatus,
                              min_frequency=1)
    stemma.build_tree(resolve_anomaly_threshold=1)
    print(stemma.display_tree())

    print("-----")
    apparatus = ListApparatus.from_csv(
        "demo.csv",
        ref_witness="MS_C",
        delimiter=",",
        exclude=["locus_2"])

    stemma = UnorientedStemma(list_apparatus=apparatus,
                              min_frequency=1)
    print(stemma.display_tree())
