"""Run the benchmark for the Froger algorithm and output the corresponding
Newick trees."""
import os
import glob
from pyfroger.list_apparatus import ListApparatus
from pyfroger.group_method import UnorientedStemma

INPUT_FOLDER = "../../datasets/synthetic/generated_trees"
OUTPUT_FOLDER = "../experiment_outputs/froger/"


for apparatus_path in glob.glob(f"{INPUT_FOLDER}/*/*.csv"):
    scenario = os.path.basename(os.path.dirname(apparatus_path))
    if scenario in ["scenario_1", "scenario_2", "scenario_3", "scenario_6"]:
        experiment_id = os.path.splitext(os.path.basename(apparatus_path))[0]
        print(experiment_id)
        apparatus = ListApparatus.from_csv(
            apparatus_path,
            ref_witness="w1",
            delimiter=",")

        stemma = UnorientedStemma(list_apparatus=apparatus,
                                  min_frequency=5)
        stemma.build_tree(resolve_anomaly_threshold=10,
                          anomaly_method="two_by_two")

        # Export the result as a stemma

        with open(f"{OUTPUT_FOLDER}/{scenario}_{experiment_id}.nwk", "w") as f:
            f.write(stemma.to_newick())
