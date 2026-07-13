"""Run the benchmark for the Froger algorithm and output the corresponding
Newick trees."""
import os
import glob
from pyfroger.list_apparatus import ListApparatus
from pyfroger.group_method import UnorientedStemma

INPUT_FOLDER = "../../datasets/synthetic/generated_trees"
OUTPUT_FOLDER = "../experiment_outputs/froger/hyperparameters"


# testing of several hyperparameters
# Factorial design of the different
# Method for resolving contamination
# None, Two by two,
# probabilistic threshold: 0.1, 0.3, 0.49
# Min frequency of group apparition: All groups (1),
# Filter out (2, 5, 10)
# Position of reference w1, w8, w15

for resolve_method in ["two_by_two", "probabilistic_cutoff"]:
    for min_frequency in [1, 2, 5, 10]:
        for reference in ["w1", "w8", "w15"]:
            for apparatus_path in glob.glob(f"{INPUT_FOLDER}/*/*.csv"):
                scenario = os.path.basename(os.path.dirname(apparatus_path))
                if scenario in ["scenario_1", "scenario_2", "scenario_3",
                                "scenario_4",
                                "scenario_5", "scenario_6"]:
                    experiment_id = os.path.splitext(
                        os.path.basename(apparatus_path))[0]
                    apparatus = ListApparatus.from_csv(
                        apparatus_path,
                        ref_witness=reference,
                        delimiter=",")

                    stemma = UnorientedStemma(apparatus=apparatus,
                                              min_frequency=min_frequency)
                    if resolve_method == "probabilistic_cutoff":
                        hyperparam_cutoff = [0.1, 0.3, 0.49]
                    else:
                        hyperparam_cutoff = [0]

                    for cuttoff in hyperparam_cutoff:
                        stemma.build_tree(
                            anomaly_method=resolve_method,
                            anomaly_threshold=cuttoff)

                        # Export the result as a stemma or
                        # an extended newick if
                        # there is
                        # some reticulation detected
                        try:
                            with open(f"{OUTPUT_FOLDER}/{scenario}_{experiment_id}_{resolve_method}_{min_frequency}_{reference}_{str(cuttoff).replace(".", "")}.nwk", "w") as f:
                                f.write(stemma.to_newick())
                        except Exception:
                            print(f"Detected contamination, dumping enewick")
                            with open(f"{OUTPUT_FOLDER}/{scenario}_{experiment_id}_{resolve_method}_{min_frequency}_{reference}_{str(cuttoff).replace(".", "")}.enwk", "w") as f:
                                f.write(stemma.to_extended_newick())

# BASE TESTING FOR GENERATION
# for apparatus_path in glob.glob(f"{INPUT_FOLDER}/*/*.csv"):
#     scenario = os.path.basename(os.path.dirname(apparatus_path))
#     if scenario in ["scenario_1", "scenario_2", "scenario_3", "scenario_4",
#                     "scenario_5", "scenario_6"]:
#         experiment_id = os.path.splitext(os.path.basename(apparatus_path))[0]
#         apparatus = ListApparatus.from_csv(
#             apparatus_path,
#             ref_witness="w1",
#             delimiter=",")

#         stemma = UnorientedStemma(apparatus=apparatus,
#                                   min_frequency=2)
#         stemma.build_tree(
#             anomaly_method="two_by_two")

#         # Export the result as a stemma or an extended newick if there is
#         # some reticulation detected
#         try:
#             with open(f"{OUTPUT_FOLDER}/{scenario}_{experiment_id}.nwk", "w") as f:
#                 f.write(stemma.to_newick())
#         except Exception:
#             with open(f"{OUTPUT_FOLDER}/{scenario}_{experiment_id}.enwk", "w") as f:
#                 f.write(stemma.to_extended_newick())

