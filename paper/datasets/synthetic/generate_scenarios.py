"""Generate trees from the 6 designed scenarios displayed
in the paper.
Seed is blocked, i.e. the generated trees are deterministic."""
import hashlib
from pathlib import Path
from pyfroger.simulator import StemmaBenchmark


def make_seed(experiment_id: str) -> int:
    """Derive a stable 32-bit seed from an experiment ID string.
    Ensures reproduction.
    """
    digest = hashlib.md5(experiment_id.encode("utf-8")).hexdigest()
    return int(digest, 16) % (2 ** 32)


def sanitize_filename(filename: str) -> str:
    """Sanitize a filename to make it valid, i.e. replace
    '.' with '_'."""
    return filename.replace(".", "_")


# ---- GLOBAL CONFIGURATION
# Unless specified otherwise, each scenario share the following base
# configuration:
N_WITNESSES = 50
N_LOCI = 1000
N_STATES = 9  # 0 is included by default as character
BASE_ERROR_RATE = 0.05
BASE_LOSS_RATE = 0
BASE_HOT_LOCI_FRACTION = 0
BASE_HOT_LOCI_STATES = None
BASE_HOT_LOCI_FACTORS = 0
BASE_CONTAMINATION_RATE = 0
BASE_CONTAMINATION_STRENGTH = 0
N_GENERATIONS = 6
MAX_CHILDREN = 4
EXTANT_ARCHETYPE = True


# --- SCENARIO 1 / NON-ANOMALOUS TRADITIONS ---
# All witnesses conserved, no hot loci and no contamination.
# Error rate varies with values taken in {0.05, 0.1, 0.2, 0.3}.

Path("generated_trees/scenario_1").mkdir(parents=True, exist_ok=True)

for error_rate in [0.05, 0.1, 0.2, 0.3]:
    for ix in range(100):
        bench = StemmaBenchmark(
            n_witnesses=N_WITNESSES,
            n_loci=N_LOCI,
            error_rate=error_rate,
            n_states=N_STATES,
            n_generations=N_GENERATIONS,
            loss_rate=BASE_LOSS_RATE,
            max_children=MAX_CHILDREN,
            hot_loci_fraction=BASE_HOT_LOCI_FRACTION,
            hot_loci_factor=BASE_HOT_LOCI_FACTORS,
            hot_loci_states=BASE_HOT_LOCI_STATES,
            contamination_rate=BASE_CONTAMINATION_RATE,
            contamination_strength=BASE_CONTAMINATION_STRENGTH,
            extant_archetype=EXTANT_ARCHETYPE,
            seed=make_seed(f"scenario_1/tree_{error_rate}_{ix}")
        ).generate()
        bench.write_all(
            sanitize_filename(f"generated_trees/scenario_1/tree_{error_rate}_{ix}"))


# --- SCENARIO 2 / TRADITION WITH MANUSCRIPT LOSS
# Loss rate varies in 0.1, 0.2, 0.5

Path("generated_trees/scenario_2").mkdir(parents=True, exist_ok=True)
for loss_rate in [0.1, 0.2, 0.5]:
    for ix in range(100):
        bench = StemmaBenchmark(
            n_witnesses=N_WITNESSES,
            n_loci=N_LOCI,
            error_rate=BASE_ERROR_RATE,
            n_states=N_STATES,
            n_generations=N_GENERATIONS,
            loss_rate=loss_rate,
            max_children=MAX_CHILDREN,
            hot_loci_fraction=BASE_HOT_LOCI_FRACTION,
            hot_loci_factor=BASE_HOT_LOCI_FACTORS,
            hot_loci_states=BASE_HOT_LOCI_STATES,
            contamination_rate=BASE_CONTAMINATION_RATE,
            contamination_strength=BASE_CONTAMINATION_STRENGTH,
            extant_archetype=EXTANT_ARCHETYPE,
            seed=make_seed(f"scenario_2/tree_{loss_rate}_{ix}")
        ).generate()
        bench.write_all(
            sanitize_filename(f"generated_trees/scenario_2/tree_{loss_rate}_{ix}"))

# --- SCENARIO 3 / HIGH POLYGENESIS
# Change in fraction of hot_loci (values in 0.1, 0.2, 0.5) and a mutation
# probability factor of 6
Path("generated_trees/scenario_3").mkdir(parents=True, exist_ok=True)
for hot_loci_fraction in [0.1, 0.2, 0.5]:
    for ix in range(100):
        bench = StemmaBenchmark(
            n_witnesses=N_WITNESSES,
            n_loci=N_LOCI,
            error_rate=BASE_ERROR_RATE,
            n_states=N_STATES,
            n_generations=N_GENERATIONS,
            loss_rate=BASE_LOSS_RATE,
            max_children=MAX_CHILDREN,
            hot_loci_fraction=hot_loci_fraction,
            hot_loci_factor=6,
            hot_loci_states=2,
            contamination_rate=BASE_CONTAMINATION_RATE,
            contamination_strength=BASE_CONTAMINATION_STRENGTH,
            extant_archetype=EXTANT_ARCHETYPE,
            seed=make_seed(f"scenario_3/tree_{hot_loci_fraction}_{ix}")
        ).generate()
        bench.write_all(
            sanitize_filename(
                f"generated_trees/scenario_3/tree_{hot_loci_fraction}_{ix}"))

# --- SCENARIO 4 / HIGH CONTAMINATION
# Contamination rate in {0.1, 0.2, 0.3, 0.7} and a contamination strength of 5%.
Path("generated_trees/scenario_4").mkdir(parents=True, exist_ok=True)
for contamination_rate in [0.1, 0.2, 0.3, 0.7]:
    for ix in range(100):
        bench = StemmaBenchmark(
            n_witnesses=N_WITNESSES,
            n_loci=N_LOCI,
            error_rate=BASE_ERROR_RATE,
            n_states=N_STATES,
            n_generations=N_GENERATIONS,
            loss_rate=BASE_LOSS_RATE,
            max_children=MAX_CHILDREN,
            hot_loci_fraction=BASE_HOT_LOCI_FRACTION,
            hot_loci_factor=BASE_HOT_LOCI_FACTORS,
            hot_loci_states=BASE_HOT_LOCI_STATES,
            contamination_rate=contamination_rate,
            contamination_strength=0.05,
            extant_archetype=EXTANT_ARCHETYPE,
            seed=make_seed(f"scenario_4/tree_{contamination_rate}_{ix}")
        ).generate()
        bench.write_all(
            sanitize_filename(
                f"generated_trees/scenario_4/tree_{contamination_rate}_{ix}")
        )

# --- SCENARIO 5 / HIGH CONTAMINATION / HIGH POLYGENESIS RISK
# hot_loci_fraction in {0.1, 0.2, 0.5}
# contamination_rate in {0.1, 0.2, 0.3, 0.7}
Path("generated_trees/scenario_5").mkdir(parents=True, exist_ok=True)
for contamination_rate in [0.1, 0.2, 0.3, 0.7]:
    for hot_loci_fraction in [0.1, 0.2, 0.5]:
        for ix in range(100):
            bench = StemmaBenchmark(
                n_witnesses=N_WITNESSES,
                n_loci=N_LOCI,
                error_rate=BASE_ERROR_RATE,
                n_states=N_STATES,
                n_generations=N_GENERATIONS,
                loss_rate=BASE_LOSS_RATE,
                max_children=MAX_CHILDREN,
                hot_loci_fraction=hot_loci_fraction,
                hot_loci_factor=6,
                hot_loci_states=2,
                contamination_rate=contamination_rate,
                contamination_strength=0.05,
                extant_archetype=EXTANT_ARCHETYPE,
                seed=make_seed(
                    f"scenario_5/tree_{contamination_rate}_{hot_loci_fraction}_{ix}")
            ).generate()
            bench.write_all(
                sanitize_filename(
                    f"generated_trees/scenario_5/tree_{contamination_rate}_{hot_loci_fraction}_{ix}")
            )

# --- SCENARIO 6 / RESOLUTION OF COLLATION
# number of loci in 50, 100, 500
Path("generated_trees/scenario_6").mkdir(parents=True, exist_ok=True)
for n_loci in [50, 100, 500]:
    for ix in range(100):
        bench = StemmaBenchmark(
            n_witnesses=N_WITNESSES,
            n_loci=n_loci,
            error_rate=BASE_ERROR_RATE,
            n_states=N_STATES,
            n_generations=N_GENERATIONS,
            loss_rate=BASE_LOSS_RATE,
            max_children=MAX_CHILDREN,
            hot_loci_fraction=BASE_HOT_LOCI_FRACTION,
            hot_loci_factor=BASE_HOT_LOCI_FACTORS,
            hot_loci_states=BASE_HOT_LOCI_STATES,
            contamination_rate=BASE_CONTAMINATION_RATE,
            contamination_strength=BASE_CONTAMINATION_STRENGTH,
            extant_archetype=EXTANT_ARCHETYPE,
            seed=make_seed(f"scenario_6/tree_{n_loci}_{ix}")
        ).generate()
        bench.write_all(
            sanitize_filename(
                f"generated_trees/scenario_6/tree_{n_loci}_{ix}")
        )
