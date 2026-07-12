# Heinrichi benchmark (MrBayes / IQ-TREE / FastME)

Straightforward companion to the simulated-scenario benchmark: one real
dataset, three tools, one score.

## Layout

    config.sh              paths, container images, analysis parameters
    scripts/prepare.py     NEXUS -> the three tool-specific inputs
    scripts/score.py       trees -> Roos average sign distance
    sbatch/00_prepare.sbatch
    sbatch/10_iqtree.sbatch
    sbatch/20_mrbayes.sbatch
    sbatch/30_fastme.sbatch
    sbatch/40_score.sbatch
    submit_all.sh          submits the chain with dependencies

## Setup

Put the two inputs in `data/`:

    data/heinrichi.nex     the 37 x 1208 matrix
    data/heinrichi.graph   the 65-node adjacency file (true stemma)

Edit `config.sh`:

  * `ENGINE` defaults to `singularity` (use `docker` on a workstation)
  * `IMG_*` — point at the images the simulated benchmark already uses; the
    defaults are guesses, not endorsements
  * add your site's `--partition` / `--account` to the sbatch headers

Then, **on a login node** (compute nodes usually have no outbound network, and
three jobs converting the same image at once will race in the cache):

    ./pull_images.sh

This builds `sif/{iqtree,mrbayes,fastme}.sif` and verifies that the expected
binary is actually on `PATH` inside each image — worth checking, because
biocontainers and staphb images do not agree on entrypoints.

Finally:

    ./submit_all.sh

`submit_all.sh` refuses to submit if the `.sif` files or the inputs are missing,
rather than letting the jobs fail one by one an hour later.

Containers run with `--cleanenv --no-home`: the host environment (conda,
`PYTHONPATH`) is not leaked into the container, and `$HOME` is not auto-mounted
over it. The job directory is bound to `/data` and is the working directory.

## Why prepare.py exists

The three programs cannot read the same matrix.

| tool     | needs                                    | consequence                       |
|----------|------------------------------------------|-----------------------------------|
| IQ-TREE  | multistate morphology, <= 32 states      | PHYLIP, symbols kept as-is        |
| MrBayes  | `datatype=standard`, **<= 10 states**    | must be recoded: 17 states won't load |
| FastME   | a distance matrix (cannot use multistate)| distances computed here           |

**MrBayes recoding.** Characters with <= 10 observed states are renumbered
0-9. Characters with more are binarised by default (one 0/1 character per
reading) — the usual stemmatological move. `--many-states drop` discards them
instead. The counts are reported in `work/prepare_report.json`; check them,
because the MrBayes matrix is genuinely not the same data as the IQ-TREE one
and any comparison has to say so.

**FastME distances.** Pairwise-complete Hamming: disagreements over loci where
both witnesses are extant, divided by that count. Missing data contributes
neither a match nor a mismatch. `prepare.py` warns about pairs resting on
fewer than 50 shared loci — in this tradition (Ca, Da, Cf, V) there will be
several, and their positions in the FastME tree are not to be trusted.

## Why score.py exists

Two settings move the number more than the choice of method does:

  * **Zero-length branches are collapsed** (`COLLAPSE_EPS`, default 1e-6),
    matching PAUP*'s `collapse=yes`. The metric counts unweighted edges, so an
    arbitrary bifurcation where the data support none destroys ties that the
    true stemma genuinely has. Setting `--collapse-eps -1` disables it; the
    difference between the two runs is worth reporting.
  * **Build on 37, score on 34.** Da, I and J are in the matrix but not in the
    true stemma, so they take part in the search and not in the scoring.
    `score.py` intersects the label sets automatically and names the excluded
    taxa in the `note` column.

Two baselines are always printed:

  * `true_stemma` — the stemma against itself. Must be 100.0. If it isn't,
    the graph file was misread and nothing else in the table means anything.
  * `star` — every triple ties, so this is the score of asserting nothing.
    It is the floor any method must beat, and it is higher than people expect.

## Reading the output

    results/scores.tsv

    method        skeleton  n_scored  n_leaves  note
    true_stemma   100.0     34        34        sanity check
    star          ~         34        34        floor: every triple ties
    iqtree        ~         34        37        not scored: Da,I,J
    mrbayes       ~         34        37        not scored: Da,I,J
    fastme        ~         34        37        not scored: Da,I,J

## Caveats worth putting in the paper

  * MrBayes convergence: the job prints the final ASDSF. If it is not below
    0.01, raise `MB_NGEN` and rerun. A consensus from an unconverged chain is
    not a result.
  * The MrBayes tree is built on binarised/recoded characters. It is not
    scoring the same matrix as IQ-TREE. Say so.
  * FastME is a distance method on a tradition with heavy, structured missing
    data. Its score should be read together with the thin-pair warnings.
  * These are single runs. For the simulated scenarios you presumably report
    a distribution; a single point estimate here is not comparable to it.
