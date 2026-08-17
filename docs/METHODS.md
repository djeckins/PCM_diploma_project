# METHODS

The system answers two questions about every frame of an ion-channel MD
trajectory: whether the configuration is ready to conduct an ion within the
next τ picoseconds and, if not, what is holding it back. Both are properties of
the frame: a static criterion answers whether the channel conducts, this one
answers whether the frame is ready.

## Label

`y(t) = 1 ⟺ the completion time of at least one crossing lies in (t, t+τ]`,
evaluated per trajectory. Right censoring: a frame with an incomplete window is
unknown, not negative (`valid = window complete OR y=1`); the labels step prints
the number of frames dropped. The "entrance" anchor is a control ablation
against circularity; an event without an observed entrance is excluded from the
entrance labels.

## Crossings

Events are collected in-house by a finite-state machine over zones along the
axis; the zones are separated by planes derived from lipid head groups, and a
crossing counts only when the full ladder is traversed inside the confinement
cylinder (key `data.pore.crossing_cylinder_radius_A`, part of the problem
definition). A jump larger than half the box is a periodic-image wrap and is
not counted as a passage. For gramicidin a pre-existing annotation also exists:
the source is declared in the config (`events.source`), the events step prints
a comparison of the sources, and the `auto` mode refuses.

## Descriptors

Ten blocks, organized by failure routes: geometry (including the orientation of
the lining carbonyls, a rearrangement route that leaves the lumen unchanged),
hydration, water wire, occupancy, named sites, electrostatics, symmetry,
fluctuations, dynamics, ion delivery. The schema is computed from the config;
every column carries its block, units, estimator, missing value and a companion
indicator (see SCHEMA.md, generated). Windows look strictly backward;
sign-carrying columns are canonicalized to the conduction direction of the
condition exactly once.

The pore profile is a 3D inscribed sphere over all volume-occupying atoms
(channel and lipids; water and ions are excluded, so the measured quantity is
the lumen available to the mobile phase), with reachability from the axis, a
"fence", and a Lipschitz repair over the envelope of all search slices. Atomic
radii follow a fallback ladder (resname+name → name → element) with sources:
Bondi 1964 (doi:10.1021/j100785a001), H — Rowland & Taylor 1996
(doi:10.1021/jp953141+), Ca — Alvarez 2013 (doi:10.1039/c3dt50599e).

## Models and evaluation

Splits are made at the level of independent assemblies (lineage), so all frames
of one assembly stay on the same side of a split. All preprocessing consists of
steps of a single pipeline fitted inside each fold.
Missing values: boosting handles them natively; the linear arm uses
[values ++ indicators] with the fold median. The linear model is ridge logistic
regression (newton-cholesky, unpenalized intercept), C selected by log-loss,
the penalty reported as effective degrees of freedom. Boosting is xgboost: the
learning rate is fixed, the number of rounds comes from early stopping with
rotation over trajectories, the leaf constraint is parameterized by the event
count (k_events), and base_score is pinned explicitly to the fold base rate.
Where fewer than two training trajectories carry both classes, the inner
rotation falls back to contiguous time blocks inside each event-carrying
trajectory, separated from their training data by a two-sided embargo of one
feature window plus one horizon; the fallback is recorded per fold as
`blocked_inner_rotation`. An arm whose training part carries a single class is
refused.
Calibration is monotone by construction (sigmoid; isotonic only when the number
of positives suffices). The head model is declared in advance: trees.

Comparison arms: the published Rao-2019 criterion on the authors' vendored
tables (it predicts local dewetting; two deviations are declared and measured,
per-frame application and the entrance radius), a linear control on the same
structural inputs, and a clock arm with no structural columns.

Metrics: average precision plus its ratio to the base rate (pooled and
per-fold, with degenerate folds excluded by name); across horizons only the
ratio-to-chance is compared. The monitor's headline quantity is event-level:
the fraction of forewarned crossings under a hysteresis alarm, with a random
null comparison. ECE is not computed; its content is carried by the slope and
intercept of the calibration line.

## Mechanism diagnosis

Contributions of the already-trained model (TreeSHAP/coefficients), aggregated
to the axes of the mechanism vocabulary: dehydration, pathway rearrangement,
ion placement, wire rupture, electrostatics, column immobility, ion delivery.
The verdict is the axis with the largest contribution toward "not ready"; a
negligible contribution yields "unclear". Axis specificity is published as an
odds ratio before the winner is chosen.

## Cross-protein evaluation

The pool is three conducting proteins: gramicidin A, MthK and the K+-conducting
KcsA-E71A arm. Leave-one-protein-out holds each of them out in turn and scores
it with a monitor fitted on the other two; the reported number is the median,
over the held protein's trajectories, of within-trajectory average precision
divided by that trajectory's own base rate. Every rotation runs twice: on the
full column set and on a physics-only subset that drops the columns whose
per-system value ranges do not overlap, since a column that never overlaps can
act as a hidden system label (21 of 100 modelled columns are dropped for this
pool).

The model trained on all three is then applied to systems it never saw: the
five KcsA negative-control arms (K+: G77A·E71A, T75A·E71A; Na+: E71A,
G77A·E71A, T75A·E71A), of which only the Na+ G77A·E71A arm leaks and so is the
one case where ranking stays definable; and connexin-43, held out as an unseen
protein. A recorded ablation refits the pool with connexin-43 added and scores
both pools on identical held-out units; it is the basis for keeping connexin-43
outside the pool.

All of this is produced by `tools/monitor_generalisation.py` under
`PCM2_GEN_TAU` (horizon), `PCM2_POOL_ADD_CX43` (the ablation) and `PCM2_INNER`
(the rejected protein-inner-selection variant); see `REPRODUCTION.md`.
