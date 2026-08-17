# Hyperparameter flatness on 104 columns (cross-protein regime, gramicidin held out, tau = 1000 ps)

18 grid points (depth 2-4 x k_events 2/4/8 x colsample 0.5/0.7), each fitted
on the MthK+KcsA pool and scored two ways: the inner-rotation selection
criterion (what hyperparameter selection sees) and the held-out gramicidin
median ratio (what cross-protein transfer delivers). Full table in
results.json.

## Findings

1. The criterion is NOT flat on the enlarged descriptor set: inner criterion
   spans 2.33-2.61, and held-out transfer spans x1.32-x2.07.
2. The two rankings are misaligned: the best inner-criterion points
   (depth 4) give the worst gramicidin transfer (x1.32-1.51), while
   depth 2 + colsample 0.5 — mediocre by the inner criterion — gives the
   best transfer (x1.99-2.07, recovering the pre-change level).
3. Mechanism consistent with dilution: at 104 columns, colsample 0.7 lets
   deep trees co-select correlated absolute+relative column pairs that fit
   the training proteins but do not transfer.

## Decision

The frozen grids STAND. Re-tuning hyperparameters by held-out-protein
performance would select on the test set — exactly the leakage this project
is built to refuse. The misalignment itself is the result: within-pool
selection criteria do not optimize cross-protein transfer, which bounds how
much any pooled monitor can be tuned for unseen architectures without
spending held-out data. The gramicidin LOPO figures are therefore reported
under the frozen, selection-blind grids.
