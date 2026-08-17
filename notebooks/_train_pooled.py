"""Train-and-cache one pooled model: python _train_pooled.py <tau_ps> <variant>."""

import sys
import warnings

warnings.filterwarnings("ignore")
import predict_lib as pl

tau, variant = float(sys.argv[1]), sys.argv[2]
b = pl.load_pooled_model(tau_ps=tau, variant=variant)
q = b["oof_quality"]
lift = q.get("ratio_to_chance", float("nan"))
print(f"tau={tau:.0f} {variant}: threshold {b['threshold']:.3f}, "
      f"positives {b['train_positives']} ({b['base_rate']:.2%}), lift x{lift:.2f}")
