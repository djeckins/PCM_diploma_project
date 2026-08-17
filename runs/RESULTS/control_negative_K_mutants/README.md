# Non-conducting K+ mutants

Negative control: filter mutants that do not conduct K+. A monitor that has learned physics must stay silent here; READY fractions above zero would be false alarms.

System `kcsa`, conditions G77AE71A, T75AE71A; complete crossings recorded: **3**.

## No in_system folder

A negative control carries (nearly) no crossings, so no in-system model can be trained on it. The reading is the one below: a monitor trained on the conducting proteins must stay silent here.

## Predicted without ever seeing these units (`zero_shot/`)

- τ = 1000 ps, negative control: `G77AE71A/1` READY 0.00 (mean P 0.006), `T75AE71A/1` READY 0.00 (mean P 0.012)
- τ = 2000 ps, negative control: `G77AE71A/1` READY 0.00 (mean P 0.027), `T75AE71A/1` READY 0.00 (mean P 0.031)
- τ = 4000 ps, negative control: `G77AE71A/1` READY 0.00 (mean P 0.080), `T75AE71A/1` READY 0.00 (mean P 0.061)

Full held-out unit reports with reality checks: `unit-report-kcsa-G77AE71A-1/`, `unit-report-kcsa-T75AE71A-1/`.
