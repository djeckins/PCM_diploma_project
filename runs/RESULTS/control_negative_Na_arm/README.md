# Non-conducting Na+ arms

Negative control: the same filter variants facing Na+, where the selective filter blocks conduction. The monitor must stay silent.

System `kcsa_na`, conditions E71A, T75AE71A; complete crossings recorded: **1**.

## No in_system folder

A negative control carries (nearly) no crossings, so no in-system model can be trained on it. The reading is the one below: a monitor trained on the conducting proteins must stay silent here.

## Predicted without ever seeing these units (`zero_shot/`)

- τ = 1000 ps, negative control: `E71A/1` READY 0.00 (mean P 0.008), `T75AE71A/1` READY 0.00 (mean P 0.012)
- τ = 2000 ps, negative control: `E71A/1` READY 0.00 (mean P 0.050), `T75AE71A/1` READY 0.00 (mean P 0.076)
- τ = 4000 ps, negative control: `E71A/1` READY 0.00 (mean P 0.103), `T75AE71A/1` READY 0.00 (mean P 0.153)

Full held-out unit reports with reality checks: `unit-report-kcsa_na-E71A-1/`, `unit-report-kcsa_na-T75AE71A-1/`.
