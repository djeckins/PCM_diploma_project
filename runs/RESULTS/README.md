# Results by role

One folder per narrative role. The five protein roles carry the same skeleton (`in_system/`, `zero_shot/`, `benchmark/`, `README.md`); the three exceptions are listed below. Raw pipeline artifacts stay in `runs/<system>-stride<N>/` — every number here traces back to them.

- [Gramicidin A](training_gramicidin/README.md)
- [MthK](training_mthk/README.md)
- [KcsA E71A (K+)](training_kcsa_E71A/README.md)
- [G77A-E71A with Na+ (zero-shot)](test_zero_shot_G77A_Na/README.md)
- [Connexin-43 (zero-shot)](test_cx43/README.md)
- [Non-conducting K+ mutants](control_negative_K_mutants/README.md)
- [Non-conducting Na+ arms](control_negative_Na_arm/README.md)
- [Cross-protein experiments](cross_protein/README.md)

The two negative controls, `control_negative_K_mutants/` and `control_negative_Na_arm/`, have no `in_system/`: a non-conducting arm has no crossings to train on, so there is no own model to report. `cross_protein/` has none of the three, because it belongs to no single protein; under its own `README.md` it holds the readable summaries of the cross-protein runs (`monitor-generalisation-tau{1000,2000,4000}/`) and of the PLS-FMA self-check (`plsfma-selfcheck/`).
