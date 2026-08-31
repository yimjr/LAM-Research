# Stage 16 — LAM candidate identity gate audit

LAM cells scored: 90,641
Normal/control cells used for calibration: 40,457
Independent positive references: 2,257
Negative references: 85,460

## Frozen evidence model

Identity anchors are PMEL/MLANA/MITF plus CORE2/CORE3 module evidence and available upstream CORE3 continuous scores. ACTA2/ESR1/VEGFD/CTSK are supportive evidence only. Competing lineage modules are continuous penalties; pericyte/VSMC is conditional on weak LAM identity and is not a standalone exclusion rule.

The formal 777-gene LAMCORE CSV was not available in the resolved input roots, so it was not fabricated or used. The report records this as an unavailable optional reference.

## Assignment counts

- LAM_core_candidate: 208
- LAM_boundary_candidate: 65,930
- non_LAM_like: 24,503
- competing-lineage exclusion: 20,820

Thresholds were calibrated without consensus_state and the existing states were joined only for the final diagnostic table. No Step 7–13 script and no scVI training was called.
