# Stage 10 biology annotation

- States analyzed independently: 20
- Each eligible state uses patient × group pseudobulk and `~ patient_id + group`.
- Each state has its own DE FDR; `global_padj` is supplemental and does not replace state-level FDR.
- No all-state multi-class or pooled binary DE model was constructed.
- Pathway/regulon outputs are explicit `not_available` placeholders unless a local state-specific result exists.
