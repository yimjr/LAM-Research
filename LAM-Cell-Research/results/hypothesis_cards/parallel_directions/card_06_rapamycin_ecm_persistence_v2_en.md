# Hypothesis Card 06 (revised): ECM/protease programs retained after rapamycin

## Definition

Retention is not defined by a non-significant post-rapamycin result. For each gene:

```text
E_pre  = TSC2-loss_vehicle - WT_vehicle
E_post = TSC2-loss_rapamycin - WT_rapamycin
suppression_fraction = 1 - E_post/E_pre
```

## Current observation

GSE179044 contains effect-size-defined retained or enhanced ECM/protease-related genes. ELANE meets replicate-concordant partial-retention criteria in both hydrogel and plastic and is the first cross-environment protease candidate; MMP2 meets the criteria in hydrogel only and is not eligible in plastic.

## Current conclusion

This remains a high-value candidate mechanism: mTOR inhibition may control growth without fully eliminating matrix-related pathology. ELANE is only a perturbation-model retained-gene candidate; this is not patient-level sirolimus persistence or a confirmed protease-resistant mechanism.

## QC robustness

GSE179044 is a processed perturbation expression matrix and does not have a single-cell baseline/strict-QC filtering layer. Model-level robustness is therefore assessed by concordance across two biological replicates and suppression-fraction sensitivity. Human-LAM cross-validation must still pass the 140-baseline versus 85-strict candidate-pool check.

## Next step

Crosswalk retained genes against human LAMCORE, spatial source attribution, proteolytic balance and lesion/niche evidence.
