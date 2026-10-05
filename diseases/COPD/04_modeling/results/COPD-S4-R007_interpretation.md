# COPD DeepFootprinting interpretation (COPD-S4-R007)

This result summarizes 2 model run(s) over 4,560 TREDNet phase-I features per model. Feature SHAP values are reported as mean absolute and mean signed contributions across the explained positive sequences, with deterministic feature-index tie breaks.

## Run summary

| Model | Explained sequences | Input positive seqlets | Retained patterns | Pattern seqlets | Significant JASPAR matches | Distinct TF names |
|---|---:|---:|---:|---:|---:|---:|
| enhancer | 1,000 | 463 | 5 | 249 | 6 | 6 |
| silencer | 1,000 | 343 | 4 | 256 | 12 | 10 |

## Interpretation boundary

TF-MoDISco discovers recurring sequence patterns in model attributions. The reported JASPAR q-values quantify sequence similarity between those patterns and reference motifs. A match nominates a TF-family hypothesis; it does not establish occupancy, direction of regulation, or COPD-specific binding.

Feature SHAP summaries explain these trained model predictions relative to the selected background controls. They are not causal effects and should not be interpreted as differential activity between people with and without COPD.
