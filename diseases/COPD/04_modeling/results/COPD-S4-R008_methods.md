# COPD allele-specific motif analysis (COPD-S4-R008)

## Tested records and denominators

The analysis tested all 15,303 candidate records with validated paired 2,001-bp REF and ALT sequences. Results remain separated by enhancer and silencer model. Each model has three explicitly reported denominators: all sequence-scorable candidates, candidates passing that model's held-out 5% FPR prediction threshold, and candidates satisfying that model's R004 predicted-causal definition. Candidate identity is the stable `candidate_record_id`.

## Motif selection and scoring

Only model-specific TF-MoDISco matches to JASPAR 2024 motifs with q <= 0.05 were tested. JASPAR probability matrices were read from the same local MEME file used to generate the R007 motif report. A pseudocount of 1e-06 was added to every base probability and each PWM row was renormalized. Log2 odds were calculated against the MEME background frequencies (A=0.25, C=0.25, G=0.25, T=0.25).

For each allele and motif, every sequence window that overlaps any base of the complete allele span was evaluated on both strands. This full-span rule covers SNVs, multinucleotide variants, insertions, and deletions. The best log-odds score was min-max normalized to the PWM's theoretical range. A relative score >= 0.8 defines a motif-compatible sequence site. A disruption is REF >= 0.8 and ALT < 0.8; a creation is REF < 0.8 and ALT >= 0.8. Score changes that do not cross the threshold are retained as quantitative deltas but are not labeled creation or disruption.

## Interpretation boundary

The 0.8 threshold is a prespecified computational heuristic, not an empirical binding probability. TF-MoDISco-to-JASPAR similarity nominates a motif or TF-family hypothesis. Neither a motif match nor an allele-specific PWM threshold crossing proves that the named TF binds in lung tissue, that binding changes in COPD, or that the variant is causal. Experimental occupancy and perturbation assays are required.
