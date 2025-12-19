# immugene1
This repository is a dynamic tool for the analysis of TCR variation in CDR3 region.

### Summary and justification of methodological choices in the `fun_rich` function

| Step / choice | What the code does | Technical justification | Biological justification |
|--------------|-------------------|-------------------------|--------------------------|
| Use of `data_list` | Processes a list of data frames (one per sample) | Enables systematic iteration across multiple samples without code duplication | Each sample represents an independent TCR repertoire |
| Use of `allVHitsWithScore` | Separates and filters sequences according to assigned V gene (TRA, TRB, TRG, TRD, TRX or IG exclusion) | Provides the most reliable V(D)J annotation, including alignment scores | TCR chain identity is defined by the V gene usage |
| Exclusion of IG sequences | Filters out sequences annotated as immunoglobulins (IG) | Prevents contamination of the analysis with non-TCR clonotypes | IG rearrangements do not belong to the TCR repertoire |
| NA / non-NA filtering after `separate()` | Uses NA values as a logical criterion for inclusion/exclusion | Simple and reproducible strategy to classify receptor chains | Ensures that only clonotypes of the correct TCR chain are analyzed |
| Focus on CDR3 (`nSeqCDR3`, `aaSeqCDR3`) | Groups clonotypes by nucleotide or amino acid CDR3 sequence | CDR3 is the standard unit for clonotype definition | The CDR3 region determines antigen specificity of the TCR |
| Grouping with `group_by()` + `summarize(n())` | Counts occurrences of each CDR3 sequence | Identifies unique clonotypes and their frequencies | Each unique CDR3 corresponds to a distinct T cell clone |
| Separation of NT vs AA (`.nt` / `.aa`) | Computes richness at nucleotide or amino acid level | Allows assessment of fine-scale diversity (nt) or convergent recombination (aa) | Different nucleotide sequences can encode the same CDR3 amino acid sequence |
| Merging samples (`full_join`) | Combines clonotypes across samples | Facilitates direct comparison between samples | Enables inter-repertoire diversity analyses |
| Richness calculation (`colSums(!is.na())`) | Counts unique clonotypes per sample | Direct implementation of the richness definition | Richness reflects the number of distinct T cell clones |
| Output formatting with `Sample` column | Returns a tidy data frame ready for downstream analysis | Ensures compatibility with statistical models and visualization tools | Enables comparison of diversity metrics across patients or groups |
