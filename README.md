# immugene1
This repository is a dynamic tool for the analysis of TCR variation in CDR3 region.

### Summary and justification of methodological choices in the `divCHAO` function

| Step / choice | What the code does | Technical justification | Biological justification |
|--------------|-------------------|--------------------------|--------------------------|
| Use of `data_list` | Processes a list of data frames, one per sample | Enables systematic iteration across multiple samples without code duplication | Each sample represents an independent TCR repertoire |
| Internal function `calcular_diversidad_chao1()` | Computes the Chao1 diversity estimator from clonal abundances | Implements the Chao1 formula using singleton and doubleton counts | Chao1 estimates true repertoire richness by correcting for unobserved rare clonotypes |
| Use of `readCount` | Uses read counts as clonal abundance | Provides a quantitative measure of clonotype frequency | Read abundance reflects clonal expansion and repertoire structure |
| Use of `allVHitsWithScore` | Separates and filters sequences by TCR chain (TRA, TRB, TRG, TRD) and excludes IG | Ensures reliable V(D)J annotation based on alignment scores | TCR chain identity is defined by V gene usage |
| Exclusion of IG sequences | Removes clonotypes annotated as immunoglobulins | Prevents contamination with non-TCR rearrangements | IG rearrangements do not belong to the TCR repertoire |
| NA / non-NA filtering | Uses NA values after separation as a logical inclusion criterion | Simple and reproducible filtering strategy | Ensures that only clonotypes from the correct TCR chain are analyzed |
| Focus on the CDR3 region | Groups clonotypes by `nSeqCDR3` or `aaSeqCDR3` | CDR3 is the standard unit for clonotype definition | The CDR3 region determines TCR antigen specificity |
| Separation of NT vs AA | Computes diversity at nucleotide or amino acid level | Allows assessment of fine-scale diversity (NT) or convergent recombination (AA) | Different nucleotide sequences can encode the same CDR3 amino acid sequence |
| Grouping and counting | Summarizes the abundance of each CDR3 sequence | Prepares the data required for correct Chao1 computation | Each unique CDR3 corresponds to a distinct T cell clone |
| Per-sample Chao1 calculation | Computes Chao1 independently for each sample | Enables direct comparison between repertoires | Reflects true differences in clonal richness across individuals |
| Handling empty repertoires | Assigns a value of zero when no valid clonotypes are present | Prevents computational errors and preserves output structure | Represents the effective absence of detectable diversity |
| Structured output format | Returns a tidy data frame including a `Sample` column | Ensures compatibility with downstream statistical analyses and visualization | Enables comparison of diversity metrics across patients or experimental groups |

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
