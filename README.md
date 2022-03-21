# CTCF-motif directionality controls CTCF-mediated chromatin interactions and correlates with topological domain structure
## R scripts and data files containing the analysis worklof and results

### R scripts:

#### topological_domains.R
contains the data analysis workflow used in the report. E.g creation of naively created loops and loops under constraints like motive directionality and length. Furthermore, a workflow to check if particular ctcf loop types cross topological domain borders.

#### within_loop_analysis.R
workflow to check for enrichment of particular tfs within ctcf mediated loops. ANOVA and Tukeyhsk

#### utility_functions.R
mainly a function like bedtools shuffle but for GRanges objects.

### R data:

#### ctcf_bs_motif_orientation.R
contains object containing the ctcf chip-seq bs but with motif and motif orentation.

#### simulated_loops.Rdata
object containing simulated loops.

#### simulated_loops_motif_matched.Rdata
object containing simulated ctcf loops but under constraints like loop length, motif and motif direcrtionality.
