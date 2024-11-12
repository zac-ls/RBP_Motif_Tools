# RNA Binding Protein Motif Tools (for rMATS)

RNA-binding proteins (RBPs) are proteins that bind to RNA molecules and are involved in many aspects of RNA processing, including alternative splicing. RBPs bind to RNA sequence motifs (about 3–8 nucleotides) via their RNA recognition motif (RRM). Many RRMs and their binding sites have been experimentally characterized, and RBP-specific consensus motif sequences have been probabilistically determined (see [ESEfinder](https://esefinder.ahc.umn.edu/tools/ESE2/) and [RBPmap](http://rbpmap.technion.ac.il)).

The tools presented here are R functions designed to provide motif scores based on nucleotide sequences and known RBP position probability and weighted matrices (PPM and PWM). **These tools are primarily intended for downstream analysis of skipped exon tables (SE.MATS.JC/JCEC) generated by [rMATS](https://rmats.sourceforge.io)**, a computational tool for detecting differential alternative splicing events in short-read RNA-seq data. **However, they can be applied to to any similar input data.** Below is an example of a data frame input containing the necessary columns for analysis:

| geneSymbol  | chr | strand | exonStart_0base | exonEnd   |
|-------|--------|-----|-----------------|-----------|
| SLK | 10     | +  | 81110914        | 81110955  |
| OPA1 | 3      | +  | 193626092       | 193626202 |
| LDHB | 12     | -  | 21657751        | 21657835  |

(Note, if using an rMATS table, consider filtering `FDR` and `IncLevelDifference`)

## Description of Tools

#### **`extract_sequence`** - *Gets Sequences from Genomic Coordinates*

This tool leverages the Bioconductor packages [BSgenome](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) and its associated genome (e.g. BSgenome.Hsapiens.UCSC.hg38) to extract nucleotide sequences based on genomic coordinates — exon start (exonStart_0base) to exon end (exonEnd). For entries where the strand is negative (-), representing a reverse orientation compared to the standard convention, the extracted sequence is reverse complemented. A new data frame will look like this:
  
| geneSymbol  | chr | strand | exonStart_0base | exonEnd   | nucleotide_sequence                                 |
|-------|--------|-----|-----------------|-----------|-----------------------------------------------------|
| SLK | 10      | +  | 81110914        | 81110955  | CCGCCGATGTGGAAGTGGCCAGATTCTGAGCCGCCTGACTAGA         |
| OPA1 | 3      | +  | 193626092       | 193626202 | GGTCTGCTTGGTGAGCTCATTCTCTTACAACAACAAATTCAAGAGCATGAAGAGGAAGCGCGCAGAGCCGCTGGCCAATATAGCACGAGCTATGCCCAACAGAAGCGCAAG           |
| LDHB | 12     | -  | 21657751        | 21657835  | AGCCGGCGCCGGAGGAGACGCACGCAGCTGACTTTGTCTTCTCCGCACGACTGTTACAGAGGTCTCCAGAGCCTTCTCTCTCCTG             |

#### **`extract_sequence_with_flanks`** - *Gets Sequences and thier Flanking Regions from Genomic Coordinates*

This function is similar to `extract_sequence` except that it extracts sequences defined by both exonStart_0base to exonEnd and a desired number of flanking seqeunces (such as 50 +/- the exon region). For easy visualizaiton and downstream analysis, the flanking (intronic) regions are in lowercase.

| geneSymbol  | chr | strand | exonStart_0base | exonEnd   | nucleotide_sequence                                 
|-------|--------|-----|-----------------|-----------|--------------------------------------------------------------------|
| SLK | 10      | +  | 81110914        | 81110955  | agtcctcagaccccatgctgcctccaactgagccttgtgtttccttgcagCGCCGATGTGGAAGTGGCCAGATTCTGAGCCGCCTGACTAGAgttagtaagttgcctggcgttctcgtgcagtcactggcctctccagtggt   |
| OPA1 | 3      | +  | 193626092       | 193626202 | attattctcctccccaatttcctcttctcctcattgtgaactcgtggcagGGTCTGCTTGGTGAGCTCATTCTCTTACAACAACAAATTCAAGAGCATGAAGAGGAAGCGCGCAGAGCCGCTGGCCAATATAGCACGAGCTATGCCCAACAGAAGCGCAAGgtgatggatggtttaagggggctaccgatacattcacactaatcagccat                                                                                                                  |
| LDHB | 12      | -  | 21657751        | 21657835  | taagaggctgcggtggttgtggggccccgccccctcctccctccttgcagAGCCGGCGCCGGAGGAGACGCACGCAGCTGACTTTGTCTTCTCCGCACGACTGTTACAGAGGTCTCCAGAGCCTTCTCTCTCCTGgtaggtttcggctcaggaccctgaatcctggcccacaggcaagcctgatg                                                                                                                                            |

### **Motif Scoring & Plots** 
#### **`find_rbp_binding_and_generate_table`** - *Finds RBP Binding Sites and Scores from PPM or PWM Input; Plots and Tabulates Motif Positions and Scores*

This is a wrapper function that contains the following:

- **`handle_ppm_or_pwm`** - *checks if the input RBP matrix is a PWM or a PPM to handle appropriately*
- **`apply_pseudocount`** - *adds a pseudocount to an input PPM to avoid zero probabilities*
- **`generate_pwm_from_ppm`** - *calculates a PWM from an input PPM*
- **`calculate_max_score`** - *finds the maximum possible score for a PWM for thresholding criteria*
- **`score_motif_with_positions_pwm`** - *scores motifs using PWM or PPM input and returns the positions of high scores and their ranges*

Here, we can generate RBP sequence motif scores based on either the position probability matrix (PPM) or the position weighted matrix (PWM) of a characterized RBP. Such matrix information can be found and downloaded from certain online databases, such as [RBPmap](http://rbpmap.technion.ac.il/download.html). The below table represents a PPM of the RBP, *RBM47*. 

|Position   |1       |2 |3 |4 |5       |6       |
|---|----------|----|----|----|----------|----------|
| A | 0.2849425 | 1  | 1  | 0  | 0.13246623 | 0.63661831 |
| C | 0.1629815 | 0  | 0  | 0  | 0.38939470 | 0.08714357 |
| G | 0.1628814 | 0  | 0  | 0  | 0.03091546 | 0.08704352 |
| T | 0.3891946 | 0  | 0  | 1  | 0.44722361 | 0.18919460 |

PPMs represent the probabilities for each nucleotide at each position. These probabilities are derived from a matrix of raw frequencies (or counts) of each nucleotide letter at each position, known as a position count matrix (PCM). The PPM normalizes such counts (note that all positions sum to 1 in the table). However, for downstream handling of probability values of zero, and to account for possible biological variation in RBP binding, we can apply a pseudocount of $p = 0.01$. Here is the formula used by `apply_pseudocount` function: 

### $\ P_{\text{pseudo}}(N) = \frac{C_N + \frac{p}{n}}{\sum C + p} \$

Where:
- $\( P_{\text{pseudo}}(N) \)$ represents the adjusted probability at position $\( N \)$ after applying a pseudocount.
- $\( C_N \)$ is the observed count for a specific nucleotide at position $\( N \)$ in the original Position Probability Matrix (PPM).
- $\( p \)$ is the pseudocount value added to avoid zero probability values.
- $\( n \)$ is the total number of nucleotide types (for RNA, $\( n = 4 \)$ for A, C, G, and U (in this case as T)).
- $\( \sum C \)$ is the sum of observed counts across all nucleotides at position $\( N \)$.

Here is the pseudocount-adjusted PPM for *RBM47*: 

|Positon   |1       |2       |3       |4       |5       |6       |
|---|----------|----------|----------|----------|----------|----------|
| A | 0.284597 | 0.992574 | 0.992574 | 0.002475 | 0.133630 | 0.632790 |
| C | 0.163843 | 0.002475 | 0.002475 | 0.002475 | 0.388015 | 0.088756 |
| G | 0.163744 | 0.002475 | 0.002475 | 0.002475 | 0.033085 | 0.088657 |
| T | 0.387816 | 0.002475 | 0.002475 | 0.992574 | 0.445271 | 0.189797 |

Next, to correct for background nucleotide frequencies, position weight matrices (PWMs) can be used. A PWM represents the "expected" motif sequence scores by incorporating background nucleotide probabilities, assuming uniform distribution (i.e., $1/4$ for each nucleotide). To compute a PWM from a PPM, the `generate_pwm_from_ppm` function applies the following log-odds scoring formula:

### $S(N) = \log_2 \left( \frac{P(C_N)}{B_N} \right)$

Where:
  - $\( S(N) \)$ is score at position $\( N \)$.
  - $\( P(C_N) \)$ is the probability of the observed counts $\( C_N \)$, which corresponds to values in the Position Probability Matrix (PPM).
  - $\( B_N \)$ is the background nucleotide frequency at position $\( N \)$, representing the expected probability of each nucleotide by chance.

Here is the PWM table generated from the above pseudocount-adjusted PPM of *RBM47*:

| Position  | 1      | 2       | 3       | 4       | 5       | 6      |
|---|---------|----------|----------|----------|----------|---------|
| A | 0.187   | 1.989    | 1.989    | -6.658   | -0.904   | 1.340   |
| C | -0.610  | -6.658   | -6.658   | -6.658   | 0.634    | -1.494  |
| G | -0.610  | -6.658   | -6.658   | -6.658   | -2.918   | -1.496  |
| T | 0.633   | -6.658   | -6.658   | 1.989    | 0.833    | -0.397  |

For a more detailed explanation of the mathematical principles behind sequence motif matrices, refer to [this vignette](https://bioconductor.org/packages/devel/bioc/vignettes/universalmotif/inst/doc/IntroductionToSequenceMotifs.pdf).

The function `score_motif_with_positions_pwm` slides a window of (i.e. scans) the motif's length across the sequence one position at a time and generates a score based on the matching nucleotides of the PWM at each window. Of course, we shouldn't accept all the motifs that are identified, so a threshold is applied. The thresholding method used here discards scores less than or equal to a percentage of the maximal possible score for a PWM (the max score is obtained by the `calculate_max_score` function).

Finally, the `find_rbp_binding_and_generate_table` function is the wrapper function that outputs results, including plots of motif scores (within a user-defined threshold) of each RBP along each sequence, with reference to the position of the exon (see below). These plots are exported as PDFs to subfolder in the working directory.

![image](https://github.com/user-attachments/assets/2dcd8d41-1376-405f-9cb9-5bc3a5acd4d6)

This output shows motif scores at a threshold of 50% maximal score for each RBP's PWM in a sequence from OPA1 gene.

The function also generates a table of the RBP motif analyses, including the gene symbol, RBP, motif sequence, the motif's score, its relative start and end positions within the sequence, and it's genomic start and end coordinates. 

| ID | gene | RBP   | motif_sequence | score     | start | end | genomic_start | genomic_end |
|----|------|-------|----------------|-----------|-------|-----|---------------|-------------|
| 1  | SLK  | RBM47 | AAATGA         | 4.576834  | 56    | 61  | 104010871     | 104010876   |
| 1  | SLK  | SRSF2 | TGCTACTA       | 3.240000  | 162   | 169 | 104010977     | 104010984   |
| 2  | OPA1 | FUS   | AGCGCGC        | 10.032832 | 107   | 113 | 193626198     | 193626204   |
| 2  | OPA1 | RBM47 | CAATTT         | 5.793409  | 15    | 20  | 193626106     | 193626111   |
| 2  | OPA1 | RBM47 | AAATTC         | 5.493475  | 85    | 90  | 193626176     | 193626181   |
| 2  | OPA1 | RBM47 | TAATCA         | 8.575169  | 201   | 206 | 193626292     | 193626297   |
| 2  | OPA1 | SRSF2 | AACTCCTG       | 4.950000  | 39    | 46  | 193626130     | 193626137   |
| 2  | OPA1 | SRSF2 | AGCCGCTG       | 5.560000  | 116   | 123 | 193626207     | 193626214   |
| 2  | OPA1 | SRSF2 | GGCCAATA       | 4.250000  | 123   | 130 | 193626214     | 193626221   |
| 2  | OPA1 | SRSF2 | GGCTACCG       | 5.100000  | 181   | 188 | 193626272     | 193626279   |
| 3  | LDHB | RBM47 | GAATCC         | 4.497426  | 158   | 163 | 21657908      | 21657913    |
| 3  | LDHB | SRSF2 | GGCCCCGG       | 3.490000  | 22    | 29  | 21657772      | 21657779    |
| 3  | LDHB | SRSF2 | AGCCGCGC       | 4.010000  | 51    | 58  | 21657801      | 21657808    |
| 3  | LDHB | SRSF2 | GTCTCCCA       | 4.020000  | 111   | 118 | 21657861      | 21657868    |
| 3  | LDHB | SRSF2 | GTCTCCAG       | 4.540000  | 112   | 119 | 21657862      | 21657869    |
| 3  | LDHB | SRSF2 | CTCTCCTG       | 3.530000  | 128   | 135 | 21657878      | 21657885    |
| 3  | LDHB | SRSF2 | GGACCCTG       | 5.430000  | 151   | 158 | 21657901      | 21657908    |
| 3  | LDHB | SRSF2 | GAATCCTG       | 4.840000  | 158   | 165 | 21657908      | 21657915    |
| 3  | LDHB | SRSF2 | GGCCCACA       | 4.080000  | 165   | 172 | 21657915      | 21657922    |

## What's Next?

For these tools, I plan to allow handling of other rMATS input formats (i.e. IR, A3SS, A5SS, and MXE tables) and eventually create an R package. In a future RBP tool, I aim to develop functions that plot motif enrichment of alternatively spliced events against their ΔPSI values. This will enable investigations into how different RBPs contribute to the differential motif enrichment of alternatively spliced events across two conditions.
  
